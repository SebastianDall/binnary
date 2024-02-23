from src import data_processing as dp
import numpy as np
import polars as pl
from multiprocessing import Pool, Queue, Process, get_context
import logging
from logging.handlers import QueueHandler, QueueListener

# Set up logging
## Create queue for logging
log_queue = Queue()

## Set up a listener to handle logs from the queue
def setup_logging_queue(queue):
    while True:
        record = queue.get()
        if record is None:  # Use None as a sentinel to stop the listener
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)

def worker_setup_logging(queue):
    q_handler = QueueHandler(queue)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(q_handler)

def define_mean_methylation_thresholds(motif_binary_compare):
    """
    Define mean methylation thresholds for bin and contig motifs
    """
    # Calculate the mean methylation value for each motif in each bin
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when(pl.col("methylation_binary") == 1)
        .then(
            pl.when(pl.col("mean_methylation") - 4 * pl.col("std_methylation_bin") > 0.1)
            .then(pl.col("mean_methylation") - 4 * pl.col("std_methylation_bin"))
            .otherwise(0.1)
        )
        .otherwise(pl.lit(None))
        .alias("methylation_mean_threshold")
    ])

    # Calculate the binary methylation value for each motif in each bin where the bin consensus is 1
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when((pl.col("methylation_binary") == 1) & 
                ((pl.col("mean") >= pl.col("methylation_mean_threshold")) | 
                (pl.col("mean") > 0.4)))
        .then(1)
        .when(pl.col("methylation_binary") == 1)
        .then(0)
        .otherwise(pl.lit(None))
        .alias("methylation_binary_compare")
    ])

    # Calculate score for bin consensus is 0
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when(pl.col("methylation_binary") == 0)
        .then(0.25)
        .otherwise(pl.col("methylation_mean_threshold"))
        .alias("methylation_mean_threshold"),

        pl.when(pl.col("methylation_binary") == 0)
        .then((pl.col("mean") >= 0.25).cast(pl.Int32))
        .otherwise(pl.col("methylation_binary_compare"))
        .alias("methylation_binary_compare")
    ])
    
    return motif_binary_compare



def compare_methylation_pattern(motif_binary_compare, choices):
    """
    Compares the methylation pattern between bin and contig motifs using Polars and calculates the motif_comparison_score.
    """
    motif_comparison_score = (
        pl.when((motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 1))
        .then(0)
        .when((motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 0))
        .then(1)
        .when((motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 1))
        .then(1)
        .when((motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 0))
        .then(0)
        .when((motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'].is_null()))
        .then(0)
        .when((motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'].is_null()))
        .then(0)
        .otherwise(pl.lit(None))
    )

    # Add the 'motif_comparison_score' column to the DataFrame
    motif_binary_compare = motif_binary_compare.with_columns(motif_comparison_score.alias("motif_comparison_score"))
    
    # Group by bin and bin_compare and calculate the sum of the motif_comparison_score and the count of non-NA values
    # contig_bin_comparison_score = motif_binary_compare \
    #     .group_by(["bin", "bin_compare"]) \
    #     .agg(
    #         pl.sum("motif_comparison_score").alias("binary_methylation_missmatch_score"),
    #         pl.count("motif_comparison_score").alias("non_na_comparisons")
    #     ) \
    #     .select(
    #         pl.col('bin_compare').str.split('_').arr.to_struct(n_field_strategy="max_width")
    #     ).unnest('bin_compare') \
    #     .drop("contig_number")
        
    contig_bin_comparison_score = motif_binary_compare \
        .group_by(["bin", "bin_compare"]) \
        .agg([
            pl.sum("motif_comparison_score").alias("binary_methylation_missmatch_score"),
            pl.count("motif_comparison_score").alias("non_na_comparisons")
        ])
    # contig_bin_comparison_score = motif_binary_compare.groupby(["bin", "bin_compare"]).agg(
    #     binary_methylation_missmatch_score=pd.NamedAgg(column="motif_comparison_score", aggfunc="sum"),
    #     non_na_comparisons=pd.NamedAgg(column="motif_comparison_score", aggfunc="count")  # Count non-NA values
    # ).reset_index()
    
    # Split bin_compare into bin and contig
    # contig_bin_comparison_score[["contig_bin", "contig", "contig_number"]] = contig_bin_comparison_score["bin_compare"].str.split("_", expand=True)
    # contig_bin_comparison_score["contig"] = (contig_bin_comparison_score["contig"] + "_" + contig_bin_comparison_score["contig_number"])
    # contig_bin_comparison_score = contig_bin_comparison_score.drop(columns=["contig_number"])
    
    
    return contig_bin_comparison_score


# def process_bin_contig(bin_contig, bin_motifs_from_motifs_scored_in_bins, motifs_scored_in_contigs, choices):
#     motif_binary_compare = pd.merge(
#         bin_motifs_from_motifs_scored_in_bins,
#         motifs_scored_in_contigs[motifs_scored_in_contigs["bin_compare"] == bin_contig],
#         on="motif_mod"
#     )

#     motif_binary_compare = define_mean_methylation_thresholds(motif_binary_compare)

#     if motif_binary_compare["methylation_binary_compare"].sum() == 0:
#         return None, bin_contig

#     contig_bin_comparison_score = dp.compare_methylation_pattern(motif_binary_compare, choices)
#     return contig_bin_comparison_score, None


def process_bin_contig(bin_contig, bin_motifs_from_motifs_scored_in_bins, motifs_scored_in_contigs, choices):
    worker_setup_logging(log_queue)
    logger = logging.getLogger(__name__)
    logger.info(f"Processing {bin_contig}")
    
    
    # Merge motifs based on 'motif_mod'
    motif_binary_compare = bin_motifs_from_motifs_scored_in_bins.join(
        motifs_scored_in_contigs.filter(pl.col("bin_compare") == bin_contig),
        on="motif_mod"
    )
    # motif_binary_compare = pd.merge(
    #     bin_motifs_from_motifs_scored_in_bins,
    #     motifs_scored_in_contigs[motifs_scored_in_contigs["bin_compare"] == bin_contig],
    #     on="motif_mod"
    # )

    # Define methylation thresholds
    motif_binary_compare = define_mean_methylation_thresholds(motif_binary_compare)

    # Calculate the comparison score regardless of methylation presence
    contig_bin_comparison_score = compare_methylation_pattern(motif_binary_compare, choices)
    
    # Check if the contig has no methylation and note it, but do not exclude it from further processing
    # contigHasNMethylation = motif_binary_compare["methylation_binary_compare"].sum()
    contigHasNMethylation = motif_binary_compare.filter(pl.col("methylation_binary_compare") == 1).height
    logger.info(f"Finished processing {bin_contig}. Contig has {contigHasNMethylation} positive methylation comparisons.")
    
    
    if contigHasNMethylation == 0:
        return contig_bin_comparison_score, bin_contig
    
    
    return contig_bin_comparison_score, None

def compare_methylation_pattern_multiprocessed(motifs_scored_in_bins, bin_consensus, choices, args, num_processes=1):
    logger = logging.getLogger(__name__)
    logger.info("Starting comparison of methylation patterns")
    
    motifs_scored_in_contigs = motifs_scored_in_bins \
        .filter(pl.col("n_motifs") >= args.n_motif_contig_cutoff) \
        .select(["bin_contig", "motif_mod", "mean"]) \
        .rename({"bin_contig": "bin_compare"})
    

    comparison_score = pl.DataFrame()
    contigs_w_no_methylation = []

    with get_context("spawn").Pool(processes=num_processes) as pool:
        results = pool.starmap(
            process_bin_contig,
            [
                (bin_contig, bin_consensus, motifs_scored_in_contigs, choices)
                for bin_contig in motifs_scored_in_contigs.select("bin_compare").unique().to_pandas()["bin_compare"].tolist()
            ]
        )
        
    for result, no_methylation in results:
        if result is not None:
            comparison_score = pl.concat([comparison_score, result])
        if no_methylation is not None:
            contigs_w_no_methylation.append(no_methylation)
        
    # for loop implementation
    # for bin_contig in motifs_scored_in_contigs.select("bin_compare").unique().to_pandas()["bin_compare"].tolist():
    #     result, no_methylation = process_bin_contig(bin_contig, bin_consensus, motifs_scored_in_contigs, choices)
    #     if result is not None:
    #         comparison_score = pl.concat([comparison_score, result])
    #     if no_methylation is not None:
    #         contigs_w_no_methylation.append(no_methylation)
        
    



    return comparison_score, contigs_w_no_methylation




# def compare_methylation_pattern(motifs_scored_in_bins, choices, args):
    
#     # Step 1 create bin_motif_from_motifs_scored_in_bins - basis for bin-contig comparison
#     bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_motifs_from_motifs_scored_in_bins(
#         motifs_scored_in_bins,
#         args
#     )
    
    
#     ## Filter motifs that are not observed more than n_motif_cutoff times
#     motifs_scored_in_contigs = motifs_scored_in_bins[motifs_scored_in_bins["n_motifs"] >= args.n_motif_contig_cutoff]   
    
#     ## Rename bin_contig to bin
#     motifs_scored_in_contigs = motifs_scored_in_contigs[["bin_contig", "motif_mod", "mean"]]
#     motifs_scored_in_contigs.rename(columns={"bin_contig": "bin_compare"}, inplace=True)
    
    
    
#     ## 
#     comparison_score = pd.DataFrame()
#     contigs_w_no_methylation = []
    
#     for bin_contig in motifs_scored_in_contigs["bin_compare"].unique():
#         motif_binary_compare = pd.merge(
#             bin_motifs_from_motifs_scored_in_bins,
#             motifs_scored_in_contigs[motifs_scored_in_contigs["bin_compare"] == bin_contig],
#             on = "motif_mod"
#         )
        
#         motif_binary_compare = define_mean_methylation_thresholds(motif_binary_compare)
        
#         # If contig has no methylation, skip
#         if motif_binary_compare["methylation_binary_compare"].sum() == 0:
#             contigs_w_no_methylation.append(bin_contig)
        
        
        
#         contig_bin_comparison_score = dp.compare_methylation_pattern(motif_binary_compare, choices)
        
#         comparison_score = pd.concat([comparison_score, contig_bin_comparison_score])
    
        
#     return comparison_score, contigs_w_no_methylation
        
        
        
        
        
    