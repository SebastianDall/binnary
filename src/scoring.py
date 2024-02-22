from src import data_processing as dp
import numpy as np
import pandas as pd
from multiprocessing import Pool, Queue, Process
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
    motif_binary_compare["methylation_mean_threshold"] = np.where(
        motif_binary_compare["methylation_binary"] == 1,
        np.maximum(motif_binary_compare["mean_methylation"] - 4 * motif_binary_compare["std_methylation_bin"], 0.1),
        np.nan
    )
    
    # Calculate the binary methylation value for each motif in each bin where the bin consensus is 1.
    # If the mean methylation is above the threshold OR the contig mean methylation is above 0.4, the motif is considered methylated in the contig.
    motif_binary_compare["methylation_binary_compare"] = np.where(
        (motif_binary_compare["methylation_binary"] == 1) & 
        ((motif_binary_compare["mean"] >= motif_binary_compare["methylation_mean_threshold"]) | 
        (motif_binary_compare["mean"] > 0.4)),
        1,  # Methylated
        np.where(motif_binary_compare["methylation_binary"] == 1, 0, np.nan)  # Unmethylated or NaN
    )
    
    # Calculate score for bin consensus is 0 
    motif_binary_compare["methylation_mean_threshold"] = np.where(
        motif_binary_compare["methylation_binary"] == 0,
        0.25,
        motif_binary_compare["methylation_mean_threshold"]
    )
    
    motif_binary_compare["methylation_binary_compare"] = np.where(
        motif_binary_compare["methylation_binary"] == 0,
        (motif_binary_compare["mean"] >= 0.25).astype(int),
        motif_binary_compare["methylation_binary_compare"]
    )
    
    return motif_binary_compare


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
    motif_binary_compare = pd.merge(
        bin_motifs_from_motifs_scored_in_bins,
        motifs_scored_in_contigs[motifs_scored_in_contigs["bin_compare"] == bin_contig],
        on="motif_mod"
    )

    # Define methylation thresholds
    motif_binary_compare = define_mean_methylation_thresholds(motif_binary_compare)

    # Calculate the comparison score regardless of methylation presence
    contig_bin_comparison_score = dp.compare_methylation_pattern(motif_binary_compare, choices)
    
    # Check if the contig has no methylation and note it, but do not exclude it from further processing
    contigHasNMethylation = motif_binary_compare["methylation_binary_compare"].sum()
    logger.info(f"Finished processing {bin_contig}. Contig has {contigHasNMethylation} positive methylation comparisons.")
    
    
    if contigHasNMethylation == 0:
        return contig_bin_comparison_score, bin_contig
    
    
    return contig_bin_comparison_score, None

def compare_methylation_pattern_multiprocessed(motifs_scored_in_bins, choices, args, num_processes=None):
    logger = logging.getLogger(__name__)
    logger.info("Starting comparison of methylation patterns")
    
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_motifs_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )

    motifs_scored_in_contigs = motifs_scored_in_bins[motifs_scored_in_bins["n_motifs"] >= args.n_motif_contig_cutoff]
    motifs_scored_in_contigs = motifs_scored_in_contigs[["bin_contig", "motif_mod", "mean"]]
    motifs_scored_in_contigs.rename(columns={"bin_contig": "bin_compare"}, inplace=True)

    comparison_score = pd.DataFrame()
    contigs_w_no_methylation = []

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(
            process_bin_contig,
            [
                (bin_contig, bin_motifs_from_motifs_scored_in_bins, motifs_scored_in_contigs, choices)
                for bin_contig in motifs_scored_in_contigs["bin_compare"].unique()
            ]
        )

    for result, no_methylation in results:
        if result is not None:
            comparison_score = pd.concat([comparison_score, result])
        if no_methylation is not None:
            contigs_w_no_methylation.append(no_methylation)

    return comparison_score, contigs_w_no_methylation




def compare_methylation_pattern(motifs_scored_in_bins, choices, args):
    
    # Step 1 create bin_motif_from_motifs_scored_in_bins - basis for bin-contig comparison
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_motifs_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    
    ## Filter motifs that are not observed more than n_motif_cutoff times
    motifs_scored_in_contigs = motifs_scored_in_bins[motifs_scored_in_bins["n_motifs"] >= args.n_motif_contig_cutoff]   
    
    ## Rename bin_contig to bin
    motifs_scored_in_contigs = motifs_scored_in_contigs[["bin_contig", "motif_mod", "mean"]]
    motifs_scored_in_contigs.rename(columns={"bin_contig": "bin_compare"}, inplace=True)
    
    
    
    ## 
    comparison_score = pd.DataFrame()
    contigs_w_no_methylation = []
    
    for bin_contig in motifs_scored_in_contigs["bin_compare"].unique():
        motif_binary_compare = pd.merge(
            bin_motifs_from_motifs_scored_in_bins,
            motifs_scored_in_contigs[motifs_scored_in_contigs["bin_compare"] == bin_contig],
            on = "motif_mod"
        )
        
        motif_binary_compare = define_mean_methylation_thresholds(motif_binary_compare)
        
        # If contig has no methylation, skip
        if motif_binary_compare["methylation_binary_compare"].sum() == 0:
            contigs_w_no_methylation.append(bin_contig)
        
        
        
        contig_bin_comparison_score = dp.compare_methylation_pattern(motif_binary_compare, choices)
        
        comparison_score = pd.concat([comparison_score, contig_bin_comparison_score])
    
        
    return comparison_score, contigs_w_no_methylation
        
        
        
        
        
    