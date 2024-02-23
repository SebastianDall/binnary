import pandas as pd
import polars as pl
import numpy as np
from src import data_processing as dp 
from src import scoring as sc
from src import utils as ut
import logging



def include_contigs(motifs_scored_in_bins, bin_consensus, contamination, args):
    """
    Takes the motifs_scored_in_bins and motifs_of_interest DataFrames and finds unbinned contigs with an exact methylation pattern as the bin.
    
    params:
        motifs_scored_in_bins: pd.DataFrame - DataFrame containing the motifs scored in each bin
        motifs_of_interest: list - List of motifs to be considered from bin_motif_binary
        args: argparse.Namespace - Namespace containing the arguments passed to the script
    
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting include_contigs analysis...")
    
    
    # Remove bins with no methylation in consensus
    bins_w_no_methylation = bin_consensus \
        .group_by("bin") \
        .agg(
            pl.sum("methylation_binary").alias("binary_sum")    
        ) \
        .filter(pl.col("binary_sum") == 0) \
        .select("bin") \
        .unique()["bin"]
    
    
    # bins_w_no_methylation = bin_consensus[
    #     bin_consensus.groupby("bin")["methylation_binary"].transform("sum") == 0
    # ]["bin"].unique()
    
    bin_consensus = bin_consensus \
        .filter(~pl.col("bin").is_in(bins_w_no_methylation))
    
    # bin_consensus = bin_consensus[
    #     ~bin_consensus["bin"].isin(bins_w_no_methylation)
    # ]
    
    # Retain only unbinned contigs or contigs in the contamination file
    contigs_for_comparison = motifs_scored_in_bins \
        .filter(
            (pl.col("bin_contig").str.contains("unbinned")) |
            (pl.col("bin_contig").is_in(contamination["bin_contig_compare"]))
        )
    
    # contigs_for_comparison = motifs_scored_in_bins[
    #     (motifs_scored_in_bins["bin_contig"].str.contains("unbinned")) |  # Retain unbinned contigs
    #     (motifs_scored_in_bins["bin_contig"].isin(contamination["bin_contig_compare"])) # Retain contigs in the contamination file
    # ]
    
    # Define the corresponding choices for each condition
    choices = [
        0,  # bin motif is methylated, contig motif is methylated
        1,  # bin motif is methylated, contig motif is not methylated
        1,  # bin motif is not methylated, contig motif is methylated
        0,  # bin motif is not methylated, contig motif is not methylated
        0,  # bin motif is methylated, contig motif is not observed
        0,  # bin motif is not methylated, contig motif is not observed
    ]

    contig_bin_comparison_score, contigs_w_no_methylation = sc.compare_methylation_pattern_multiprocessed(
        contigs_for_comparison,
        bin_consensus,
        choices,
        args,
        num_processes=args.threads
    )
    
    contig_bin_comparison_score = ut.split_bin_contig(contig_bin_comparison_score)
    
    dp.generate_output(contig_bin_comparison_score.to_pandas(), args.out, "motif_binary_comparison.tsv")
    
    # Filter contigs where motif comparisons are less than args.min_motif_comparisons
    contigs_of_interest = contig_bin_comparison_score \
        .filter(
            pl.col("non_na_comparisons") >= args.min_motif_comparisons,
            (~pl.col("bin_compare").is_in(contigs_w_no_methylation)) &  # Remove contigs with no methylation
            (pl.col("binary_methylation_missmatch_score") == 0)#,        # Retain contigs with no methylation missmatch
            # pl.col("contig").is_unique()
        ) \
        .sort("bin","bin_compare")
        
    
    
    # contig_bin_comparison_score = contig_bin_comparison_score[
    #     contig_bin_comparison_score["non_na_comparisons"] >= args.min_motif_comparisons
    # ]
    
    
    logger.info("Assigning contigs to bins...")
    
    # contigs_of_interest = contig_bin_comparison_score[
    #     (~contig_bin_comparison_score["bin_compare"].isin(contigs_w_no_methylation)) &  # Remove contigs with no methylation
    #     (contig_bin_comparison_score["binary_methylation_missmatch_score"] == 0)        # Retain contigs with no methylation missmatch   
    # ].drop_duplicates(subset=["contig"], keep=False)                        # Remove duplicate contigs
    
    # sort by bin
    # contigs_of_interest = contigs_of_interest.sort_values(by=["bin", "bin_compare"])
    logger.info("Finished include_contigs analysis.")
    return contigs_of_interest