import pandas as pd
import numpy as np
from src import data_processing as dp 
from src import scoring as sc
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
    bins_w_no_methylation = bin_consensus[
        bin_consensus.groupby("bin")["methylation_binary"].transform("sum") == 0
    ]["bin"].unique()
    
    bin_consensus = bin_consensus[
        ~bin_consensus["bin"].isin(bins_w_no_methylation)
    ]
    
    # Retain only unbinned contigs or contigs in the contamination file
    contigs_for_comparison = motifs_scored_in_bins[
        (motifs_scored_in_bins["bin_contig"].str.contains("unbinned")) |  # Retain unbinned contigs
        (motifs_scored_in_bins["bin_contig"].isin(contamination["bin_contig_compare"])) # Retain contigs in the contamination file
    ]
    
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
    
    
    dp.generate_output(contig_bin_comparison_score, args.out, "motif_binary_comparison.tsv")
    
    # # Remove bins with no methylation in consensus
    # bins_w_no_methylation = motif_binary_compare[
    #     motif_binary_compare.groupby("bin")["methylation_binary"].transform("sum") == 0
    # ]["bin"].unique()
    
    # # Remove bins with no methylation from the comparison
    # motif_binary_compare = motif_binary_compare[
    #     ~motif_binary_compare["bin"].isin(bins_w_no_methylation)
    # ]
    
    # Remove comparisons between bins and contigs with less than args.min_motif_comparisons from the comparison
    # motif_binary_compare = motif_binary_compare.groupby(['bin', 'bin_compare']).filter(lambda x: x['mean'].count() >= args.min_motif_comparisons)
    
        
    # Define the corresponding choices for each condition
    # choices = [
    #     0,  # bin motif is methylated, contig motif is methylated
    #     1,  # bin motif is methylated, contig motif is not methylated
    #     1,  # bin motif is not methylated, contig motif is methylated
    #     0,  # bin motif is not methylated, contig motif is not methylated
    #     0,  # bin motif is methylated, contig motif is not observed
    #     0,  # bin motif is not methylated, contig motif is not observed
    # ]

    # contig_bin_comparison_score = dp.compare_methylation_pattern(motif_binary_compare, choices)
    
    # # Find contigs with no methylation
    # contigs_w_no_methylation = motif_binary_compare[
    #     motif_binary_compare.groupby("bin_compare")["methylation_binary_compare"].transform("sum") == 0
    # ]["bin_compare"].unique()
    
    # Filter contigs where motif comparisons are less than args.min_motif_comparisons
    contig_bin_comparison_score = contig_bin_comparison_score[
        contig_bin_comparison_score["non_na_comparisons"] >= args.min_motif_comparisons
    ]
    
    
    # contigs_of_interest = contig_bin_comparison_score[
    #     (contig_bin_comparison_score["bin_compare"].str.contains("unbinned")) |  # Retain unbinned contigs
    #     (contig_bin_comparison_score["bin_compare"].isin(contamination["bin_contig_compare"])) # Retain contigs in the contamination file   
    # ]
    
    logger.info("Assigning contigs to bins...")
    
    contigs_of_interest = contig_bin_comparison_score[
        (~contig_bin_comparison_score["bin_compare"].isin(contigs_w_no_methylation)) &  # Remove contigs with no methylation
        (contig_bin_comparison_score["binary_methylation_missmatch_score"] == 0)        # Retain contigs with no methylation missmatch   
    ].drop_duplicates(subset=["contig"], keep=False)                        # Remove duplicate contigs
    
    # sort by bin
    contigs_of_interest = contigs_of_interest.sort_values(by=["bin", "bin_compare"])
    logger.info("Finished include_contigs analysis.")
    return contigs_of_interest