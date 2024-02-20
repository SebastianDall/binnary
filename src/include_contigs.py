import pandas as pd
import numpy as np
from src import data_processing as dp 




def include_contigs(motifs_scored_in_bins, contamination, args):
    """
    Takes the motifs_scored_in_bins and motifs_of_interest DataFrames and finds unbinned contigs with an exact methylation pattern as the bin.
    
    params:
        motifs_scored_in_bins: pd.DataFrame - DataFrame containing the motifs scored in each bin
        motifs_of_interest: list - List of motifs to be considered from bin_motif_binary
        args: argparse.Namespace - Namespace containing the arguments passed to the script
    
    """
    
    motif_binary_compare = dp.calculate_binary_motif_comparison_matrix(
        motifs_scored_in_bins,
        args
    )
    
    dp.generate_output(motif_binary_compare, args.out, "motif_binary_comparison.tsv")
    
    # Remove bins with no methylation in consensus
    bins_w_no_methylation = motif_binary_compare[
        motif_binary_compare.groupby("bin")["methylation_binary"].transform("sum") == 0
    ]["bin"].unique()
    
    # Remove contigs with no methylation from the comparison
    motif_binary_compare = motif_binary_compare[
        ~motif_binary_compare["bin"].isin(bins_w_no_methylation)
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

    contig_bin_comparison_score = dp.compare_methylation_pattern(motif_binary_compare, choices)
    
    # Find contigs with no methylation
    contigs_w_no_methylation = motif_binary_compare[
        motif_binary_compare.groupby("bin_compare")["methylation_binary_compare"].transform("sum") == 0
    ]["bin_compare"].unique()
    
    
    contigs_of_interest = contig_bin_comparison_score[
        (contig_bin_comparison_score["bin_compare"].str.contains("unbinned")) |  # Retain unbinned contigs
        (contig_bin_comparison_score["bin_compare"].isin(contamination["bin_contig_compare"])) # Retain contigs in the contamination file   
    ]
    
    
    contigs_of_interest = contigs_of_interest[
        (~contigs_of_interest["bin_compare"].isin(contigs_w_no_methylation)) &  # Remove contigs with no methylation
        (contigs_of_interest["binary_methylation_missmatch_score"] == 0)        # Retain contigs with no methylation missmatch   
    ].drop_duplicates(subset=["contig"], keep=False)                        # Remove duplicate contigs
    
        
    return contigs_of_interest