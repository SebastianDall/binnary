import pandas as pd
import numpy as np
from src import data_processing as dp 




def include_contigs(motifs_scored_in_bins, motifs_of_interest, args):
    """
    
    
    """
    
    # Step 1 create bin_motif_from_motifs_scored_in_bins - basis for bin-contig comparison
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_motifs_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        motifs_of_interest,
        args
    )
    
    # Create contig motifs binary
    contig_motif_binary = dp.construct_contig_motif_binary(
        motifs_scored_in_bins, 
        motifs_of_interest, 
        args.mean_methylation_cutoff-0.15, 
        args.n_motif_cutoff
    )
    
    # Merge bin_motifs_from_motifs_scored_in_bins and contig_motif_binary    
    motif_binary_compare = pd.merge(
        bin_motifs_from_motifs_scored_in_bins, contig_motif_binary, on="motif_mod"
    )
    
    
    # match pattern between bin and contig
    # Define the conditions
    conditions = [
        (motif_binary_compare["methylation_binary"] == 1)
        & (motif_binary_compare["methylation_binary_compare"] == 1),
        (motif_binary_compare["methylation_binary"] == 1)
        & (motif_binary_compare["methylation_binary_compare"] == 0),
        (motif_binary_compare["methylation_binary"] == 0)
        & (motif_binary_compare["methylation_binary_compare"] == 1),
        (motif_binary_compare["methylation_binary"] == 0)
        & (motif_binary_compare["methylation_binary_compare"] == 0),
        (motif_binary_compare["methylation_binary"] == 1)
        & (motif_binary_compare["methylation_binary_compare"].isna()),
        (motif_binary_compare["methylation_binary"] == 0)
        & (motif_binary_compare["methylation_binary_compare"].isna()),
    ]

    # Define the corresponding choices for each condition
    choices = [
        0,  # bin motif is methylated, contig motif is methylated
        1,  # bin motif is methylated, contig motif is not methylated
        1,  # bin motif is not methylated, contig motif is methylated
        0,  # bin motif is not methylated, contig motif is not methylated
        1,  # bin motif is methylated, contig motif is not observed
        1,  # bin motif is not methylated, contig motif is not observed
    ]

    # Use numpy.select to apply these conditions and choices
    motif_binary_compare["motif_comparison_score"] = np.select(
        conditions, choices, default=np.nan
    )
    
    # sum motif_comparison_score by bin
    contig_bin_comparison_score = (
        motif_binary_compare.groupby(["bin", "bin_compare"])["motif_comparison_score"]
        .sum()
        .reset_index(name="binary_methylation_missmatch_score")
    )

    
    
    # Split bin_compare into bin and contig
    contig_bin_comparison_score[["contig_bin", "contig", "contig_number", "length"]] = contig_bin_comparison_score["bin_compare"].str.split("_", expand=True)
    contig_bin_comparison_score["contig"] = (contig_bin_comparison_score["contig"] + "_" + contig_bin_comparison_score["contig_number"])
    contig_bin_comparison_score = contig_bin_comparison_score.drop(columns=["contig_number"])
    
    # Find contigs with no methylation
    contigs_w_no_methylation = contig_motif_binary[contig_motif_binary.groupby("bin_compare")["methylation_binary_compare"].transform("sum") == 0]["bin_compare"].unique()
    
    
    # 
    new_contig_bins = contig_bin_comparison_score[contig_bin_comparison_score["bin_compare"].str.contains("unbinned")]
    
    new_contig_bins = new_contig_bins[
        (~new_contig_bins["bin_compare"].isin(contigs_w_no_methylation)) &  # Remove contigs with no methylation
        (new_contig_bins["binary_methylation_missmatch_score"] == 0)        # Retain contigs with no methylation missmatch   
    ].drop_duplicates(subset=["contig"], keep=False)                        # Remove duplicate contigs
    
    return new_contig_bins