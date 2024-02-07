import pandas as pd
import numpy as np
from src import data_processing as dp 

def detect_contamination(motifs_scored_in_bins, motifs_of_interest, args):
    """
    Takes the bin_motif_binary and motifs_scored_in_bins DataFrames and performs the contamination detection analysis.
    Firstly bin_motif_binary is used to create a binary representation of the methylation status of each motif in each bin.
    This is the bases for the comparison. Only motifs that has a mean methylation value of at least 0.75 in bin_motif_binary are considered.
    
    Motifs_scored_in_bins is then used to create a binary representation of the methylation status of each motif in each contig.
    This is then compared to the binary representation of the bin to identify any mismatches.
    Only motifs that are observed more than 6 times are considered and must have a mean methylation value of at least 0.75-0.15.
    
    
    params:
        motifs_scored_in_bins: pd.DataFrame - DataFrame containing the motifs scored in each bin
        motifs_of_interest: list - List of motifs to be considered from bin_motif_binary
        args: argparse.Namespace - Namespace containing the arguments passed to the script
    
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
    # Remove unbinned motifs
    contig_motif_binary = contig_motif_binary[~contig_motif_binary["bin_compare"].str.contains("unbinned")]

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
        0,  # bin motif is methylated, contig motif is not observed
        0,  # bin motif is not methylated, contig motif is not observed
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
    
    
    # Filter contig_bin == bin and contig_bin_comparison_score > 0
    contamination_contigs = contig_bin_comparison_score[
        # NOTE: This line also removes all contigs from bins with no methylation
        (contig_bin_comparison_score["bin"] == contig_bin_comparison_score["contig_bin"]) &
        (contig_bin_comparison_score["binary_methylation_missmatch_score"] > 0)
    ]
    
    contigs_w_no_methylation = contig_motif_binary[contig_motif_binary.groupby("bin_compare")["methylation_binary_compare"].transform("sum") == 0]["bin_compare"].unique()
    
    # TODO: Find alternative bin for contamination contigs where binary_methylation_missmatch_score != 0

    # Find alternative bin for contamination contigs
    ## Must have a perfect match
    contamination_contigs_alternative_bin = contig_bin_comparison_score[
        # This line removes all bin - contig mathces where the bin is the same as the contig
        (contig_bin_comparison_score["bin"] != contig_bin_comparison_score["contig_bin"]) &
        # This line has a side consequence that all contigs from bins with no methylation are removed
        (contig_bin_comparison_score["binary_methylation_missmatch_score"] == 0) & 
        (~contig_bin_comparison_score["bin_compare"].isin(contigs_w_no_methylation))
    ]
    contamination_contigs_alternative_bin = contamination_contigs_alternative_bin[
        ["contig", "bin", "binary_methylation_missmatch_score"]
    ].rename(
        columns={
            "bin": "alternative_bin",
            "binary_methylation_missmatch_score": "alternative_bin_binary_methylation_missmatch_score",
        }
    )

    contamination_contigs = pd.merge(
        contamination_contigs,
        contamination_contigs_alternative_bin,
        on="contig",
        how="left",
    )
    
    # Remove redundant columns
    contamination_contigs = contamination_contigs.drop(columns=["contig_bin"])
    # Rename bin_compare
    contamination_contigs = contamination_contigs.rename(
        columns={"bin_compare": "bin_contig_compare"}
    )

    return contamination_contigs
