import pandas as pd
import numpy as np


def detect_contamination(motifs_scored_in_bins, bin_motifs, args):
    """ """
    # Create find bin motifs binary
    ## Step 1: Group by bin and motif_mod and calculate mean methylation
    bin_motif_binary = (
        bin_motifs.groupby(["bin", "motif_mod"])["mean"]
        .mean()
        .reset_index(name="mean_methylation")
    )

    ## Step 2: Convert mean methylation values to binary
    bin_motif_binary["methylation_binary"] = (
        bin_motif_binary["mean_methylation"] >= args.mean_methylation_cutoff
    ).astype(int)

    # Create contig motifs binary
    ## Filter motifs that are not in bin_motif_binary
    contig_motif_binary = motifs_scored_in_bins[
        motifs_scored_in_bins["motif_mod"].isin(bin_motif_binary["motif_mod"])
    ]

    ## Filter motifs that are not observed more than n_motif_cutoff times
    contig_motif_binary = contig_motif_binary[
        contig_motif_binary["n_motifs"] >= args.n_motif_cutoff
    ]

    ## Convert mean methylation values to binary
    contig_motif_binary["methylation_binary"] = (
        contig_motif_binary["mean"] >= args.mean_methylation_cutoff
    ).astype(int)

    ## Remove unbinned contigs
    contig_motif_binary = contig_motif_binary[contig_motif_binary["bin"] != "unbinned"]

    ## Rename bin_contig to bin
    contig_motif_binary = contig_motif_binary[
        ["bin_contig", "motif_mod", "methylation_binary"]
    ]
    contig_motif_binary.rename(columns={"bin_contig": "bin"}, inplace=True)

    # Pivot the DataFrame
    contig_motif_binary_pivoted = contig_motif_binary.pivot_table(
        index="bin", columns="motif_mod", values="methylation_binary", fill_value=None
    )
    # Unpivot the DataFrame back to long format
    contig_motif_binary = contig_motif_binary_pivoted.reset_index().melt(
        id_vars=["bin"], var_name="motif_mod", value_name="methylation_binary"
    )

    ## Remove contigs with no methylated motifs
    ### Fill na == 0 before
    contig_motif_binary["methylation_binary"] = contig_motif_binary[
        "methylation_binary"
    ].fillna(0)

    contig_motif_binary = contig_motif_binary[
        contig_motif_binary.groupby("bin")["methylation_binary"].transform("sum") > 0
    ]

    # Combine bin_motif_binary and contig_motif_binary
    contig_motif_binary = contig_motif_binary.rename(
        columns={
            "bin": "bin_compare",
            "methylation_binary": "methylation_binary_compare",
        }
    )
    motif_binary_compare = pd.merge(
        bin_motif_binary, contig_motif_binary, on="motif_mod"
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
    contig_bin_comparison_score[
        ["contig_bin", "contig", "contig_number", "length"]
    ] = contig_bin_comparison_score["bin_compare"].str.split("_", expand=True)
    contig_bin_comparison_score["contig"] = (
        contig_bin_comparison_score["contig"]
        + "_"
        + contig_bin_comparison_score["contig_number"]
    )
    contig_bin_comparison_score = contig_bin_comparison_score.drop(
        columns=["contig_number"]
    )

    # Filter contig_bin == bin and contig_bin_comparison_score > 0
    contamination_contigs = contig_bin_comparison_score[
        (
            contig_bin_comparison_score["bin"]
            == contig_bin_comparison_score["contig_bin"]
        )
        & (contig_bin_comparison_score["binary_methylation_missmatch_score"] > 0)
    ]
    print(contamination_contigs)

    # TODO: Find alternative bin for contamination contigs where binary_methylation_missmatch_score != 0

    # Find alternative bin for contamination contigs
    ## Must have a perfect match
    contamination_contigs_alternative_bin = contig_bin_comparison_score[
        (
            contig_bin_comparison_score["bin"]
            != contig_bin_comparison_score["contig_bin"]
        )
        & (contig_bin_comparison_score["binary_methylation_missmatch_score"] == 0)
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

    return contamination_contigs
