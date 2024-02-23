import polars as pl
import pandas as pd
from src import data_processing as dp
from src import scoring as sc
import logging

def detect_contamination(motifs_scored_in_bins, bin_consensus, args):
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
    logger = logging.getLogger(__name__)
    logger.info("Starting contamination detection analysis...")
    motifs_scored_in_bins_wo_unbinned = motifs_scored_in_bins \
        .filter(~pl.col("bin_contig").str.contains("unbinned"))
    
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
        motifs_scored_in_bins_wo_unbinned,
        bin_consensus,
        choices,
        args,
        num_processes=args.threads
    )

    logger.info("Finding contamination in bins")
    
    bin_compare_column = contig_bin_comparison_score.select("bin_compare").to_pandas()
    bin_compare_column[["contig_bin", "contig", "contig_number"]] = bin_compare_column["bin_compare"].str.split("_", expand=True)
    bin_compare_column = bin_compare_column.drop(columns=["contig"]).drop_duplicates()
    bin_compare_column_pl = pl.DataFrame(bin_compare_column)
    
    contig_bin_comparison_score = contig_bin_comparison_score \
        .join(bin_compare_column_pl, on="bin_compare", how="left") \
        .with_columns(("contig_" + pl.col("contig_number")).alias("contig")) \
        .drop("contig_number")
    # 
    
    # split bin_compare column
    # print(contig_bin_comparison_score.columns)
    
    # contig_bin_comparison_score[["contig_bin", "contig", "contig_number"]] = contig_bin_comparison_score["bin_compare"].str.split("_", expand=True)
    # contig_bin_comparison_score["contig"] = (contig_bin_comparison_score["contig"] + "_" + contig_bin_comparison_score["contig_number"])
    # contig_bin_comparison_score = contig_bin_comparison_score.drop(columns=["contig_number"])
    
    # Filter contig_bin == bin and contig_bin_comparison_score > 0
    contamination_contigs = contig_bin_comparison_score \
        .filter(
            (pl.col("bin") == pl.col("contig_bin")) &
            (pl.col("binary_methylation_missmatch_score") > 0)
        )
    
    
    # contamination_contigs = contig_bin_comparison_score[
    #     # NOTE: This line also removes all contigs from bins with no methylation
    #     (contig_bin_comparison_score["bin"] == contig_bin_comparison_score["contig_bin"]) &
    #     (contig_bin_comparison_score["binary_methylation_missmatch_score"] > 0)
    # ]

    logger.info("Finding alternative bin for contamination contigs")
    # Find alternative bin for contamination contigs
    ## Must have a perfect match
    contamination_contigs_alternative_bin = contig_bin_comparison_score \
        .filter(
            (pl.col("bin") != pl.col("contig_bin")) &
            (pl.col("binary_methylation_missmatch_score") == 0) &
            (~pl.col("bin_compare").is_in(contigs_w_no_methylation))
        ) \
        .select(["contig", "bin", "binary_methylation_missmatch_score"]) \
        .rename(
            {
                "bin": "alternative_bin",
                "binary_methylation_missmatch_score": "alternative_bin_binary_methylation_missmatch_score"
            }
        )
    
    
    
    # contamination_contigs_alternative_bin = contig_bin_comparison_score[
    #     # This line removes all bin - contig mathces where the bin is the same as the contig
    #     (contig_bin_comparison_score["bin"] != contig_bin_comparison_score["contig_bin"]) &
    #     # This line has a side consequence that all contigs from bins with no methylation are removed
    #     (contig_bin_comparison_score["binary_methylation_missmatch_score"] == 0) & 
    #     (~contig_bin_comparison_score["bin_compare"].isin(contigs_w_no_methylation))
    # ]
    # contamination_contigs_alternative_bin = contamination_contigs_alternative_bin[
    #     ["contig", "bin", "binary_methylation_missmatch_score"]
    # ].rename(
    #     columns={
    #         "bin": "alternative_bin",
    #         "binary_methylation_missmatch_score": "alternative_bin_binary_methylation_missmatch_score",
    #     }
    # )

    contamination_contigs = contamination_contigs \
        .join(contamination_contigs_alternative_bin, on="contig", how="left") \
        .drop("contig_bin") \
            .rename(
                {
                    "bin_compare": "bin_contig_compare"
                }
            ) \
        .sort("bin", "bin_contig_compare")

    # contamination_contigs = pd.merge(
    #     contamination_contigs,
    #     contamination_contigs_alternative_bin,
    #     on="contig",
    #     how="left",
    # )
    
    # Remove redundant columns
    # contamination_contigs = contamination_contigs.drop(columns=["contig_bin"])
    # # Rename bin_compare
    # contamination_contigs = contamination_contigs.rename(
    #     columns={"bin_compare": "bin_contig_compare"}
    # )

    # # sort by bin
    # contamination_contigs = contamination_contigs.sort_values(by=["bin", "bin_contig_compare"])
    
    logger.info("Contamination detection complete")
    
    return contamination_contigs
