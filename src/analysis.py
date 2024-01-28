import pandas as pd
import numpy as np
from .data_processing import prepare_motifs_scored_in_bins


def perform_analysis(motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file, args):
    """
    Performs the core analysis of the tool. This function is called from the main entry point.
    """
    # Find relevant motifs in bins and contigs
    motifs_scored_in_bins = prepare_motifs_scored_in_bins(motifs_scored, bin_motifs, contig_bins, assembly_stats)
    
    # Calculate belonging score
    belonging_score = calculate_belonging_score(motifs_scored_in_bins, args)
    
    
    return belonging_score





def calculate_belonging_score(motifs_scored_in_bins, args):
    """
    Finds the closest match between contig and bin based on methylation pattern.
    """
    # Step 1: Group by bin and motif_mod and calculate mean methylation
    bin_motif_binary = motifs_scored_in_bins[motifs_scored_in_bins['bin'] != "unbinned"].groupby(['bin', 'motif_mod'])['mean'].mean().reset_index(name='mean_methylation')
    
    # Step 2: Convert mean methylation values to binary
    bin_motif_binary['methylation_binary'] = (bin_motif_binary['mean_methylation'] >= args.mean_methylation_cutoff).astype(int)

    
    # Create contig_motif_binary
    ## Filter motifs that are not observed in bins
    contig_motif_binary = motifs_scored_in_bins[motifs_scored_in_bins["motif_mod"].isin(bin_motif_binary["motif_mod"])]
    
    ## Filter motifs that are not observed more than n_motif_cutoff times
    contig_motif_binary = contig_motif_binary[contig_motif_binary["n_motifs"] >= args.n_motif_cutoff]
    
    ## Convert mean methylation values to binary
    contig_motif_binary["methylation_binary"] = (contig_motif_binary["mean"] >= args.mean_methylation_cutoff).astype(int)
    
    ## Rename bin_contig to bin
    contig_motif_binary = contig_motif_binary[["bin_contig", "motif_mod", "methylation_binary"]]
    contig_motif_binary.rename(columns={"bin_contig": "bin"}, inplace=True)
    
    # Pivot the DataFrame
    contig_motif_binary_pivoted = contig_motif_binary.pivot_table(
        index='bin', columns='motif_mod', values='methylation_binary', fill_value=None
    )
    # Unpivot the DataFrame back to long format
    contig_motif_binary = contig_motif_binary_pivoted.reset_index().melt(
        id_vars=['bin'], var_name='motif_mod', value_name='methylation_binary'
    )
    
    # Check if there is a bin with no methylated motifs
    bin_with_no_methylations_exists = bin_motif_binary.groupby("bin")["methylation_binary"].sum().min() == 0
    
    if bin_with_no_methylations_exists == False:
        print("No bin with no methylated motifs exists. Removing contigs with no methylated motifs...")
        # Remove contigs with no methylated motifs
        contig_motif_binary = contig_motif_binary[contig_motif_binary.groupby("bin")["methylation_binary"].transform("sum") > 0]

    
    # Combine bin_motif_binary and contig_motif_binary
    contig_motif_binary = contig_motif_binary.rename(columns={
        'bin': 'bin_compare',
        'methylation_binary': 'methylation_binary_compare'
    })
    
    motif_binary_compare = pd.merge(bin_motif_binary, contig_motif_binary, on='motif_mod')
    
    # Define the conditions
    conditions = [
        (motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 1),
        (motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 0),
        (motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 1),
        (motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 0),
        (motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'].isna()),
        (motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'].isna())
    ]

    # Define the corresponding choices for each condition
    choices = [1, -1, 0, 0, 0, 0]

    # Use numpy.select to apply these conditions and choices
    motif_binary_compare['motif_comparison_score'] = np.select(conditions, choices, default=np.nan)
    
    
    # Calculate belonging score
    # ## Sum motif_comparison_score by bin and bin_compare
    # belonging_score = motif_binary_compare.groupby(["bin","bin_compare"])["motif_comparison_score"].sum().reset_index(name="belonging_score")
    # ## Retain the highest belonging score for each bin_compare (contig)
    # belonging_score = belonging_score.groupby("bin_compare")["belonging_score"].max().reset_index(name="belonging_score")
    # print(belonging_score)
    
    # ## Add belonging_bins column
    # belonging_score["belonging_bins"] = belonging_score.groupby("bin_compare")["bin_compare"].transform("count")
    
    # Calculate the sum of motif_comparison_score for each combination of bin and bin_compare
    belonging_score = motif_binary_compare.groupby(["bin", "bin_compare"])["motif_comparison_score"].sum().reset_index(name="belonging_score")
    
    # Find the max belonging_score for each bin_compare
    max_scores = belonging_score.groupby("bin_compare")["belonging_score"].max().reset_index(name="max_belonging_score")

    # Merge the max_scores back with the original belonging_score DataFrame to get the corresponding bin values
    # This merge operation ensures that we only retain the rows with the max belonging_score for each bin_compare
    belonging_score_with_max = pd.merge(belonging_score, max_scores, how='inner', left_on=["bin_compare", "belonging_score"], right_on=["bin_compare", "max_belonging_score"])

    # Count the number of best matches for each bin_compare
    belonging_score_with_max["belonging_bins"] = belonging_score_with_max.groupby("bin_compare")["bin_compare"].transform("count")
    
    
    # Drop the max_belonging_score column as it's redundant now
    belonging_score_final = belonging_score_with_max.drop(columns=["max_belonging_score"])

    return belonging_score_final


