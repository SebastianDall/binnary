from src import data_processing as dp
import numpy as np
import pandas as pd

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