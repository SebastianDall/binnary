import pandas as pd
from Bio import SeqIO
import numpy as np
import os

def load_data(args):
    """
    Load and preprocess data for the analysis.
    """
    # Load the data from the provided files
    motifs_scored = pd.read_csv(args.motifs_scored, delimiter = "\t")
    bin_motifs = pd.read_csv(args.bin_motifs, delimiter = "\t")
    contig_bins = pd.read_csv(args.contig_bins, delimiter = "\t", header = None)
    assembly_stats = pd.read_csv(args.assembly_stats, delimiter = "\t")
    
    # Perform any additional preprocessing steps here
    # Change colnames of contig_bins
    contig_bins.columns = ["contig", "bin"]
    
    # Rename the first column to 'contig'
    assembly_stats.rename(columns={assembly_stats.columns[0]: "contig"}, inplace=True)
    
    return motifs_scored, bin_motifs, contig_bins, assembly_stats



def generate_output(output_df, output_path):
    """
    Generate the output files for the analysis.
    """
    # If the output directory does not exist, create it
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir) and output_dir != "":
        os.makedirs(output_dir)
        
    # Generate the output files
    output_df.to_csv(output_path, sep="\t", index=False)



def prepare_bin_motifs_binary(bin_motifs, args):
    """
    Prepares the bin_motifs_binary DataFrame by calculating the mean methylation per bin and motif_mod and converting it to binary.    
    """
    # Combine 'motif' and 'mod_type' into 'motif_mod'
    bin_motifs["motif_mod"] = bin_motifs["motif"] + "_" + bin_motifs["mod_type"]

    # Calculate total motifs and mean methylation
    bin_motifs["n_motifs"] = bin_motifs["n_mod_bin"] + bin_motifs["n_nomod_bin"]
    bin_motifs["mean_methylation"] = bin_motifs["n_mod_bin"] / bin_motifs["n_motifs"]

    # Create 'methylation_binary' column based on mean methylation cutoff
    bin_motifs["methylation_binary"] = (bin_motifs["mean_methylation"] >= args.mean_methylation_cutoff).astype(int)

    # filter motifs that are not observed more than n_motif_bin_cutoff times
    bin_motifs = bin_motifs[bin_motifs["n_motifs"] >= args.n_motif_bin_cutoff]
    
    # Filter for rows where 'methylation_binary' is 1 and select relevant columns
    bin_motif_binary = bin_motifs[bin_motifs["methylation_binary"] == 1][["bin", "motif_mod", "mean_methylation", "methylation_binary"]]
    
    ## Step 4: Create a binary index for bin and motif_mod
    bin_motif_binary_index = pd.MultiIndex.from_product([bin_motif_binary["bin"].unique(), bin_motif_binary["motif_mod"].unique()], names=['bin', 'motif_mod'])

    bin_motif_binary = bin_motif_binary.set_index(["bin", "motif_mod"])

    bin_motif_binary = bin_motif_binary.reindex(bin_motif_binary_index, fill_value=0).reset_index()
    
    return bin_motif_binary


def prepare_motifs_scored_in_bins(motifs_scored, motifs_of_interest, contig_bins, assembly_stats):
    """
    Prepares the motifs_scored_in_bins DataFrame by merging with bin motifs, contig bins, and assembly stats,
    and calculates additional metrics like number of motifs and mean methylation per contig.
    """
       
    # Filter and enhance motifs_scored based on motifs_in_bins
    motifs_scored["motif_mod"] = motifs_scored["motif"] + "_" + motifs_scored["mod_type"]
    motifs_scored_in_bins = motifs_scored[motifs_scored["motif_mod"].isin(motifs_of_interest)]
    
    # Merge with contig_bins and assembly_stats
    motifs_scored_in_bins = motifs_scored_in_bins.merge(contig_bins, on="contig", how="left")
    motifs_scored_in_bins = motifs_scored_in_bins.merge(assembly_stats[["contig", "length"]], on="contig", how="left")
    
    # Handle NA bins as 'unbinned'
    motifs_scored_in_bins["bin"] = motifs_scored_in_bins["bin"].fillna("unbinned")
    
    # Add bin_contig identifier
    motifs_scored_in_bins["bin_contig"] = motifs_scored_in_bins["bin"] + "_" + motifs_scored_in_bins["contig"] + "_" + motifs_scored_in_bins["length"].astype(str)
    
    # Calculate n_motifs and mean methylation
    motifs_scored_in_bins["n_motifs"] = motifs_scored_in_bins["n_mod"] + motifs_scored_in_bins["n_nomod"]
    motifs_scored_in_bins["mean"] = motifs_scored_in_bins["n_mod"] / motifs_scored_in_bins["n_motifs"]
    
    
    # Remove complement columns:
    # Identify columns that contain the word "complement"
    complement_columns = [col for col in motifs_scored_in_bins.columns if 'complement' in col]

    # Drop these columns from the DataFrame
    motifs_scored_in_bins = motifs_scored_in_bins.drop(columns=complement_columns)

    
    return motifs_scored_in_bins


def construct_bin_motifs_from_motifs_scored_in_bins(motifs_scored_in_bins, motifs_of_interest, args):
    """
    Constructs the bin_motifs_from_motifs_scored_in_bins DataFrame by filtering motifs that are not in bin_motif_binary,
    """
    # create a comparison dataframe from motifs_scored_in_bins and bin_motif_binary
    bin_motifs_from_motifs_scored_in_bins = motifs_scored_in_bins[(motifs_scored_in_bins["bin"] != "unbinned") & motifs_scored_in_bins["motif_mod"].isin(motifs_of_interest) ]\
        .groupby(["bin", "motif_mod"])[["n_mod", "n_nomod"]]\
        .sum()\
        .reset_index()

    bin_motifs_from_motifs_scored_in_bins["n_motifs_bin"] = bin_motifs_from_motifs_scored_in_bins["n_mod"] + bin_motifs_from_motifs_scored_in_bins["n_nomod"]

    # Filter motifs that are not observed more than n_motif_bin_cutoff times
    bin_motifs_from_motifs_scored_in_bins = bin_motifs_from_motifs_scored_in_bins[bin_motifs_from_motifs_scored_in_bins["n_motifs_bin"] > args.n_motif_bin_cutoff] 
    
    # Calculate mean methylation
    bin_motifs_from_motifs_scored_in_bins["mean_methylation"] = bin_motifs_from_motifs_scored_in_bins["n_mod"] / bin_motifs_from_motifs_scored_in_bins["n_motifs_bin"]

    # Calculate standard deviation of methylation per bin and motif_mod
    bin_motifs_mean_and_sd = motifs_scored_in_bins[
        (motifs_scored_in_bins["bin"] != "unbinned") & 
        motifs_scored_in_bins["motif_mod"].isin(motifs_of_interest) & 
        (motifs_scored_in_bins["mean"] > 0.1) &                                 # TODO: Remove this line if the negative cases should be used to determine methylation pattern.
        (motifs_scored_in_bins["n_motifs"] > args.n_motif_contig_cutoff)
        ] \
        .groupby(["bin", "motif_mod"]).agg(
            mean_methylation_bin = pd.NamedAgg(column="mean", aggfunc="mean"),
            std_methylation_bin = pd.NamedAgg(column="mean", aggfunc="std"),
            n_contigs = pd.NamedAgg(column="contig", aggfunc="count"),
        ).reset_index()

    # Fill std_methylation_bin with 0.15/4 for bins with only one contig, which will allow the score to be mean-0.15 at most
    bin_motifs_mean_and_sd["std_methylation_bin"] = bin_motifs_mean_and_sd["std_methylation_bin"].fillna(0.15/4)
    
    
    # Merge with bin_motifs_from_motifs_scored_in_bins
    bin_motifs_from_motifs_scored_in_bins = bin_motifs_from_motifs_scored_in_bins.merge(bin_motifs_mean_and_sd, on=["bin", "motif_mod"], how="left")
    
    bin_motifs_from_motifs_scored_in_bins.to_csv("bin_motifs_from_motifs_scored_in_bins.tsv", index=False, sep="\t")
    
    # TODO: Does this make sense? All motifs with mean_methylation above 0.25 is called as methylated.
    ## Convert mean methylation values to binary
    bin_motifs_from_motifs_scored_in_bins["methylation_binary"] = (
        bin_motifs_from_motifs_scored_in_bins["mean_methylation"] >= 0.14
    ).astype(int)


    # NOTE This makes us unable to find some contigs, which should be included in the analysis
    # ## Remove bins that has no methylated motifs
    # bin_motifs_from_motifs_scored_in_bins = bin_motifs_from_motifs_scored_in_bins.groupby("bin").filter(lambda x: x["methylation_binary"].sum() > 0)
    
    return bin_motifs_from_motifs_scored_in_bins
    


def construct_contig_motif_binary(motifs_scored_in_bins, motifs_of_interest, contig_methylation_cutoff, n_motif_cutoff):
    """
    Constructs the contig_motif_binary DataFrame by filtering motifs that are not in bin_motif_binary,
    filtering motifs that are not observed more than n_motif_cutoff times, and converting mean methylation values to binary.
    
    params:
        motifs_scored_in_bins: DataFrame
        motifs_of_interest: list
        contig_methylation_cutoff: float
        n_motif_cutoff: int
    """
    # Create contig motifs binary
    ## Filter motifs that are not in bin_motif_binary
    contig_motif_binary = motifs_scored_in_bins[motifs_scored_in_bins["motif_mod"].isin(motifs_of_interest)]
    
    ## Filter motifs that are not observed more than n_motif_cutoff times
    contig_motif_binary = contig_motif_binary[contig_motif_binary["n_motifs"] >= n_motif_cutoff]

    ## Convert mean methylation values to binary
    contig_motif_binary["methylation_binary"] = (contig_motif_binary["mean"] >= (contig_methylation_cutoff)).astype(int)

    ## Rename bin_contig to bin
    contig_motif_binary = contig_motif_binary[["bin_contig", "motif_mod", "methylation_binary"]]
    
    # TODO: This makes no sense with two renamings in a row
    
    contig_motif_binary.rename(columns={"bin_contig": "bin"}, inplace=True)

    # Combine bin_motif_binary and contig_motif_binary
    contig_motif_binary = contig_motif_binary.rename(
        columns={
            "bin": "bin_compare",
            "methylation_binary": "methylation_binary_compare",
        }
    )
    
    return contig_motif_binary


def compare_methylation_pattern(motif_binary_compare, choices):
    """
    Compares the methylation pattern between bin and contig motifs and calculates the motif_comparison_score.
    """
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
    
    
    return contig_bin_comparison_score