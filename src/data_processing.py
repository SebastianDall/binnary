import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import gzip
import os
import logging

def load_data(args):
    """
    Load data for the analysis.
    """
    logger = logging.getLogger(__name__)
    logger.info("Loading data.")
    # Load the data from the provided files
    # Load the data from the provided files using Polars
    
    motifs_scored = pl.read_csv(args.motifs_scored, separator="\t")
    
    bin_motifs = pl.read_csv(args.bin_motifs, separator="\t")
    
    # Polars automatically infers headers; if args.contig_bins file doesn't have headers, you need to specify it
    contig_bins = pl.read_csv(args.contig_bins, separator="\t", has_header=False) \
        .rename({
            "column_1": "contig",
            "column_2": "bin"
        })

    logger.info("Data loaded.")
    return motifs_scored, bin_motifs, contig_bins

def read_fasta(file_path):
    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")    
    
    # Check if the file has a valid FASTA extension
    valid_extensions = ['.fasta', '.fa', '.fna', '.gz']
    if not any(file_path.endswith(ext) for ext in valid_extensions):
        raise ValueError(f"Unsupported file extension. Please provide a FASTA file with one of the following extensions: {', '.join(valid_extensions)}")

    # Check if the file is a gzipped FASTA file
    if file_path.endswith('.gz'):
        with gzip.open(file_path, "rt") as handle:  # "rt" mode for reading as text
            return {record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")}
    else:
        # Read a regular (uncompressed) FASTA file
        with open(file_path, "r") as handle:
            return {record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")}


def write_bins_from_contigs(new_contig_bins, assembly_dict, output_dir):
    """
    Creates new bin files based on contig assignments in new_contig_bins DataFrame.

    Args:
    - new_contig_bins (DataFrame): DataFrame with contig and bin assignments.
    - assembly_dict (dict): Dictionary with contig IDs as keys and sequences as values, from the assembly file.
    - output_dir (str): Path to the directory where the new bin files will be saved.
    """
    logger = logging.getLogger(__name__)
    
    # check if the output directory exists
    logger.info(f"Writing bins to {output_dir}/bins/...")
    output_bins_dir = os.path.join(output_dir, "bins")
    if not os.path.exists(output_bins_dir):
        os.makedirs(output_bins_dir)
    
    # Group the DataFrame by 'bin'
    grouped = new_contig_bins.groupby('bin')

    for bin_name, group in grouped:
        # Initialize a list to hold SeqRecord objects for the current bin
        bin_records = []

        for contig in group['contig']:
            # Create a SeqRecord object for each contig in the bin, if it exists in the assembly dictionary
            if contig in assembly_dict:
                seq_record = SeqRecord(Seq(assembly_dict[contig]), id=contig, description="")
                bin_records.append(seq_record)

        # Define the output file name for the current bin
        output_file = f"{bin_name}.fa"
        output_path = os.path.join(output_bins_dir, output_file)

        # Write the SeqRecord objects to a FASTA file
        with open(output_path, "w") as output_handle:
            SeqIO.write(bin_records, output_handle, "fasta")

        logger.info(f"Written {len(bin_records)} contigs to {output_file}")



def generate_output(output_df, outdir, filename, header=True):
    """
    Generate the output files for the analysis.
    """
    # If the output directory does not exist, create it
    if not os.path.exists(outdir) and outdir != "":
        os.makedirs(outdir)
    
    file_path = os.path.join(outdir, filename)
    
    # Generate the output files
    output_df.to_csv(file_path, sep="\t", index=False, header=header)


def prepare_bin_consensus(bin_motifs, args):
    """
    Prepares the bin_consensus_from_bin_motifs DataFrame by calculating the mean methylation per bin and motif_mod and converting it to binary.    
    """
    # Combine 'motif' and 'mod_type' into 'motif_mod'
    bin_motifs = bin_motifs.with_columns((pl.col("motif") + "_" + pl.col("mod_type")).alias("motif_mod"))

    # Calculate total motifs and mean methylation
    bin_motifs = bin_motifs.with_columns([
        (pl.col("n_mod_bin") + pl.col("n_nomod_bin")).alias("n_motifs"),
        (pl.col("n_mod_bin") / (pl.col("n_mod_bin") + pl.col("n_nomod_bin"))).alias("mean_methylation")
    ])

    # Create 'methylation_binary' column based on mean methylation cutoff
    bin_motifs = bin_motifs.with_columns((pl.col("mean_methylation") >= args.mean_methylation_cutoff).cast(int).alias("methylation_binary"))

    # Filter motifs that are not observed more than n_motif_bin_cutoff times
    bin_motifs = bin_motifs.filter(pl.col("n_motifs") >= args.n_motif_bin_cutoff)

    # Filter for rows where 'methylation_binary' is 1 and select relevant columns
    bin_motif_binary = bin_motifs.filter(pl.col("methylation_binary") == 1)[["bin", "motif_mod", "mean_methylation", "methylation_binary"]]

    # Create a binary index for bin and motif_mod
    # Assuming unique_bins and unique_motif_mods are lists of unique values
    unique_bins = bin_motif_binary["bin"].unique().to_list()
    unique_motif_mods = bin_motif_binary["motif_mod"].unique().to_list()

    # Create a DataFrame with all combinations of "bin" and "motif_mod"
    bin_motif_combinations = pl.DataFrame({
        "bin": [bin for bin in unique_bins for _ in unique_motif_mods],
        "motif_mod": [motif_mod for _ in unique_bins for motif_mod in unique_motif_mods]
    })

    # Join the combinations DataFrame with bin_motif_binary
    bin_motif_binary_enriched = bin_motif_combinations.join(
        bin_motif_binary, on=["bin", "motif_mod"], how="left"
    )

    # Fill missing values if necessary and adjust columns as needed
    bin_motif_binary_enriched = bin_motif_binary_enriched.with_columns([
        pl.col("mean_methylation").fill_null(0).alias("mean_methylation"),
        pl.col("methylation_binary").fill_null(0).alias("methylation_binary")
    ])
    

    return bin_motif_binary_enriched


def prepare_motifs_scored_in_bins(motifs_scored, motifs_of_interest, contig_bins):
    """
    Prepares the motifs_scored_in_bins DataFrame by merging with bin motifs, contig bins, and assembly stats,
    and calculates additional metrics like number of motifs and mean methylation per contig.
    """
    
    # Filter and enhance motifs_scored based on motifs_in_bins
    motifs_scored_in_bins = motifs_scored \
        .with_columns(
            (pl.col("motif") + "_" + pl.col("mod_type")).alias("motif_mod"),
            (pl.col("n_mod") + pl.col("n_nomod")).alias("n_motifs")
        ) \
        .filter(pl.col("motif_mod").is_in(motifs_of_interest)) \
        .with_columns((pl.col("n_mod") / pl.col("n_motifs")).alias("mean"))
    
    # Merge with contig_bins
    motifs_scored_in_bins = motifs_scored_in_bins.join(contig_bins, on="contig", how="left") \
        .with_columns(pl.col("bin").fill_null("unbinned")) \
        .with_columns((pl.col("bin") + "_" + pl.col("contig")).alias("bin_contig"))
    
    
    # Remove complement columns:
    # Identify columns that contain the word "complement"
    complement_columns = [col for col in motifs_scored_in_bins.columns if 'complement' in col]

    # Drop these columns from the DataFrame
    motifs_scored_in_bins = motifs_scored_in_bins.drop(complement_columns)

    
    return motifs_scored_in_bins


def remove_ambiguous_motifs_from_bin_consensus(motifs_scored_in_bins, args):
    # Remove motifs in bins where the majority of the mean methylation of motifs is in the range of 0.05-0.4
    contig_motif_mean_density = motifs_scored_in_bins \
        .filter(pl.col("bin") != "unbinned") \
        .with_columns(
            ((pl.col("mean") > 0.05) & (pl.col("mean") < 0.4)).alias("is_ambiguous")
        )

    # Count the number of ambiguous motifs per bin
    # Group by 'bin' and 'motif_mod' and calculate the sum of 'is_ambiguous' and the total count in each group
    bin_consensus_ambiguous_motifs = contig_motif_mean_density.group_by(["bin", "motif_mod"]) \
        .agg(
            pl.col("is_ambiguous").sum().alias("total_ambiguous"),
            pl.col("is_ambiguous").count().alias("n_contigs_with_motif")
        ) \
        .with_columns((pl.col("total_ambiguous") / pl.col("n_contigs_with_motif")).alias("percentage_ambiguous")) \
        .filter(pl.col("percentage_ambiguous") <= args.ambiguous_motif_percentage_cutoff) \
        .select(["bin", "motif_mod"])
    
    
    return bin_consensus_ambiguous_motifs


# TODO: rename - calculate_bin_consensus_from_contigs
def construct_bin_consensus_from_motifs_scored_in_bins(motifs_scored_in_bins, args):
    """
    Constructs the bin_motifs_from_motifs_scored_in_bins DataFrame by filtering motifs that are not in bin_motif_binary,
    """
    
    # Remove motifs in bins where the majority of the mean methylation of motifs is in the range of 0.05-0.4
    bin_consensus_without_ambiguous_motifs = remove_ambiguous_motifs_from_bin_consensus(motifs_scored_in_bins, args)
    
    
    # Find n_motifs in bin TODO: rename bin_motifs_from_motifs_scored_in_bins to bin_consensus_from_motifs_scored_in_bins
    bin_motifs_from_motifs_scored_in_bins = motifs_scored_in_bins[(motifs_scored_in_bins["bin"] != "unbinned")]\
        .groupby(["bin", "motif_mod"])[["n_mod", "n_nomod"]]\
        .sum()\
        .reset_index()

    bin_motifs_from_motifs_scored_in_bins["n_motifs_bin"] = bin_motifs_from_motifs_scored_in_bins["n_mod"] + bin_motifs_from_motifs_scored_in_bins["n_nomod"]

    # Filter motifs that are not observed more than n_motif_bin_cutoff times
    bin_motifs_from_motifs_scored_in_bins = bin_motifs_from_motifs_scored_in_bins[bin_motifs_from_motifs_scored_in_bins["n_motifs_bin"] > args.n_motif_bin_cutoff] 
    
    # Calculate mean methylation
    bin_motifs_from_motifs_scored_in_bins["mean_methylation"] = bin_motifs_from_motifs_scored_in_bins["n_mod"] / bin_motifs_from_motifs_scored_in_bins["n_motifs_bin"]

    # retain only motifs found in bin_consensus_without_ambiguous_motifs
    bin_motifs_from_motifs_scored_in_bins = bin_motifs_from_motifs_scored_in_bins.merge(bin_consensus_without_ambiguous_motifs, on=["bin", "motif_mod"], how="inner")

    # Calculate standard deviation of methylation per bin and motif_mod
    bin_motifs_mean_and_sd = motifs_scored_in_bins[
        (motifs_scored_in_bins["bin"] != "unbinned") &  
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
    
    # TODO: Does this make sense? All motifs with mean_methylation above 0.25 is called as methylated.
    ## Convert mean methylation values to binary
    bin_motifs_from_motifs_scored_in_bins["methylation_binary"] = (
        bin_motifs_from_motifs_scored_in_bins["mean_methylation"] >= args.mean_methylation_cutoff
    ).astype(int)
    
    return bin_motifs_from_motifs_scored_in_bins
    

def calculate_binary_motif_comparison_matrix(motifs_scored_in_bins, args):
    # Step 1 create bin_motif_from_motifs_scored_in_bins - basis for bin-contig comparison
    bin_motifs_from_motifs_scored_in_bins = construct_bin_motifs_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    ## Filter motifs that are not observed more than n_motif_cutoff times
    motifs_scored_in_contigs = motifs_scored_in_bins[motifs_scored_in_bins["n_motifs"] >= args.n_motif_contig_cutoff]   
    
    ## Rename bin_contig to bin
    motifs_scored_in_contigs = motifs_scored_in_contigs[["bin_contig", "motif_mod", "mean"]]
    motifs_scored_in_contigs.rename(columns={"bin_contig": "bin_compare"}, inplace=True)
    

    # Merge bin_motifs_from_motifs_scored_in_bins and motifs_scored_in_contigs    
    motif_binary_compare = pd.merge(
        bin_motifs_from_motifs_scored_in_bins, motifs_scored_in_contigs, on="motif_mod"
    )
    
    
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



def load_contamination_file(contamination_file):
    """
    Load the contamination file from the provided path.
    """
    contamination = pd.read_csv(contamination_file, delimiter = "\t")
    
    # Check if the file contains the required columns
    # bin	bin_contig_compare	binary_methylation_missmatch_score	contig	alternative_bin	alternative_bin_binary_methylation_missmatch_score
    required_columns = ["bin", "bin_contig_compare", "binary_methylation_missmatch_score", "contig", "alternative_bin", "alternative_bin_binary_methylation_missmatch_score"]
    if not all(column in contamination.columns for column in required_columns):
        raise ValueError("The contamination file does not contain the required columns.")
    
    # Check if the file contains any rows
    if len(contamination) == 0:
        raise ValueError("The contamination file is empty.")
    
    return contamination


def create_contig_bin_file(contig_bins, include, contamination):
    """
    Create a new contig_bin file based on the analysis results and contamination file.
    """
    # Remove contigs in the contamination file from the contig_bins
    contig_bins = contig_bins[~contig_bins["contig"].isin(contamination["contig"])]
    
    # Add the contigs in the include DataFrame to the contig_bins
    contig_bins = pd.concat([contig_bins, include[["contig", "bin"]]], ignore_index=True)
    
    # Sort the contig_bins by bin and contig
    contig_bins = contig_bins.sort_values(by=["bin", "contig"])
    
    return contig_bins