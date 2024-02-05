import pandas as pd
from Bio import SeqIO

def load_data(args):
    """
    Load and preprocess data for the analysis.
    """
    # Load the data from the provided files
    motifs_scored = pd.read_csv(args.motifs_scored, delimiter = "\t")
    bin_motifs = pd.read_csv(args.bin_motifs, delimiter = "\t")
    contig_bins = pd.read_csv(args.contig_bins, delimiter = "\t", header = None)
    assembly_stats = pd.read_csv(args.assembly_stats, delimiter = "\t")
    
    
    # print("Loading assembly file...")
    # assembly_file = read_fasta(args.assembly_file)
    # print("Assembly file loaded.")    

    # Perform any additional preprocessing steps here
    # Change colnames of contig_bins
    contig_bins.columns = ["contig", "bin"]
    
    # Rename the first column to 'contig'
    assembly_stats.rename(columns={assembly_stats.columns[0]: "contig"}, inplace=True)
    
    return motifs_scored, bin_motifs, contig_bins, assembly_stats


def read_fasta(file_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}



def generate_output(output_df, output_path):
    """
    Generate the output files for the analysis.
    """
    # Generate the output files
    output_df.to_csv(output_path, sep="\t", index=False)


def prepare_motifs_scored_in_bins(motifs_scored, bin_motif_binary, contig_bins, assembly_stats):
    """
    Prepares the motifs_scored_in_bins DataFrame by merging with bin motifs, contig bins, and assembly stats,
    and calculates additional metrics like number of motifs and mean methylation per contig.
    """
    # Create find bin motifs
    motifs_in_bins = bin_motif_binary["motif_mod"].unique()
    
    # Filter and enhance motifs_scored based on motifs_in_bins
    motifs_scored["motif_mod"] = motifs_scored["motif"] + "_" + motifs_scored["mod_type"]
    motifs_scored_in_bins = motifs_scored[motifs_scored["motif_mod"].isin(motifs_in_bins)]
    
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


def prepare_bin_motifs_binary(bin_motifs, args):
    """
    Prepares the bin_motifs_binary DataFrame by calculating the mean methylation per bin and motif_mod and converting it to binary.    
    """
    # Alter bin_motifs to include motif_mod and mean
    bin_motifs["motif_mod"] = bin_motifs["motif"] + "_" + bin_motifs["mod_type"]
    # Calculate n_motifs and mean methylation
    bin_motifs["n_motifs"] = bin_motifs["n_mod_bin"] + bin_motifs["n_nomod_bin"]
    bin_motifs["mean"] = bin_motifs["n_mod_bin"] / bin_motifs["n_motifs"]
    
    
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
    
    ## Step 3: Only keep motifs with methylation binary == 1
    bin_motif_binary = bin_motif_binary[bin_motif_binary["methylation_binary"] == 1]

    ## Step 4: Create a binary index for bin and motif_mod
    bin_motif_binary_index = pd.MultiIndex.from_product([bin_motif_binary["bin"].unique(), bin_motif_binary["motif_mod"].unique()], names=['bin', 'motif_mod'])

    bin_motif_binary = bin_motif_binary.set_index(["bin", "motif_mod"])

    bin_motif_binary = bin_motif_binary.reindex(bin_motif_binary_index, fill_value=0).reset_index()
    
    return bin_motif_binary
