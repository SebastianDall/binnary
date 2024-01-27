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
    
    print("Loading assembly file...")
    assembly_file = read_fasta(args.assembly_file)

    # Perform any additional preprocessing steps here
    # Change colnames of contig_bins
    contig_bins.columns = ["contig", "bin"]
    
    # Rename the first column to 'contig'
    assembly_stats.rename(columns={assembly_stats.columns[0]: "contig"}, inplace=True)
    
    return motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file


def read_fasta(file_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}


def prepare_motifs_scored_in_bins(motifs_scored, bin_motifs, contig_bins, assembly_stats):
    """
    Prepares the motifs_scored_in_bins DataFrame by merging with bin motifs, contig bins, and assembly stats,
    and calculates additional metrics like number of motifs and mean methylation per contig.
    """
    # Create find bin motifs
    bin_motifs["motif_mod"] = bin_motifs["motif"] + "_" + bin_motifs["mod_type"]
    motifs_in_bins = bin_motifs["motif_mod"].unique()
    
    # Filter and enhance motifs_scored based on motifs_in_bins
    motifs_scored["motif_mod"] = motifs_scored["motif"] + "_" + motifs_scored["mod_type"]
    motifs_scored_in_bins = motifs_scored[motifs_scored["motif_mod"].isin(motifs_in_bins)]
    
    # Merge with contig_bins and assembly_stats
    motifs_scored_in_bins = motifs_scored_in_bins.merge(contig_bins, on="contig")
    motifs_scored_in_bins = motifs_scored_in_bins.merge(assembly_stats[["contig", "length"]], on="contig")
    
    # Handle NA bins as 'unbinned'
    motifs_scored_in_bins["bin"] = motifs_scored_in_bins["bin"].fillna("unbinned")
    
    # Add bin_contig identifier
    motifs_scored_in_bins["bin_contig"] = motifs_scored_in_bins["bin"] + "_" + motifs_scored_in_bins["contig"] + "_" + motifs_scored_in_bins["length"].astype(str)
    
    # Calculate n_motifs and mean methylation
    motifs_scored_in_bins["n_motifs"] = motifs_scored_in_bins["n_mod"] + motifs_scored_in_bins["n_nomod"]
    motifs_scored_in_bins["mean"] = motifs_scored_in_bins["n_mod"] / motifs_scored_in_bins["n_motifs"]
    
    return motifs_scored_in_bins
