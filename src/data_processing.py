import pandas as pd
from Bio import SeqIO

def load_data(args):
    """
    Load and preprocess data for the analysis.
    """
    # Load the data from the provided files
    motifs_scored = read_csv(args.motifs_scored, delimeter = "\t")
    bin_motifs = read_csv(args.bin_motifs, delimeter = "\t")
    contig_bins = read_csv(args.contig_bins, delimeter = "\t", header = None)
    assembly_stats = read_csv(args.assembly_stats, delimeter = "\t")
    
    print("Loading assembly file...")
    assembly_file = read_fasta(args.assembly_file)

    # Perform any additional preprocessing steps here
    # Change colnames of contig_bins
    contig_bins.columns = ["contig", "bin"]

    return motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file



def read_fasta(file_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}
