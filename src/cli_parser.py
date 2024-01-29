import argparse

def get_parser():
    """Constructs and returns the argument parser for the application."""
    parser = argparse.ArgumentParser(description="DNA Methylation Pattern Analysis")

    # Adding command-line options
    parser.add_argument("--motifs_scored", type=str, help="Path to motifs-scored.tsv file")
    parser.add_argument("--bin_motifs", type=str, help="Path to bin-motifs.tsv file")
    parser.add_argument("--contig_bins", type=str, help="Path to bins.tsv file for contig bins")
    parser.add_argument("--assembly_stats", type=str, help="Path to assembly_info.txt file")
    parser.add_argument("--assembly_file", type=str, help="Path to assembly.fasta file")
    parser.add_argument("--mean_methylation_cutoff", type=float, default=0.25, help="Cutoff value for considering a motif as methylated")
    parser.add_argument("--n_motif_cutoff", type=int, default=6, help="Number of motifs observed before considering a motif valid")
    parser.add_argument("--kmer_window_size", type=int, default=4, help="kmer window size")
    parser.add_argument("--out", type=str, default="contig-bin-association.tsv", help="Path to output filename")

    return parser
