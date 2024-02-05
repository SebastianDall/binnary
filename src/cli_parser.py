import argparse


def add_common_arguments(subparser):
    """Function to add common arguments to a subparser."""
    subparser.add_argument(
        "--motifs_scored", type=str, help="Path to motifs-scored.tsv from nanomotif"
    )
    subparser.add_argument("--bin_motifs", type=str, help="Path to bin-motifs.tsv file")
    subparser.add_argument(
        "--contig_bins", type=str, help="Path to bins.tsv file for contig bins"
    )
    subparser.add_argument(
        "--assembly_stats", type=str, help="Path to assembly_info.txt file"
    )
    # subparser.add_argument("--assembly_file", type=str, help="Path to assembly.fasta file")
    subparser.add_argument(
        "--mean_methylation_cutoff",
        type=float,
        default=0.75,
        help="Cutoff value for considering a motif as methylated",
    )
    subparser.add_argument(
        "--n_motif_cutoff",
        type=int,
        default=6,
        help="Number of motifs observed before considering a motif valid",
    )
    subparser.add_argument("--out", type=str, help="Path to output filename")


def get_parser():
    # Create the top-level parser
    parser = argparse.ArgumentParser(description="DNA Methylation Pattern Analysis")
    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")

    # Create subparsers for each command
    commands = ["detect_contamination", "bin_unbinned", "bin_contamination"]
    for command in commands:
        subparser = subparsers.add_parser(command, help=f"{command} help")
        add_common_arguments(subparser)  # Add common arguments to each subparser

    return parser
