#!/bin/env python3
import sys
from src import data_processing, analysis
from src.cli_parser import get_parser

# Import other necessary libraries here


def main(args):
    """
    Main entry point for the DNA Methylation Pattern Analysis tool.
    Orchestrates the workflow of the tool based on the provided arguments.
    """
    print("Starting Binnary Analysis...")

    # Step 1: Load and preprocess data
    # These functions would be defined in your data_processing module
    (
        motifs_scored,
        bin_motifs,
        contig_bins,
        assembly_stats,
        assembly_file,
    ) = data_processing.load_data(args)

    # Step 2: Perform core analysis
    # Functions from the analysis module
    analysis_results = analysis.perform_analysis(
        motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file, args
    )

    # Step 3: Post-analysis processing and output generation
    # data_processing.generate_output(analysis_results, args.out)
    data_processing.generate_output(analysis_results, args.out)
    print("Analysis Completed. Results are saved to:", args.out)


if __name__ == "__main__":
    parser = get_parser()  # Get the configured argument parser

    # Check if no arguments were provided (only the script name is present)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)  # Print help message to stderr
        sys.exit(
            1
        )  # Exit with a non-zero status to indicate that the program was not executed as intended

    args = parser.parse_args()  # Parse the command-line arguments
    main(args)
