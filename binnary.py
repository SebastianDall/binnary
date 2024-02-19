#!/usr/bin/env python3
import sys
from src import data_processing, detect_contamination, include_contigs
from src.cli_parser import get_parser

# Import other necessary libraries here


def main(args):
    """
    Main entry point for the DNA Methylation Pattern Analysis tool.
    Orchestrates the workflow of the tool based on the provided arguments.
    """
    print("Starting Binnary ", args.command, " analysis...")

    # Step 1: Load and preprocess data
    # These functions would be defined in your data_processing module
    (
        motifs_scored,
        bin_motifs,
        contig_bins,
        assembly_stats,
        # assembly_file,
    ) = data_processing.load_data(args)

    

    # Step 2: create motifs_scored_in_bins and bin_motif_binary
    bin_motif_binary = data_processing.calculate_binary_methylation_bin_consensus_from_bin_motifs(bin_motifs, args)
    
    motifs_in_bin_consensus = bin_motif_binary["motif_mod"].unique()
    
    motifs_scored_in_bins = data_processing.prepare_motifs_scored_in_bins(
        motifs_scored,
        motifs_in_bin_consensus,
        contig_bins,
        assembly_stats,
    )

    # Functions from the analysis module
    if args.command == "detect_contamination" or (args.command == "include_contigs" and args.run_detect_contamination):
        contamination = detect_contamination.detect_contamination(
            motifs_scored_in_bins, args
        )
        data_processing.generate_output(contamination, args.out, "bin_contamination.tsv")
        

    if args.command == "include_contigs":
        # User provided contamination file
        if args.contamination_file:
            print("Loading contamination file...")
            contamination = data_processing.load_contamination_file(args.contamination_file)
        
        # Run the include_contigs analysis    
        include_contigs_df = include_contigs.include_contigs(
            motifs_scored_in_bins, contamination, args
        )
        
        # Create a new contig_bin file
        data_processing.create_contig_bin_file(contig_bins, include_contigs_df, contamination, args.out)
        
        # Save the include_contigs_df results
        data_processing.generate_output(include_contigs_df, args.out, "include_contigs.tsv")
    
    
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
