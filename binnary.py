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
    
    # Conditional check for --write_bins and --assembly_file
    if args.command == "include_contigs":
        if args.write_bins and not args.assembly_file:
            print("Error: --assembly_file must be specified when --write_bins is used.")
            sys.exit(1)
        elif not args.write_bins and args.assembly_file:
            print("Error: --assembly_file can only be used when --write_bins is specified.")
            sys.exit(1)
    
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
        
        # Save the include_contigs_df results
        data_processing.generate_output(include_contigs_df, args.out, "include_contigs.tsv")
        
        # Create a new contig_bin file
        new_contig_bins = data_processing.create_contig_bin_file(contig_bins, include_contigs_df, contamination)
        data_processing.generate_output(new_contig_bins, args.out, "decontaminated_contig_bins_with_include.tsv")
        
        if args.write_bins:
            print("Loading assembly file...")
            assembly = data_processing.read_fasta(args.assembly_file)
            
            print(f"Writing bins to {args.out}/bins/...")
            data_processing.write_bins_from_contigs(new_contig_bins, assembly, args.out)
            
            
            
            
    
    
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
