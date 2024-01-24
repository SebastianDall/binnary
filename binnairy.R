#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(readr)
library(dplyr)


# Define the options
option_list <- list(
  make_option(c("--motifs"), type = "character", default = NULL, help = "Path to motifs.tsv file", metavar = "FILE"),
  make_option(c("--motifs_scored"), type = "character", default = NULL, help = "Path to motifs-scored.tsv file", metavar = "FILE"),
  make_option(c("--bin_motifs"), type = "character", default = NULL, help = "Path to bin-motifs.tsv file", metavar = "FILE"),
  make_option(c("--contig_bins"), type = "character", default = NULL, help = "Path to bins.tsv file for contig bins", metavar = "FILE"),
  make_option(c("--bin_stats"), type = "character", default = NULL, help = "Path to bin-stats.tsv file", metavar = "FILE"),
  make_option(c("--assembly_stats"), type = "character", default = NULL, help = "Path to assembly_info.txt file", metavar = "FILE"),
  make_option(c("--assembly_file"), type = "character", default = NULL, help = "Path to assembly.fasta file", metavar = "FILE")
)

# Create a parser
parser <- OptionParser(option_list = option_list)

# Parse the arguments
arguments <- parse_args(parser)

# Read and process the files based on the provided arguments
if (!is.null(arguments$motifs)) {
  motifs <- read_tsv(arguments$motifs)
}

if (!is.null(arguments$motifs_scored)) {
  motifs_scored <- read_tsv(arguments$motifs_scored)
}

if (!is.null(arguments$bin_motifs)) {
  bin_motifs <- read_tsv(arguments$bin_motifs) %>%
    mutate(bin = remove_sample_barcode(bin))  # Modify or remove this line based on your needs
}

if (!is.null(arguments$contig_bins)) {
  contig_bins <- read_tsv(arguments$contig_bins, col_names = c("contig", "bin")) %>%
    mutate(bin = remove_sample_barcode(bin))  # Modify or remove this line based on your needs
}

if (!is.null(arguments$bin_stats)) {
  bin_stats <- read_tsv(arguments$bin_stats) %>%
    mutate(bin = remove_sample_barcode(bin))  # Modify or remove this line based on your needs
}

if (!is.null(arguments$assembly_stats)) {
  assembly_stats <- read_tsv(arguments$assembly_stats) %>%
    rename(contig = `#seq_name`)
}

if (!is.null(arguments$assembly_file)) {
  assembly_stats <- read_tsv(arguments$assembly_stats) %>%
    rename(contig = `#seq_name`)
}

# Add further processing or output here
