#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
suppressPackageStartupMessages(library(tidyverse))
# source("./R/setup.R")

# Install the here package if not already installed
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

# Use the here package
library(here)
# Now here() will always point to your project's top-level directory

here()

# Construct the path
setup_path <- here("subprojects", "3.methylation_calling", "tools","binnairy" ,"R", "setup.R")

# Source the setup.R file
source(setup_path)

common_path <- here("subprojects", "3.methylation_calling", "tools","binnairy" ,"R", "common.R")
source(common_path)

# Define the options
option_list <- list(
  make_option(c("--mode"), type = "character", default = NULL, help = "Either unbinned or contamination", metavar = "MODE"),
  # make_option(c("--motifs"), type = "character", default = NULL, help = "Path to motifs.tsv file", metavar = "FILE"),
  make_option(c("--motifs_scored"), type = "character", default = NULL, help = "Path to motifs-scored.tsv file", metavar = "FILE"),
  make_option(c("--bin_motifs"), type = "character", default = NULL, help = "Path to bin-motifs.tsv file", metavar = "FILE"),
  make_option(c("--contig_bins"), type = "character", default = NULL, help = "Path to bins.tsv file for contig bins", metavar = "FILE"),
  # make_option(c("--bin_stats"), type = "character", default = NULL, help = "Path to bin-stats.tsv file", metavar = "FILE"),
  make_option(c("--assembly_stats"), type = "character", default = NULL, help = "Path to assembly_info.txt file", metavar = "FILE"),
  make_option(c("--assembly_file"), type = "character", default = NULL, help = "Path to assembly.fasta file", metavar = "FILE"),
  make_option(c("--mean_methylation_cutoff"), type = "double", default = 0.25, help = "Cutoff value for considering a motif as methylated", metavar = "FLOAT"),
  make_option(c("--n_motif_cutoff"), type = "integer", default = 6, help = "Number of motifs observed before considering a motif valid", metavar = "INTEGER")
)

# Create a parser
parser <- OptionParser(option_list = option_list)

# Parse the arguments
arguments <- parse_args(parser)

# Check mode
if (!is.null(arguments$mode)) {
  mode <- arguments$mode
  
  if (!mode %in% c("contamination", "unbinned")) {
    stop(paste0("Mode should be either 'unbinned' or 'contamination', got '", mode, "'"), call. = FALSE)
  }
}


# Read and process the files based on the provided arguments
if (!is.null(arguments$motifs)) {
  motifs <- read_tsv(arguments$motifs)
}

if (!is.null(arguments$motifs_scored)) {
  motifs_scored <- read_tsv(arguments$motifs_scored)
}

if (!is.null(arguments$bin_motifs)) {
  bin_motifs <- read_tsv(arguments$bin_motifs) 
}

if (!is.null(arguments$contig_bins)) {
  contig_bins <- read_tsv(arguments$contig_bins, col_names = c("contig", "bin")) 
}

# if (!is.null(arguments$bin_stats)) {
#   bin_stats <- read_tsv(arguments$bin_stats) %>%
#     mutate(bin = remove_sample_barcode(bin))  # Modify or remove this line based on your needs
# }

if (!is.null(arguments$assembly_stats)) {
  assembly_stats <- read_tsv(arguments$assembly_stats) %>%
    rename(contig = `#seq_name`)
}

if (!is.null(arguments$assembly_file)) {
  assembly_file <- seqinr::read.fasta(arguments$assembly_file, seqtype = "DNA")
}

# Add further processing or output here
## Setup
motifs_scored_in_bins <- setup_motifs_scored_in_bins(
  bin_motifs = bin_motifs,
  motifs_scored = motifs_scored,
  contig_bins = contig_bins,
  assembly_stats = assembly_stats
)


belonging_score <- calculate_belonging_score(
  motifs_scored_in_bins = motifs_scored_in_bins,
  MEAN_METHYLATION_CUTOFF = arguments$mean_methylation_cutoff,
  N_MOTIFS_CUTOFF = arguments$n_motif_cutoff
)





# if (mode == "unbinned") {
#   # run_unbinned(
#   #   # # motifs,
#   #   # motifs_scored,
#   #   # bins_motifs,
#   #   # contig_bins
#   # )
# }



