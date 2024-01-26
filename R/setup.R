# source("./R/common.R")

#### Setup 
setup_motifs_scored_in_bins <- function (
  bin_motifs = bin_motifs,
  motifs_scored = motifs_scored,
  contig_bins = contig_bins,
  assembly_stats = assembly_stats
) {
  
  # Get motifs associated with bins. Only methylations found in bins are considered as they are the basis for scoring contamination and binning
  motifs_in_bins <- bin_motifs %>%
    mutate(motif_mod = paste0(motif,"_",mod_type)) %>% 
    pull(motif_mod) %>%
    unique()
  
  # Calculate motif methylation degree for all contigs.
  motifs_scored_in_bins <- motifs_scored %>%
    mutate(motif_mod = paste0(motif,"_",mod_type)) %>% 
    
    # Only consider motifs scored in bins
    filter(motif_mod %in% motifs_in_bins) %>%
    
    # Add assembly stats
    left_join(contig_bins, by = "contig") %>%
    left_join(assembly_stats %>% select(contig, length), by = "contig") %>%
    
    # If bin is NA it is considered as unbinned
    mutate(
      bin = if_else(is.na(bin), "unbinned", bin),
      bin_contig = paste(bin, contig, length, sep = "_")
    ) %>%
    
    # Calculate mean
    mutate(
      n_motifs = n_mod + n_nomod
    ) %>% 
    mutate(
      mean = n_mod / n_motifs
    )
  
  return(motifs_scored_in_bins)
  
}



calculate_belonging_score <- function(
  motifs_scored_in_bins = motifs_scored_in_bins,
  MEAN_METHYLATION_CUTOFF = 0.25,
  N_MOTIFS_CUTOFF = 6
) {
  
  # Convert bin motifs to binary
  bin_motif_binary <- motifs_scored_in_bins %>% 
    filter(bin != "unbinned") %>% 
    group_by(bin, motif_mod) %>% 
    convert_to_binary(MEAN_METHYLATION_CUTOFF) %>% 
    ungroup()
  
  # Convert contig methylation to binary
  contig_motif_binary <- motifs_scored_in_bins %>% 
    # Consider only motifs found in bins
    filter(motif_mod %in% bin_motif_binary$motif_mod %>% unique()) %>% 
    group_by(bin_contig, motif_mod) %>% 
    # Place a filter for number of motifs before. Few motifs can by noise have a huge impact
    filter(n_motifs > N_MOTIFS_CUTOFF) %>%
    # Convert to binary
    convert_to_binary(MEAN_METHYLATION_CUTOFF) %>% 
    rename(bin = bin_contig) %>% 
    # Add NA values
    pivot_wider(names_from = motif_mod, values_from = methylation_binary) %>% 
    pivot_longer(-bin, names_to = "motif_mod", values_to = "methylation_binary") %>% 
    # If there is no bin with no methylation pattern, then remove contigs with no methylation pattern
    find_bins_w_no_methylation(bin_motif_binary)
  
  
  # Compute scoring metrics:
  ## TODO: Implement more sophisticated scoring system taking n_motifs into account.
  ## TODO: Deferentiate scoring system based on no methylations. If there is no bins with no methylation patterns, then remove all contigs with no methylation pattern. 
  
  
  
  motif_binary_compare <- bin_motif_binary %>%
    left_join(
      contig_motif_binary  %>% rename(bin_compare = bin, methylation_binary_compare = methylation_binary),
      relationship = "many-to-many"
    ) %>%
    mutate(
      # motif_comparison = methylation_binary == methylation_binary_compare,
      # motif_exists_where = case_when(
      #   methylation_binary == methylation_binary_compare ~ "both",
      #   methylation_binary == 1 & methylation_binary_compare == 0 ~ "bin",
      #   methylation_binary == 0 & methylation_binary_compare == 1 ~ "other",
      #   methylation_binary == 0 & methylation_binary_compare == 0 ~ "neither",
      #   methylation_binary == 1 & is.na(methylation_binary_compare) ~ "unknown",
      #   methylation_binary == 0 & is.na(methylation_binary_compare) ~ "unknown"
      # ),
      motif_comparison_score = case_when(
        methylation_binary == 1 & methylation_binary_compare == 1 ~ 1,    # Methylation on both contig and bin +1
        methylation_binary == 1 & methylation_binary_compare == 0 ~ -1,   # Methylation missing on contig -1
        methylation_binary == 0 & methylation_binary_compare == 1 ~ 0,    # Methylation only found on contig 0 (can be due to noise)
        methylation_binary == 0 & methylation_binary_compare == 0 ~ 0,    # Methylation missing on both contig and bin 0
        methylation_binary == 1 & is.na(methylation_binary_compare) ~ 0,   # Methylation missing on contig due to NA 0 (No penalty means benefit of the doubt)
        methylation_binary == 0 & is.na(methylation_binary_compare) ~ 0   # Methylation missing on contig due to NA 0 (No penalty means benefit of the doubt)
      )
    )
  
  
  
  belonging_score <- motif_binary_compare %>% 
    # Summarize the scores
    group_by(bin, bin_compare) %>% 
    summarise(n = sum(motif_comparison_score), .groups = "drop") %>%
    # Find the max per bin_compare
    group_by(bin_compare) %>%
    filter(n == max(n)) %>%
    mutate(belonging_bins = n()) %>% 
    separate(bin_compare, into = c("bin_id", "contig", "contig_number", "length"), sep = "_", remove = FALSE) %>% 
    mutate(
      contig = paste0("contig_", contig_number),
      length = as.numeric(length)
    ) %>% 
    select(!contig_number)
  
  return(belonging_score)
}






