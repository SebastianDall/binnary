convert_to_binary <- function(df, MEAN_METHYLATION_CUTOFF){
  df <- df %>% 
    summarise(
      mean_methylation = mean(mean)
    ) %>% 
    mutate(
      methylation_binary = if_else(mean_methylation >= MEAN_METHYLATION_CUTOFF, 1, 0)
    ) %>% 
    select(!mean_methylation)
  
  return(df)
}


find_bins_w_no_methylation <- function(df, bin_motif_binary){
  no_methylation_bin_present <- bin_motif_binary %>% 
    group_by(bin) %>% 
    summarise(
      sum_methylation = sum(methylation_binary, na.rm = TRUE), .groups = "drop"
    ) %>% 
    filter(sum_methylation == 0)
  
  if (length(no_methylation_bin_present$bin) == 0) {
    bin_contigs_w_0_or_NA_only <- df %>% 
      group_by(bin) %>% 
      summarise(
        sum_methylation = sum(methylation_binary, na.rm = TRUE), .groups = "drop"
      ) %>% 
      filter(sum_methylation == 0)
    
    df <- df %>% 
      filter(!bin %in% bin_contigs_w_0_or_NA_only$bin)
    
    return(df)
    
  } else {
    return(df)
  }
}