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


# Suppressing both messages and warnings
quiet_read_tsv <- function(file_path) {
  suppressWarnings(
    suppressMessages(
      read_tsv(file_path, progress = FALSE)
    )
  )
}