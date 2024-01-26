

determine_contamination <- function(belonging_score) {
  contamination_df <- belonging_score %>% 
    mutate(
      group = case_when(
        bin == bin_id & belonging_bins == 1 ~ "correct",
        bin != bin_id & belonging_bins == 1 ~ "single_contamination",
        belonging_bins > 1 ~ "multi_contamination"
      )
    )
  
  return(contamination_df)
}




extract_kmer_dist <- function(dist_matrix) {
  dist_matrix %>% 
    # Convert to dataframe
    as.matrix() %>% 
    as.data.frame() %>% 
    # Isolate contig distance
    rownames_to_column("contig") %>% 
    filter(str_detect(contig, "contig")) %>% 
    select(!contains("contig_")) %>% 
    pivot_longer(-contig) %>% 
    # Summarise distance
    mutate(
      value = round(value, 5),
      matching_bins = paste0(name,"|",value)
    ) %>% 
    select(contig,matching_bins) %>% 
    group_by(contig) %>%
    summarize(matching_bins = paste(matching_bins, collapse = ";"), .groups = "drop")
}


calculate_kmer_frequency <- function(df, fasta, kmer_window) {
  df <- df %>% 
    mutate(
      contig_seq = map(.x = contig, ~fasta[[paste0(.x)]]),
      kmer_count = map(.x = contig_seq, ~(as.data.frame(seqinr::count(.x, wordsize = kmer_window)) %>% pivot_wider(names_from = Var1, values_from = Freq)))
    ) %>% 
    unnest(kmer_count) %>% 
    select(!contig_seq)
}



# TODO: Implement logic if no contigs could be assigned to multiple bins

assign_contamination_from_kmer <- function(
  belonging_score,
  fasta,
  contig_bins,
  kmer_window
) {
  print("Calculating Bin kmer frequency")
  kmer_bin_df <- contig_bins %>% 
    calculate_kmer_frequency(fasta, kmer_window)
    # mutate(
    #   contig_seq = map(.x = contig, ~fasta[[paste0(.x)]]),
    #   kmer_count = map(.x = contig_seq, ~(as.data.frame(seqinr::count(.x, wordsize = kmer_window)) %>% pivot_wider(names_from = Var1, values_from = Freq)))
    # ) %>% 
    # unnest(kmer_count) %>% 
    # select(!contig_seq)
  
  # Remove unbinned contigs as kmer frequency does not make sense.
  contamination_multi <- belonging_score %>% 
    filter(!str_detect(bin_compare, "unbinned")) %>% 
    determine_contamination() %>% 
    filter(
      group == "multi_contamination"
    )
  
  if (length(contamination_multi$bin_compare) == 0) {
    print("No contigs could be assigned to multiple bins")
    return(tibble(contig = c(), matching_bins = c()))
  }
  print("Calculating contig kmer frequency")
  kmer_contig_df <- tibble(contig = contamination_multi %>% pull(contig) %>% unique()) %>% 
    # calculate_kmer_frequency(fasta)
    left_join(kmer_bin_df) %>%
    select(!bin)
  
  # print(kmer_contig_df)
  
  cat("Calculating kmer distance to bins\n")
  contig_bins_out_tmp <- contamination_multi %>% 
    ungroup() %>% 
    mutate(
      # Extract kmer bins, where contig could belong (remove contig if found in bin)
      bin_kmer = map2(
        .x = contig, .y = bin, 
        ~(
          kmer_bin_df %>% 
            filter(contig != .x, bin == .y) %>% 
            mutate(contig = .y) %>% 
            select(-bin)
        )
      ),
      # Extract kmer for contig
      bin_contig_kmer = map2(
        .x = bin_kmer, .y = contig, 
        ~(
          kmer_contig_df %>% 
            filter(contig == .y) %>%
            bind_rows(.x)
        )
      ),
      # Summarise kmer counts for the bins and convert to frequency
      bin_contig_kmer = map(
        .x = bin_contig_kmer, 
        ~(
          .x %>% 
            pivot_longer(-contig, names_to = "kmer", values_to = "count") %>% 
            group_by(contig, kmer) %>% 
            summarise(count = sum(count), .groups = "drop") %>% 
            # Calculate total count per contig
            mutate(total_count = sum(count)) %>%
            # Calculate frequency
            mutate(frequency = count / total_count) %>%
            # Remove the total_count column
            select(-c(total_count, count)) %>%
            pivot_wider(names_from = "kmer", values_from = "frequency")
        )
      )
    ) %>% 
    ungroup() %>% 
    select(contig, bin_contig_kmer) %>% 
    arrange(contig)
  
  
  contig_bins_out <- contig_bins_out_tmp %>% 
    rename(unbinned = contig) %>% 
    unnest(bin_contig_kmer) %>%
    group_by(unbinned) %>% 
    # Remove multiple entries for contig
    distinct(contig, .keep_all = TRUE) %>% 
    ungroup() %>% 
    # Print matching bins and distance
    nest(.by = unbinned) %>% 
    mutate(
      df_matrix = map(
        .x = data, 
        ~(.x %>% column_to_rownames("contig"))
      ),
      dist_matrix = map(
        .x = df_matrix, 
        ~(.x %>% dist())
      ),
      kmer_dist = map(.x = dist_matrix, .y = unbinned, ~extract_kmer_dist(.x))
    ) %>% 
    select(kmer_dist) %>% 
    unnest(kmer_dist)
  
  return(contig_bins_out)
}

# 
# construct_contig_bin_out <- function(
#   belonging_score
# )



