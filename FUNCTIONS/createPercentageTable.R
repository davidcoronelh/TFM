createPercentageTable <- function(data, start_position, desired_positions = 20, cancer_type) {
  
  require(purrr)
  require(dplyr)
  require(stringi)
  require(tidyr)
  
  if ("CDS" %in% unique(data$class)) {
    class <- "CDS"
  } else {
    class <- "non_CDS"
  }
  
  valid_amino_acids <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
  
  if (start_position == "first") { 
    
    data <- data %>%
      filter(class != "uoORF")
    
    position_freq <- vector("list", desired_positions)
    names(position_freq) <- 1:desired_positions
    
  } else if (start_position == "last") {
    
    data <- data %>%
      mutate(seq = stringi::stri_reverse(seq)) %>% 
      filter(class != "doORF")
    
    position_freq <- vector("list", desired_positions)
    names(position_freq) <- (1:desired_positions - 1) * -1
    
  }
  
  for (i in 1:desired_positions) {
    position_freq[[i]] <- setNames(rep(0, length(valid_amino_acids)), valid_amino_acids)
  }
  
  # Collapse all sequences into a single string to calculate total frequencies
  all_seqs <- paste(data$seq, collapse = "")
  total_freq <- as.data.frame(table(strsplit(all_seqs, "")[[1]]))
  colnames(total_freq) <- c("amino_acid", "total_count")
  total_freq <- total_freq %>% filter(amino_acid %in% valid_amino_acids)
  
  # Calculate the global total count of all amino acids
  global_total <- sum(total_freq$total_count)
  
  # Process each truncated sequence
  for (seq in data$seq) {
    seq_length <- min(nchar(seq), desired_positions)
    for (i in 1:seq_length) {
      amino_acid <- substr(seq, i, i)
      if (amino_acid %in% valid_amino_acids) {
        position_freq[[i]][amino_acid] <- position_freq[[i]][amino_acid] + 1
      }
    }
  }
  
  # Convert frequencies to data frame and calculate percentages
  range <- if (start_position == "first") {
    1:desired_positions
  } else if (start_position == "last") {
    (1:desired_positions - 1) * -1
  } else {stop()}
  
  position_percentages <- map_df(seq_along(range), function(pos) {
    total <- sum(position_freq[[pos]])
    freq_df <- as.data.frame(as.table(position_freq[[pos]]))
    colnames(freq_df) <- c("amino_acid", "count")
    freq_df$position <- names(position_freq)[pos]
    freq_df$percentage <- if (total > 0) (freq_df$count / total) * 100 else 0
    freq_df
  })
  
  # Prepare the total count data
  total_percentages <- total_freq %>%
    mutate(position = "total", percentage = (total_count / global_total) * 100) %>%
    rename(count = total_count)
  
  # Combine the position and total data
  combined_percentages <- bind_rows(position_percentages, total_percentages)
  
  # Order the data to interleave positions and totals
  combined_percentages <- combined_percentages %>%
    mutate(position = factor(position, levels = c(as.character(range), "total"))) %>%
    arrange(amino_acid, position) %>%
    mutate(group = cancer_type, class = class)
  
  return(combined_percentages)
}
