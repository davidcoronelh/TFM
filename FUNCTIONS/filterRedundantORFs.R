filterRedundantORFs <- function(table, initial_sequence_len = 5) {
  filtered_table <- table %>%
    mutate(initial_sequence = substr(seq, 1, initial_sequence_len)) %>%
    group_by(gene) %>%
    arrange(-nchar(seq)) %>%
    group_by(initial_sequence, gene, class) %>%
    slice_head() %>%
    ungroup() %>% 
    select(-initial_sequence)
  
  return(filtered_table)
}
