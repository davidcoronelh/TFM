aminoAcidEnrichment <- function(orf_percentage_data, cds_percentage_data) {
  
  require(dplyr)
  
  if (unique(cds_percentage_data$group) != unique(orf_percentage_data$group)) {
    stop("The tables belong to different groups")
  }
  
  calc_log2FC <- function(n) {
    if (is.nan(n) || n == 0) {
      return(0)
    }
    log2_fold_change <- log2(n)
    return(log2_fold_change)
  }
  
  # Subset the data to exclude position 1
  orf_subset <- subset(orf_percentage_data, position != 1)
  cds_subset <- subset(cds_percentage_data, position != 1)
  
  # Ensure the data is aligned by amino acid and position before calculating the difference
  merged_data <- merge(orf_subset, cds_subset, by = c("amino_acid", "position", "group"), suffixes = c("_orf", "_cds"))
  
  result <- merged_data %>%
    select(-class_orf, -class_cds) %>% 
    mutate(perc_diff = percentage_orf / percentage_cds,
           log2FC = sapply(perc_diff, calc_log2FC))
                       
  return(result)
}
