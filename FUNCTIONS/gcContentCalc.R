gcContentCalc <- function(sequence) {
  # Split the sequence into nucleotides
  nucleotides <- unlist(strsplit(sequence, ""))
  
  # Calculate the length of the sequence
  sequence_length <- length(nucleotides)
  
  # Calculate the frequency table of nucleotides
  nucleotide_counts <- table(nucleotides)
  
  # Calculate the percentage of G and C
  g_percentage <- ifelse("G" %in% names(nucleotide_counts), nucleotide_counts["G"] / sequence_length * 100, 0)
  c_percentage <- ifelse("C" %in% names(nucleotide_counts), nucleotide_counts["C"] / sequence_length * 100, 0)
  
  # Sum the percentages of G and C to get the GC content
  gc_content <- g_percentage + c_percentage
  
  return(gc_content)
}
