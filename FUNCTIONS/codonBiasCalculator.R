codonBiasCalculator <- function(transcript) {
  
  genetic_code <- data.frame(codon = c("TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG",
                                       "ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG",
                                       "TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG",
                                       "ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG",
                                       "TAT","TAC","TAA","TAG","CAT","CAC","CAA","CAG",
                                       "AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG",
                                       "TGT","TGC","TGA","TGG","CGT","CGC","CGA","CGG",
                                       "AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG"),
                             amino_acid = c("F","F","L","L","L","L","L","L",
                                            "I","I","I","M","V","V","V","V",
                                            "S","S","S","S","P","P","P","P",
                                            "T","T","T","T","A","A","A","A",
                                            "Y","Y","*","*","H","H","Q","Q",
                                            "N","N","K","K","D","D","E","E",
                                            "C","C","*","W","R","R","R","R",
                                            "S","S","R","R","G","G","G","G"))
  
  # concatenate the sequences
  concatenated_seqs <- paste(transcript, collapse = "")
  

  nucleotides <- nchar(concatenated_seqs)
  
  # vector to store codons
  codon <- vector("character", ceiling(nucleotides / 3))
  
  # extract codons
  for (i in seq(1, nucleotides, by = 3)) {
    codon[(i - 1) / 3 + 1] <- substring(concatenated_seqs, i, min(i + 2, nucleotides))
  }
  
 
  codon_table <- table(codon)
  
  codon_count <- as.data.frame(codon_table)
  colnames(codon_count) <- c("codon", "number")
  
  # merge with the genetic code
  bias <- merge(codon_count, genetic_code, by = "codon")
  
  # calculate the fraction and frequency per thousand
  bias <- within(bias, {
    fraction <- ave(number, amino_acid, FUN = function(x) x / sum(x))
    total_number <- sum(number)
    frequency_thousand <- number / total_number * 1000
  })
  
  # order by amino acid
  bias <- bias[order(bias$amino_acid), ]
  
  return(bias)
}
