library(showtext)
library(dplyr)
library(ggplot2)
library(readxl)

showtext_auto()

font_add(family = "CMU Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")

load("DATA/datasets/ORFs_HCC_2024.RData")
load("DATA/datasets/CDSs_HCC_2024.RData")


codonBiasCalculator <- readRDS("FUNCTIONS/codonBiasCalculator.rds")
filterRedundantORFs <- readRDS("FUNCTIONS/filterRedundantORFs.rds")
gcContentCalc <- readRDS("FUNCTIONS/gcContentCalc.rds")


non_canonical_codon_bias <- codonBiasCalculator(ORFs_HCC_2024$cDNA)

public_codon_bias <- read_excel("DATA/codon bias.xlsx") %>%
  select(codon = 1, number = 5, amino_acid = 2, frequency_thousand = 4, fraction = 3) %>%
  mutate(across(c(frequency_thousand, number, fraction), as.numeric)) %>%
  arrange(amino_acid, codon)

# canonical_codon_bias <- codonBiasCalculator(CDSs_HCC_2024$cDNA)

## BLAND-ALTMAN PLOT
plot_bland_altman_codon_frequencies <- function(non_canonical_codon_bias, public_codon_bias) {
  require(ggrepel)
 
  y <- non_canonical_codon_bias %>%
    mutate(source = "non_canonical") %>%
    select(codon, amino_acid, frequency_thousand, source) %>% 
    mutate(codon = paste0(codon, "_", amino_acid))
  
  
  x <- public_codon_bias %>%
    mutate(source = "public_data") %>% 
    select(codon, amino_acid, frequency_thousand, source) %>%
    mutate(codon = paste0(codon, "_", amino_acid))
  
  # Combine the tables by the codon value
  combined_data <- x %>%
    left_join(y, by = c("codon", "amino_acid"), suffix = c("_x", "_y"))
  
  
  filtered_data <- combined_data %>%
    filter(!is.na(frequency_thousand_x) & !is.na(frequency_thousand_y)) %>% 
    rowwise() %>% 
    mutate(gc_content = gcContentCalc(substr(codon, 1, 3)),
           label = case_when(
             gc_content > 66.6 ~ "high_gc",
             gc_content < 33.4 ~ "low_gc",
             TRUE ~ "neutral"))
  
  # Calculate means and differences 
  filtered_data <- filtered_data %>%
    mutate(mean_bias = (frequency_thousand_x + frequency_thousand_y) / 2,
           difference = frequency_thousand_y - frequency_thousand_x)
  
  
   bland_altman_plot <- ggplot(filtered_data, aes(x = mean_bias, y = difference)) +
    geom_point(aes(fill = label), shape = 21, size = 2) + 
    geom_hline(yintercept = mean(filtered_data$difference), linetype = "dashed", color = "blue3") +
    geom_hline(yintercept = mean(filtered_data$difference) + 1.96 * sd(filtered_data$difference), linetype = "dotted", color = "red") +
    geom_hline(yintercept = mean(filtered_data$difference) - 1.96 * sd(filtered_data$difference), linetype = "dotted", color = "red") +
    geom_text(aes(x = max(mean_bias), y = mean(filtered_data$difference)), 
              label = paste("Media = ", round(mean(filtered_data$difference), 2)), 
              vjust = -1.5, hjust = 1, family = "CMU Serif", color = "blue3", size = 3.5) +
    geom_text(aes(x = max(mean_bias), y = mean(filtered_data$difference) + 1.96 * sd(filtered_data$difference)), 
              label = paste("+1.96 SD = ", round(mean(filtered_data$difference) + 1.96 * sd(filtered_data$difference), 2)), 
              vjust = -1.5, hjust = 1, family = "CMU Serif", color = "red", size = 3.5) +
    geom_text(aes(x = max(mean_bias), y = mean(filtered_data$difference) - 1.96 * sd(filtered_data$difference)), 
              label = paste("-1.96 SD = ", round(mean(filtered_data$difference) - 1.96 * sd(filtered_data$difference), 2)), 
              vjust = -1.5, hjust = 1, family = "CMU Serif", color = "red", size = 3.5) +
    geom_text_repel(aes(label = codon), size = 2.5, max.overlaps = 40, family = "CMU Serif") +
    labs(x = "Promedio de las frecuencias",
         y = "Diferencia de las frecuencias (No canónica - Canónica)",
         fill = "Tipo de codón") +
    theme_bw() +
    theme(legend.position = "top",
          text = element_text(family = "CMU Serif")) +
    scale_fill_manual(values = c("high_gc" = "indianred2", "low_gc" = "skyblue2"),
                      labels = c("high_gc" = "Alto en GC", "low_gc" = "Alto en AT"))
  
  return(bland_altman_plot)
}

plot_bland_altman_codon_frequencies(non_canonical_codon_bias, public_codon_bias)

ggsave(plot = last_plot(),
       filename = "DOCUMENTS/Memoria/Images/bland_altman.pdf",
       width = 16,
       height = 16,
       units = "cm",
       device = "pdf")
