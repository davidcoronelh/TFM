library(showtext)
library(dplyr)
library(ggplot2)
library(ggseqlogo)
library(cowplot)

showtext_auto()

font_add(family = "CMU_Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


load("DATA/datasets/CDSs_HCC_2024.RData")
load("DATA/datasets/ORFs_HCC_2024.RData")

# Filter the data sets to avoid redundancies in the analysis 
five_prime_analysis <- list(CDSs_HCC_2024, ORFs_HCC_2024) %>%
  lapply(function(x) {
    positive <- x %>% 
      filter(sense == "+") %>% 
      group_by(chromosome, posicion_inicial) %>% 
      arrange(desc(len)) %>%
      slice_head()

    negative <- x %>% 
      filter(sense == "-") %>% 
      group_by(chromosome, posicion_final) %>%
      arrange(desc(len)) %>%
      slice_head()
    combined <- rbind(positive, negative)
    return(combined)
  }) %>%
  do.call(rbind, .)


# Function to generate sequence logos
generateSeqLogo <- function(dataset, target_class, start_codon = "any", y_limit = 0.2, full_title = FALSE) {
  showtext_auto()
  require(dplyr)
  
  target_class <- as.character(target_class)
  start_codon <- as.character(start_codon)
  
  y_limit <- ifelse("CDS" %in% target_class, 0.4, y_limit)
  
  filtered_sequences <- dataset %>%
    filter(class %in% target_class) %>%
    pull(exp_seq) %>%
    substring(91, 113)

  if (start_codon == "non_atg") {
    filtered_sequences <- filtered_sequences[substr(filtered_sequences, 11, 13) != "ATG"]
  } else if (start_codon != "any") {
    filtered_sequences <- filtered_sequences[substr(filtered_sequences, 11, 13) %in% start_codon]
  }
  
  
  # Step 3: Transform the sequences by inserting "---" at the desired position
  kozak_region <- sapply(filtered_sequences, function(seq) {
    paste0(substr(seq, 1, 10), "---", substr(seq, 14, nchar(seq)))
  })
  
  plot_title <- 
    paste0(
      paste(target_class, collapse = ", "),
      "; ",
      " N = ", length(kozak_region),
      # if (start_codon != "any" | start_codon == "non_atg") paste0("Codón de inicio: ", start_codon, "; ") else ""
      if (start_codon == "non_atg") paste0("; Codón de inicio: No ATG") else if (start_codon != "any") paste0("; Codón de inicio: ", start_codon) else ""
      
      )
 
    plot <- ggseqlogo(kozak_region) +
    ggtitle(plot_title) +
    ylim(c(0, y_limit)) +
    scale_x_continuous(breaks = c(1:23),
                       labels = if ("CDS" %in% target_class) {
                         c(-10:-1, "A", "T", "G", 4:13) 
                       } else if (start_codon %in% c("any", "non_atg")) {
                         c(-10:-1, 1:13)
                       } else {
                         c(-10:-1, unlist(strsplit(start_codon, "")), 4:13)
                       }) +
    theme(text = element_text(family = "CMU_Serif"),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),    
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  return(plot)
}

# Generate sequence logos
plots <- list(
  generateSeqLogo(dataset = five_prime_analysis, target_class = "CDS"),
  generateSeqLogo(dataset = five_prime_analysis, target_class = c("uORF", "uoORF", "lncRNA-ORF"), y_limit = 0.4)
)

for (class in c("uORF", "uoORF", "lncRNA-ORF")) {
  plots[[length(plots) + 1]] <- generateSeqLogo(dataset = five_prime_analysis, target_class = class)
}

non_canonical_atg <- five_prime_analysis %>% 
  filter(substr(DNA_sequence, 1, 3) == "ATG" & class != "CDS")

non_canonical_ntg <- five_prime_analysis %>% 
  filter(substr(DNA_sequence, 1, 3) != "ATG" & class != "CDS")

plots[[length(plots) + 1]] <- generateSeqLogo(dataset = non_canonical_atg, target_class = c("uORF", "uoORF", "lncRNA-ORF"), start_codon = "ATG")
plots[[length(plots) + 1]] <- generateSeqLogo(dataset = non_canonical_ntg, target_class = c("uORF", "uoORF", "lncRNA-ORF"), start_codon = "non_atg")

# Combine plots
combined_plot <- plot_grid(
  plotlist = plots,
  labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"),
  ncol = 1,
  label_fontfamily = "CMU_Serif",
  label_size = 9
)

# Save plot
path <- "./DOCUMENTS/Memoria/Images/kozak_sequences.pdf"

ggsave(plot = combined_plot,
       filename = path,
       width = 12,
       height = 24,
       units = "cm")
