library(showtext)
library(dplyr)
library(ggplot2)
library(ggbreak)

showtext_auto()

font_add(family = "CMU Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


load("DATA/datasets/datasets_list.RData")

filterRedundantORFs <- readRDS("FUNCTIONS/filterRedundantORFs.rds")

filtered_datasets <- lapply(datasets_list, filterRedundantORFs)

all_filtered_ORFs <- bind_rows(filtered_datasets, .id = "group") %>% 
  mutate(group = case_when(
    grepl("BRE", group) ~ "BRE",
    grepl("HCC_2024", group) ~ "HCC_2024",
    grepl("HCC", group) ~ "HCC",
    TRUE ~ group
  ))


yticks <- c(0, 0.10, 0.20, 0.30, 0.40, 1.00)

start_codon_analysis <- all_filtered_ORFs %>%
  mutate(
    start_codon = substr(DNA_sequence, 1, 3),
    amino_acid = substr(seq, 1, 1),
    class = ifelse(class == "CDS", "CDS", "non_CDS")
  ) %>%
  group_by(group, class) %>%
  mutate(total_start_codons = n()) %>%
  group_by(group, class, start_codon, amino_acid) %>%
  summarise(start_codon_frequency = n() / dplyr::first(total_start_codons) * 100) %>%
  ungroup() %>% 
  filter(start_codon_frequency > 0.12)

group_labels <- c("BRE" = "BRE",
                  "HCC" = "HCC",
                  "HCC_2024" = "HCC 2024")

legend_lables <- c("CDS" = "CDS", "non_CDS" = "No CDS")


start_codon_analysis <- start_codon_analysis %>%
  mutate(interaction_label = paste(amino_acid, start_codon, sep = "_"))


ggplot(start_codon_analysis, aes(x = interaction_label, y = start_codon_frequency / 100, fill = class)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw() +
  facet_grid(~group,
             labeller = labeller(group = group_labels),
             ) +
  scale_fill_grey(labels = legend_lables) +
  labs(x = "(c) Codón de inicio", y = "Frecuencia", fill = "Clase") +
  scale_y_continuous(limits = c(0, 1.1), breaks = yticks) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    text = element_text(family = "CMU Serif"),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 8)
  ) +
  scale_y_cut(breaks=c(0.4), which = c(1,2), scales = c(0.1, 1), space = 0.25)

ggsave(plot = last_plot(),
       filename = "DOCUMENTS/Memoria/Images/start_codon.pdf",
       width = 16,
       height = 8,
       units = "cm",
       device = "pdf")



