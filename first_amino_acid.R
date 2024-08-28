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


gcContentCalc <- readRDS("FUNCTIONS/gcContentCalc.rds")
filterRedundantORFs <- readRDS("FUNCTIONS/filterRedundantORFs.rds")

# Filtrar los ORFs redundantes
filtered_datasets <- lapply(datasets_list, filterRedundantORFs)

all_filtered_ORFs <- bind_rows(filtered_datasets, .id = "group") %>% 
  mutate(group = case_when(
    grepl("BRE", group) ~ "BRE",
    grepl("HCC_2024", group) ~ "HCC_2024",
    grepl("HCC", group) ~ "HCC",
    TRUE ~ group
  )) %>%
  mutate(first_aa = substr(seq, 1, 1)) %>%
  mutate(class = ifelse(class == "CDS", "CDS", "non_CDS")) %>%
  group_by(group, class) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(first_aa, group, class) %>%
  summarize(n = n(), percentage = (n / first(total))) %>%
  ungroup() %>%
  filter(n > 4)

group_labels <- c("BRE" = "BRE",
                  "HCC" = "HCC",
                  "HCC_2024" = "HCC 2024")

legend_lables <- c("CDS" = "CDS", "non_CDS" = "No CDS")

yticks <- c(0, 0.10, 0.20, 0.30, 0.40, 1.00)

ggplot(all_filtered_ORFs, aes(x = first_aa, y = percentage, fill = class)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw() +
  facet_grid(~group,
             labeller = labeller(group = group_labels)) +
  scale_fill_grey(labels = legend_lables) +
  labs(x = "(b) Primer aminoácido", y = "Frecuencia", fill = "Clase") +
  scale_y_continuous(limits = c(0, 1.1), breaks = yticks) +
  theme(
    legend.position = "none",
    text = element_text(family = "CMU Serif"),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 8)
  ) +
  scale_y_cut(breaks=c(0.4), which = c(1,2), scales = c(0.1, 1), space = 0.25)


ggsave(plot = last_plot(),
       filename = "DOCUMENTS/Memoria/Images/first_amino_acid.pdf",
       width = 16,
       height = 7,
       units = "cm",
       device = "pdf")