library(showtext)

showtext_auto()

font_add(family = "CMU Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gtools)
library(cowplot)

load("~/Documents/UNAV/MC2/TFM/DATA/datasets/ORFs_HCC_2024.RData")
counts_HCC_Navarra <- read_excel("DATA/counts_HCC_Navarra.xlsx")

pacientes <- names(counts_HCC_Navarra)[7:16]

counts_HCC_Navarra_normalized <- counts_HCC_Navarra %>%
  filter(complete.cases(.)) %>% 
  mutate(Length_per_kb = Length / 1000) %>%
  relocate(Length_per_kb, .after = Length) %>%
  mutate(across(all_of(pacientes), ~ .x / Length_per_kb)) %>% 
  mutate(across(all_of(pacientes), ~ .x / (sum(.x) / 1e6)) ) %>% 
  dplyr::select(Geneid, all_of(pacientes)) %>% 
  rename(gene = Geneid)

ORFs_normalized <- ORFs_HCC_2024 %>% 
  rename(
    L1_S1 = "X285_BH25MNDSX7_S1..hcc_hcc_20230212",
    L2_S2 = "X286_BH25MNDSX7_S2..hcc_hcc_20230212",
    L3_S3 = "X287_BH25MNDSX7_S3..hcc_hcc_20230212",
    L4_S4 = "X288_BH25MNDSX7_S4..hcc_hcc_20230212",
    L5_S5 = "X289_BH25MNDSX7_S5..hcc_hcc_20230212",
    L6_S6 = "X290_BH25MNDSX7_S6..hcc_hcc_20230212",
    L32_S17 = "X316_BH25MNDSX7_S30..hcc_hcc_20230212",
    L33_S18 = "X317_BH25MNDSX7_S31..hcc_hcc_20230212",
    L34_S19 = "X318_BH25MNDSX7_S32..hcc_hcc_20230212",
    L35_S20 = "X319_BH25MNDSX7_S33..hcc_hcc_20230212"
  ) %>% 
  mutate(len_per_kb = len / 1000) %>%
  relocate(len_per_kb, .after = len) %>%
  mutate(across(all_of(pacientes), ~ .x / len)) %>% 
  mutate(across(all_of(pacientes), ~ .x / (sum(.x) / 1e6)) ) %>% 
  dplyr::select(orf, gene, class, all_of(pacientes), DNA_sequence)

# Combinar tablas de psite intensity y counts
combined_table <- left_join(ORFs_normalized, counts_HCC_Navarra_normalized,
                            by = "gene",
                            suffix = c("_Psites_int", "_counts_per_nt"))

# Idenitificar las columnas que son de psites, y las que son de counts
counts_cols <- grep("_counts_per_nt$", names(combined_table), value = TRUE)
psites_cols <- sub("_counts_per_nt$", "_Psites_int", counts_cols)



filtered_data_for_plot <- combined_table %>%
  pivot_longer(
    cols = all_of(c(psites_cols, counts_cols)),
    names_to = c("patient_id", ".value"),
    names_pattern = "(.*)_(Psites_int|counts_per_nt)"
  ) %>% 
  mutate(counts_per_nt_mod = counts_per_nt + 1e-1,
         Psites_int_mod = Psites_int + 1e-1,
         translation_index_mod = Psites_int_mod / counts_per_nt_mod,
         translation_index = Psites_int / counts_per_nt) %>%
  relocate(Psites_int_mod, .after = Psites_int) %>% 
  mutate(patient_id = factor(patient_id, levels = mixedsort(unique(patient_id)))) %>%
  mutate(start_codon = ifelse(substr(DNA_sequence, 1, 3) == "ATG", "ATG", "No ATG")) %>% 
  dplyr::select(-DNA_sequence) %>% 
  filter(Psites_int > 0,
         counts_per_nt > 0) %>% 
  ungroup()

plot1 <- ggplot(filtered_data_for_plot, aes(x = class, y = translation_index)) +
  geom_violin() +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  geom_boxplot(width = 0.25, outliers = F) +
  labs(x = "Clase",
       y = "Log10(Índice de traducción)") +
  scale_x_discrete(labels = c("lncRNA", "uoORF", "uORF"))

plot2 <- ggplot(filtered_data_for_plot, aes(x = patient_id, y = translation_index)) +
  geom_violin() +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  geom_boxplot(width = 0.25, outliers = F) +
  labs(x = "ID del paciente",
       y = "Log10(Índice de traducción)")

# create a grid for the first plot to adjust its width
plot1_wrapped <- plot_grid(NULL, plot1, NULL,
                           ncol = 3,
                           rel_widths = c(0.33,0.33,0.33),
                           labels = c("", "(a)", ""),
                           label_size = 10,
                           label_fontfamily = "CMU_Serif",
                           label_x = -0.05)


plot_grid(
  plot1_wrapped, plot2,
  ncol = 1,
  labels = c("", "(b)"),
  label_fontfamily = "CMU_Serif",
  label_size = 10,
  label_x = -0.02
)

ggsave(plot = last_plot(),
       filename = "DOCUMENTS/Memoria/Images/translation_index/translation_index_combined.pdf",
       width = 16,
       height = 14,
       units = "cm",
       device = "pdf")
