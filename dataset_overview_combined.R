library(showtext)
library(dplyr)
library(ggplot2)
library(cowplot)


showtext_auto()
font_add(family = "CMU_Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


load("DATA/datasets/all_datasets.RData")


BRE <- bind_rows(CDSs_BRE, ORFs_BRE) %>% mutate(group = "BRE")
HCC <- bind_rows(CDSs_HCC, ORFs_HCC) %>% mutate(group = "HCC")
HCC_2024 <- bind_rows(CDSs_HCC_2024, ORFs_HCC_2024) %>% mutate(group = "HCC_2024")
all_ORFs <- bind_rows(BRE, HCC, HCC_2024)


tabla_resumen <- all_ORFs %>% 
  mutate(class = recode(class, `lncRNA-ORF` = "lncRNA", `ncRNA-ORF` = "ncRNA")) %>% 
  group_by(group, class) %>% 
  summarise(n = n(), .groups = 'drop') %>% 
  mutate(class = factor(class, levels = c("CDS", "uORF", "uoORF", "intORF", "doORF", "dORF", "lncRNA", "ncRNA")),
         combined = interaction(class, group))


custom_colors <- c("BRE" = "#FFFFFF", "HCC" = "#A9A9A9", "HCC_2024" = "slategray2")
group_labels <- c("BRE" = "BRE", "HCC" = "HCC", "HCC_2024" = "HCC 2024")

# Plot 1: Esquema
open_reading_frames <- data.frame(
  orf = c("uORF", "uoORF", "CDS", "intORF", "doORF", "dORF", "lncRNA", "ncRNA"),
  y = c(1, 1, 1.5, 1, 1, 1, 1, 1),
  inicio = c(1, 6, 7.5, 11, 16, 21, 30, 35),
  fin = c(4, 9, 17.5, 14, 19, 24, 33, 38)
)
open_reading_frames$midpoint <- (open_reading_frames$inicio + open_reading_frames$fin) / 2

esquema <- ggplot(open_reading_frames, aes(x = inicio, y = y, xend = fin, yend = y, color = orf)) +
  annotate("segment", x = 1, xend = 24, y = 1.5, yend = 1.5, color = "black", linewidth = 0.2) +
  geom_segment(linewidth = 3) +
  geom_text(aes(x = midpoint, y = y - 0.2, label = orf), color = "black", size = 2.5, family = "CMU_Serif") +
  geom_text(aes(x = 0, y = 1.5, label = "5'"), size = 3, color = "black", family = "CMU_Serif") +
  geom_text(aes(x = 40, y = 1.5, label = "3'"), size = 3, color = "black", family = "CMU_Serif") +
  ylim(0.7, 2) + xlim(0, 40) +
  theme_minimal() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("uORF" = "#FF9999", "uoORF" = "#99CCFF", "CDS" = "#99FF99", "intORF" = "#FFCC99", "doORF" = "#FF99CC", "dORF" = "#CCCCFF", "lncRNA" = "#FFD700", "ncRNA" = "#8A2BE2"))

# Plot 2: Resumen
resumen <- ggplot(tabla_resumen, aes(x = combined, y = n)) + 
  geom_col(width = 0.5) +
  geom_text(aes(label = n), vjust = -0.5, size = 2.5, family = "CMU_Serif") +
  facet_grid(~group, scales = "free_x", space = "free_x", labeller = labeller(group = group_labels)) +
  theme_bw() +
  ylim(0, 31000) +
  labs(x = "Clase de marco abierto de lectura", y = "Número") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), text = element_text(family = "CMU_Serif")) +
  scale_x_discrete(drop = TRUE, labels = function(x) sub("\\..*", "", x))

# Plot 3: Length distribution
length_plot <- ggplot(all_ORFs, aes(x = group, y = len, fill = group)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  facet_grid(~class, scales = "free_x", space = "free_x", switch = "x", labeller = labeller(class = group_labels)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(angle = 35, size = 8),
        panel.spacing = unit(0.5, "lines"),
        text = element_text("CMU_Serif"),
        legend.position = "top",
        legend.text = element_text(size = 8)
  ) +
  scale_fill_manual(values = custom_colors, labels = group_labels) +
  labs(fill = "Grupo", y = "Log10(Longitud)") +
  scale_y_log10()

# Combine plots
combined_plot <- plot_grid(
  esquema, resumen, length_plot,
  ncol = 1,
  rel_heights = c(0.15, 0.5, 0.5),
  labels = c("(a)", "(b)", "(c)"),
  label_fontfamily = "CMU_Serif",
  label_size = 10
)

# Save plot
path <- "./DOCUMENTS/Memoria/Images/datasets_overview_combined.pdf"
  
ggsave2(plot = combined_plot,
        filename = path,
        width = 16,
        height = 16,
        units = "cm",
        device = "pdf")
