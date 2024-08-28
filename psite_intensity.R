library(showtext)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(rstatix)

showtext_auto()

font_add(family = "CMU Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")

load("~/Documents/UNAV/MC2/TFM/DATA/datasets/all_datasets.RData")


process_data <- function(dataset) {
  dataset %>%
    mutate(class = as.character(class)) %>%
    mutate(class = ifelse(class %in% c("lncRNA-ORF", "ncRNA-ORF"), 
                          sub("-ORF", "", class), class))
}

BRE <- bind_rows(CDSs_BRE, ORFs_BRE) %>% process_data()
HCC <- bind_rows(CDSs_HCC, ORFs_HCC) %>% process_data()
HCC_2024 <- bind_rows(CDSs_HCC_2024, ORFs_HCC_2024) %>% process_data()


set_class_factor <- function(data, levels) {
  data$class <- factor(data$class, levels = levels)
  data
}

BRE <- set_class_factor(BRE, c("CDS", "uORF", "uoORF", "intORF", "doORF", "dORF", "lncRNA", "ncRNA"))
HCC <- set_class_factor(HCC, c("CDS", "uORF", "uoORF", "intORF", "doORF", "dORF", "lncRNA", "ncRNA"))
HCC_2024 <- set_class_factor(HCC_2024, c("CDS", "uORF", "uoORF", "lncRNA"))

# Function to plot p-site intensity
plot_psite_int <- function(table, title) {
  ggplot(table, aes(x = class, y = mean_norm_psite_int_per_aa)) +
    geom_boxplot(width = 0.4, outlier.alpha = 0.5, outlier.size = 0.5) +
    theme_bw() +
    theme(text = element_text(family = "CMU Serif"),
          axis.text.x = element_text(angle = 45, vjust = 0.5),
    ) +
    scale_y_log10(limits = c(1e-4, 1e+4)) +
    ggtitle(title)
}

# Function to compare expression
comparar_expresion <- function(tabla, title) {
  plot_data_x <- tabla %>% 
    filter(class == "CDS") %>% 
    group_by(gene) %>% 
    arrange(desc(mean_norm_psite_int_per_aa)) %>% 
    slice_head(n = 1) %>% 
    dplyr::select(class, orf, gene, mean_norm_psite_int_per_aa)
  
  plot_data_y <- tabla %>% 
    filter(class == "uORF") %>% 
    group_by(gene) %>% 
    arrange(desc(mean_norm_psite_int_per_aa)) %>% 
    slice_head(n = 1) %>% 
    dplyr::select(class, orf, gene, mean_norm_psite_int_per_aa)
  
  final_plot_data <- inner_join(plot_data_x, plot_data_y, by = "gene")
  
  pearson_coeff <- cor(final_plot_data$mean_norm_psite_int_per_aa.x,
                       final_plot_data$mean_norm_psite_int_per_aa.y,
                       method = "pearson")
  
  spearman_coeff <- cor(final_plot_data$mean_norm_psite_int_per_aa.x,
                        final_plot_data$mean_norm_psite_int_per_aa.y,
                        method = "spearman")
  
  
  p <- ggplot(final_plot_data, aes(x = mean_norm_psite_int_per_aa.x, y = mean_norm_psite_int_per_aa.y)) +
    geom_point(size = 0.1, alpha = 0.4) +
    geom_smooth(method = "lm", color = "firebrick3", linetype = "dashed", linewidth = 0.5) +
    theme_bw() +
    theme(
      text = element_text(family = "CMU Serif"),
      axis.text.x = element_text(angle = 45, hjust = 1),
    ) +
    scale_y_log10(limits = c(1e-4, 1e4)) + 
    scale_x_log10(limits = c(1e-3, 1e4)) +
    geom_text(aes(label = sprintf("r: %.2f", pearson_coeff), x = Inf, y = 0), 
              hjust = 1.1, vjust = -1, family = "CMU Serif") +
    geom_text(aes(label = paste0("ρ: ", sprintf("%.2f", spearman_coeff)), x = Inf, y = 0),
              hjust = 1.1, vjust = -4, family = "CMU Serif") +
    ggtitle(title)
  
  if (title == "BRE") {
    p <- p + 
      xlab("") +
      ylab("Int. del sitio P (uORF)")
  } else if (title == "HCC") {
    p <- p + 
      xlab("Int. del sitio P (CDS)") + 
      ylab("")
  } else {
    p <- p +
      xlab("") +
      ylab("")
  }
    return(p)
}


# Generate the individual plots
p1 <- plot_psite_int(BRE, "BRE") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

p2 <- plot_psite_int(HCC, "HCC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p3 <- plot_psite_int(HCC_2024, "HCC 2024") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

row_1 <- plot_grid(
  p1, p2, p3,
  nrow = 1,
  rel_widths = c(1, 0.87, 0.87))


p4 <- comparar_expresion(BRE, "BRE") +
  theme(
    axis.title.x = element_blank())

p5 <- comparar_expresion(HCC, "HCC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p6 <- comparar_expresion(HCC_2024, "HCC 2024") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

row_2 <- plot_grid(
  p4, p5, p6,
  nrow = 1,
  rel_widths = c(1, 0.8, 0.8),
  axis = 'l')

combined_plot <- plot_grid(
  row_1, row_2,
  nrow = 2,
  labels = c("(a)", "(b)"),
  label_fontfamily = "CMU Serif",
  label_size = 10)


print(combined_plot)


# ggsave(plot = combined_plot,
#        filename = "./DOCUMENTS/Memoria/Images/psite_intensity/psite_intensity.png",
#        device = "png",
#        width = 604,
#        height = 500,
#        units = "px",
#        dpi = 300
# )

## ANOVA DE WELCH ---

HCC <- HCC %>%
  # mutate(group_interaction = interaction(class, group)) %>% 
  #  rowwise()) %>% 
  mutate(log_10_ps = log10(mean_norm_psite_int_per_aa)) %>% 
  ungroup()

# Perform Welch's ANOVA ----
welch_anova <- HCC %>%
  welch_anova_test(log_10_ps ~ class)

print(welch_anova)

# Perform the Games-Howell post-hoc test ----
games_howell_results <- HCC %>% 
  games_howell_test(log_10_ps ~ class)

print(games_howell_results, n = 28)
