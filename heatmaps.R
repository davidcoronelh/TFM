library(showtext)

showtext_auto()

font_add(family = "CMU_Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")

library(dplyr)
library(cowplot)

load("DATA/datasets/all_datasets.RData")
filterRedundantORFs <- readRDS("./FUNCTIONS/filterRedundantORFs.rds")
createPercentageTable <- readRDS("./FUNCTIONS/createPercentageTable.rds")
aminoAcidEnrichment <- readRDS("./FUNCTIONS/aminoAcidEnrichment.rds")


my_heatmap <- function(CDS_table, ORF_table, start_position, cancer_type) {
  
  number_of_positions <- ifelse(start_position == "first", 21, 20)
  
  CDS_table <- filterRedundantORFs(CDS_table)
  ORF_table <- filterRedundantORFs(ORF_table)
  
  CDS_info <- createPercentageTable(data = CDS_table,
                                    start_position = start_position,
                                    desired_positions = number_of_positions,
                                    cancer_type = cancer_type) 
  
  ORF_info <- createPercentageTable(data = ORF_table,
                                    start_position = start_position,
                                    desired_positions = number_of_positions,
                                    cancer_type = cancer_type)
  
  fold_change <- aminoAcidEnrichment(orf_percentage_data = ORF_info,
                                     cds_percentage_data = CDS_info)
  
  heatmap <- plotHeatmap(fold_change_data = fold_change)
  
  return(heatmap)
}

plotHeatmap <- function(fold_change_data) {
  
  require(ggplot2)
  require(gtools)
  require(dplyr)
  
  datos_filtrados <- fold_change_data %>%
    group_by(amino_acid) %>%
    mutate(promedio = mean(log2FC[!(position %in% c("total", "end"))], na.rm = TRUE),
           amino_acid = case_when(
             amino_acid == "A" ~ "Ala",
             amino_acid == "R" ~ "Arg",
             amino_acid == "N" ~ "Asn",
             amino_acid == "D" ~ "Asp",
             amino_acid == "C" ~ "Cys",
             amino_acid == "Q" ~ "Gln",
             amino_acid == "E" ~ "Glu",
             amino_acid == "G" ~ "Gly",
             amino_acid == "H" ~ "His",
             amino_acid == "I" ~ "Ile",
             amino_acid == "L" ~ "Leu",
             amino_acid == "K" ~ "Lys",
             amino_acid == "M" ~ "Met",
             amino_acid == "F" ~ "Phe",
             amino_acid == "P" ~ "Pro",
             amino_acid == "S" ~ "Ser",
             amino_acid == "T" ~ "Thr",
             amino_acid == "W" ~ "Trp",
             amino_acid == "Y" ~ "Tyr",
             amino_acid == "V" ~ "Val",
             TRUE ~ NA_character_
           )
           ) %>%
    arrange(desc(promedio)) %>% 
    mutate(position = as.character(position),
           position = ifelse(position == "total", "T", position))
  

  
  positions <- mixedsort(unique(datos_filtrados$position))

  
  if (unique(fold_change_data$group) == "BRE") {
    
    titulo <- "BRE"
    } else if (unique(fold_change_data$group) == "HCC") {
      titulo <- "HCC"
    } else {
      titulo <- "HCC 2024"
    }
  
  # Create the heatmap with the ordered data
  heatmap <- ggplot(datos_filtrados, aes(
    x = factor(position, levels = positions),
    y = factor(amino_acid, levels = rev(unique(amino_acid))),
    fill = log2FC)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac",
                         mid = "white",
                         high = "#b2182b",
                         limits = c(-1.75, 1.75),
                         midpoint = 0) +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = titulo
    ) +
    theme(
      text = element_text(family = "CMU_Serif"),
      axis.text = element_text(size = 5),
      axis.title  = element_blank(),
      title = element_text(size = 6)
    )
  
  return(heatmap)
}

# List of plots without the legends
plots <- list(
  my_heatmap(CDSs_BRE, ORFs_BRE, "first", "BRE") + theme(legend.position = "none"),
  my_heatmap(CDSs_BRE, ORFs_BRE, "last", "BRE") + theme(legend.position = "none"),
  my_heatmap(CDSs_HCC, ORFs_HCC, "first", "HCC") + theme(legend.position = "none"),
  my_heatmap(CDSs_HCC, ORFs_HCC, "last", "HCC") + theme(legend.position = "none"),
  my_heatmap(CDSs_HCC_2024, ORFs_HCC_2024, "first", "HCC_2024") + theme(legend.position = "none"),
  my_heatmap(CDSs_HCC_2024, ORFs_HCC_2024, "last", "HCC_2024") + theme(legend.position = "none")
)

# Create one plot with the legend to extract it
legend_plot <- my_heatmap(CDSs_BRE, ORFs_BRE, "last", "BRE")

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

# Combine all plots 
combined_plot <- cowplot::plot_grid(
  plotlist = plots, 
  ncol = 2, 
  labels = c("(a)", "", "", "", "", ""),
  label_size = 9, 
  label_x = -0.1, 
  label_y = 0.99, 
  label_fontfamily = "CMU_Serif"
)


final_plot_with_labels <- cowplot::ggdraw() +
  cowplot::draw_plot(combined_plot, x = 0.05, y = 0.05,  width = 0.95, height =  0.95) +  # Adjust plot position and size to add padding
  cowplot::draw_label("Posición desde el inicio", x = 0.3, y = 0.05, vjust = 1, size = 9, fontfamily = "CMU_Serif") +  # X-axis label for first column
  cowplot::draw_label("Posición desde el final", x = 0.77, y = 0.05, vjust = 1, size = 9, fontfamily = "CMU_Serif") +  # X-axis label for second column
  cowplot::draw_label("Aminoácido", x = 0.05, y = 0.5, angle = 90,  size = 10, fontfamily = "CMU_Serif")  # Common Y-axis label with padding

# Combine the plot with the legend
final_plot <- cowplot::plot_grid(
  final_plot_with_labels, legend, ncol = 2, rel_widths = c(1, 0.1)
)


# Save plot
path <- "./DOCUMENTS/Memoria/Images/heatmaps/all_heatmaps.pdf"
ggsave(plot = final_plot,
       filename = path,
       width = 16,
       height = 20,
       units = "cm",
       device = "pdf")
