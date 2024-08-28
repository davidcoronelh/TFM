library(showtext)

showtext_auto()

font_add(family = "CMU_Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(dplyr)
library(cowplot)
library(gtools)

# Cargar los datos de los pacientes
load("~/Documents/UNAV/MC2/TFM/DATA/datasets/all_datasets.RData")

BRE = bind_rows(CDSs_BRE, ORFs_BRE)
HCC = bind_rows(CDSs_HCC, ORFs_HCC)
HCC_2024 = bind_rows(CDSs_HCC_2024, ORFs_HCC_2024)

createPercentageTable <- readRDS("FUNCTIONS/createPercentageTable.rds")
filterRedundatORFs <- readRDS("FUNCTIONS/filterRedundantORFs.rds")

combined_plot <- function(amino_acid, ylims, return_legend = FALSE) {
  
  require(gtools)
  individual_amino_acid_enrichment <- function(aa, cancer_type, position, orf_data, cds_data, ylims) { 
    
    number_of_positions <- ifelse(position == "first", 21, 20)
    
    amino_acid_ORF <- createPercentageTable(data = orf_data,
                                            start_position = position,
                                            desired_positions = number_of_positions,
                                            cancer_type = cancer_type) %>%
      filter(amino_acid == aa) %>%
      mutate(class = "ORF")
    
    amino_acid_CDS <- createPercentageTable(data = cds_data,
                                            start_position = position,
                                            desired_positions = number_of_positions,
                                            cancer_type = cancer_type) %>%
      filter(amino_acid == aa) %>%
      mutate(class = "CDS")
    
    amino_acid_data <- bind_rows(amino_acid_CDS, amino_acid_ORF) %>% 
      filter(position != 1) %>%
      mutate(position = as.character(position),
             position = ifelse(position == "total", "T", position))
    
    
    amino_acid_data$position <- factor(
      amino_acid_data$position, 
      levels = unique(mixedsort(amino_acid_data$position))
    )
    
    amino_acid_dict <- list(
      A = "Alanina",
      R = "Arginina",
      N = "Asparagina",
      D = "Aspartato",
      C = "Cisteína",
      E = "Glutamato",
      Q = "Glutamina",
      G = "Glicina",
      H = "Histidina",
      I = "Isoleucina",
      L = "Leucina",
      K = "Lisina",
      M = "Metionina",
      F = "Fenilalanina",
      P = "Prolina",
      S = "Serina",
      T = "Treonina",
      W = "Triptófano",
      Y = "Tirosina",
      V = "Valina"
    )
    
    aa_label <- amino_acid_dict[[aa]]
    
    p <- ggplot(amino_acid_data, aes(x = position, y = percentage, fill = class)) +
      geom_bar(stat = "identity", position = position_dodge(), width = 1.5) +
      facet_wrap(~ position, scales = "free_x", ncol = 21) +
      scale_fill_grey(name = "Clase") +
      theme(
        panel.grid = element_blank(),  
        strip.text = element_blank(),  
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),  
        plot.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        text = element_text(family = "CMU_Serif"),
        panel.background = element_blank(),
        legend.position = "none"
      ) +
      ggtitle(if (position == "first") aa_label else "") +
      scale_y_continuous(limits = ylims)
    
    return(p)
  }
  
  plots <- list()
  i <- 1
  
  for (aa in amino_acid) { 
    for (pos in c("first", "last")) {
      plots[[i]] <- individual_amino_acid_enrichment(
        aa = aa,
        cancer_type = "HCC",
        position = pos,
        orf_data = ORFs_HCC,
        cds_data = CDSs_HCC,
        ylims = ylims
      )
      i <- i + 1
    }
  }
  
  get_legend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  legend <- get_legend(plots[[1]] +
                         theme(legend.position = "top",
                               text = element_text(family = "CMU_Serif", size = 12)))
  
  combined_plot <- grid.arrange(
    do.call(arrangeGrob, c(plots[1:2], ncol = 1)),  # Combine only the first two plots
    bottom = textGrob("Posición", gp = gpar(fontsize = 8, fontfamily = "CMU_Serif")),  # X-axis label
    left = textGrob("Porcentaje", gp = gpar(fontsize = 8, fontfamily = "CMU_Serif"), rot = 90),  # Y-axis label
    nrow = 2, 
    heights = c(4, 0.05)
  )
  
  if(return_legend){return(legend)}
  
  return(combined_plot)
}


R <- combined_plot("R", ylims = c(0,15))
P <- combined_plot("P", ylims = c(0,10))
A <- combined_plot("A", ylims = c(0,25))
L <- combined_plot("L", ylims = c(0,15))

legend <- combined_plot("L", ylims = c(0,15), return_legend = T)

plot_grid(
  legend,
  plot_grid(
    R, P, A, L,
    ncol = 2,
    labels = c("(b)", "", "", ""),
    label_fontfamily = "CMU_Serif",
    label_size = 11,
    label_x = -0.01
  ),
  ncol = 1,
  rel_heights = c(0.05, 1))

ggsave(last_plot(),
       file = "./DOCUMENTS/Memoria/Images/individual_aa/combined_aa.pdf",
       width = 8.21,
       height = 6,
       device = "pdf")
