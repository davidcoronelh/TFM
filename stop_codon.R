library(showtext)
library(dplyr)
library(ggplot2)

showtext_auto()

font_add(family = "CMU Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


load("DATA/datasets/ORFs_HCC_2024.RData")
load("DATA/datasets/CDSs_HCC_2024.RData")
gcContentCalc <- readRDS("FUNCTIONS/gcContentCalc.rds")


three_prime_analysis <- list(CDSs_HCC_2024, ORFs_HCC_2024) %>% 
  lapply(function(x) {

    positive <- x[x$sense == "+", ]
    negative <- x[x$sense == "-", ]
    

    positive_filtrados <- positive[!duplicated(positive$posicion_final), ]
    negative_filtrados <- negative[!duplicated(negative$posicion_inicial), ]
    
    combined <- rbind(positive_filtrados, negative_filtrados)
    return(combined)
  }) %>%
  do.call(rbind, .) %>% 
  rowwise() %>% 
  mutate(three_prime_utr = gcContentCalc(substr(exp_seq,
                                                nchar(exp_seq) - 100,
                                                nchar(exp_seq)))
  ) %>% 
  ungroup()

stop_codon_frequency <- three_prime_analysis %>%
  mutate(stop_codon = substr(cDNA,
                             nchar(cDNA) - 2,
                             nchar(cDNA)),
         class = ifelse(class == "CDS", "CDS", "non-CDS")) %>%
  group_by(class, stop_codon) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  group_by(class) %>% 
  mutate(total = sum(count),
         percentage = count / total * 100) %>% 
  select(class, stop_codon, count, percentage) %>% 
  filter(percentage > 0.1)


custom_colors <- c("CDS" = "white", "non-CDS" = "#696969")

labels <- c("CDS" = "CDS", "non-CDS" = "No CDS")

ggplot(stop_codon_frequency, aes(x = stop_codon,
                                 y = percentage / 100,
                                 fill = class)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           width =  0.5) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(family = "CMU Serif")) +
  labs(x = "Codón de terminación",
       y = "Frecuencia",
       fill = "Clase") +
  scale_fill_grey(labels = labels) +
  ylim(0, 1)

ggsave(plot = last_plot(),
       filename = "DOCUMENTS/Memoria/Images/stop_codons.pdf",
       width = 8,
       height = 8,
       units = "cm",
       device = "pdf")

