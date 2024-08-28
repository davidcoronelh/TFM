library(showtext)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)  


showtext_auto()

font_add(family = "CMU Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


load("DATA/datasets/ORFs_HCC_2024.RData")
load("JULY/lncRNA/lncRNA_transcripts.RData")
gcContentCalc <- readRDS("FUNCTIONS/gcContentCalc.rds")

# Prepare data for plot 1
lncRNA_HCC_2024 <- ORFs_HCC_2024 %>%
  filter(class == "lncRNA-ORF") %>%
  rename(open_reading_frame = cDNA) %>%
  left_join(lncRNA_transcripts, by = "orf") %>%
  rowwise() %>%
  mutate(
    gc_content_orf = gcContentCalc(open_reading_frame),
    gc_content_cdna = gcContentCalc(cdna),
    first_codon = ifelse(substr(open_reading_frame, 1, 3) == "ATG", "ATG", "non_ATG"),
    codon_type = case_when(
      substr(open_reading_frame, 1, 3) == "ATG" ~ "ATG",
      substr(open_reading_frame, 1, 3) %in% c("ACG", "AGG", "CTG", "GTG") ~ "gc_rich",
      TRUE ~ "at_rich"
    )
  )

sans_atg <- lncRNA_HCC_2024 %>% 
  filter(codon_type != "ATG") %>% 
  dplyr::select(orf, codon_type, gc_content_orf, gc_content_cdna, open_reading_frame, cdna) %>% 
  drop_na() %>% 
  pivot_longer(cols = c(gc_content_orf, gc_content_cdna),
               names_to = "region",
               values_to = "gc_content")

region_labels <- c(
  "gc_content_orf" = "Open Reading Frame",
  "gc_content_cdna" = "cDNA"
)

# Calculate n values for the first plot
n_plot1 <- sans_atg %>% 
  group_by(codon_type, region) %>%
  summarise(n = n())

########### CALCULAR LA DIFERENCIA ESTADISTICA CON UN TEST DE PERMUTACION
permutation_test <- function(data, observations_col, groups_col, groupA, groupB, P = 100000, seed = 1234) {
  
  observations <- data[[observations_col]]
  groups <- data[[groups_col]]
  
  
  data <- data[groups %in% c(groupA, groupB), ]
  observations <- data[[observations_col]]
  groups <- data[[groups_col]]
  
  
  group_means <- with(data, tapply(observations, groups, mean))
  
  
 
  test.stat1 <- abs(diff(group_means))
  
  
  # PERMUTATION TEST
  set.seed(seed)
  n <- length(observations)
  
  PermSamples <- matrix(0, nrow = n, ncol = P)
  
  
  PermSamples <- replicate(P, sample(observations, size = n, replace = FALSE))
  
  
  
  Perm.test.stat1 <- apply(PermSamples, 2, function(column) {
    abs(mean(column[data[[groups_col]] == groupA]) - median(column[data[[groups_col]] == groupB]))
  })
  
  
  
  p_value1 <- mean(Perm.test.stat1 >= test.stat1)
  
  if("region" %in% names(data)) {
    return(data.frame(
      region = unique(data$region),
      pvalue_mean = p_value1
    ))
  } else {
    return(data.frame(
      pvalue_mean = p_value1
    ))
  }
}

perm_results_df_plot1 <- data.frame()

regions <- c("gc_content_cdna", "gc_content_orf")

for (re in regions) {
  
  data_permutation_1 <- sans_atg %>% filter(region == re)
  
  perm_result <- permutation_test(
    data = data_permutation_1,
    observations_col = "gc_content",
    groups_col = "codon_type",
    groupA = "at_rich",
    groupB = "gc_rich"
  )
  
  perm_results_df_plot1 <- rbind(perm_results_df, perm_result)
}

p_value_labels_plot1 <- perm_results_df_plot1 %>% 
  rowwise() %>% 
  mutate(pvalue = ifelse(pvalue_mean < 0.0001, "< 0.0001", paste("=", as.character( round( pvalue_mean, 4) ) ) ))

pval_y_position <- 12  

plot1 <- ggplot(sans_atg, aes(x = codon_type, y = gc_content)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  labs(x = "Tipo de codón de inicio", y = "Porcentaje de GC") +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_wrap(~ region, labeller = labeller(region = region_labels)) +
  ylim(0, 100) +
  geom_text(data = n_plot1, aes(x = codon_type, y = 95, label =n),
            position = position_dodge(width = 0.75), family = "CMU Serif", size = 4.5) +
  geom_text(data = p_value_labels_plot1,
            aes(x = 1.5, y = pval_y_position, label = paste("p", pvalue)),
            family = "CMU Serif", size = 3.5) +
  scale_x_discrete(labels = c("at_rich" = "Alto en AT", "gc_rich" = "Alto en GC"))


# Prepare data for plot 2
lncRNA_HCC_2024 <- lncRNA_HCC_2024 %>%
  mutate(cdna_len = nchar(cdna),
         start_codon = ifelse(substr(open_reading_frame, 1, 3) == "ATG", "ATG", "non-ATG")) %>%
  rowwise() %>%
  mutate(gc_last_100_cdna = gcContentCalc(substr(cdna, nchar(cdna) - 99, nchar(cdna)))) %>%
  ungroup() %>%
  na.omit()

n <- lncRNA_HCC_2024 %>%
  group_by(start_codon) %>%
  summarise(n = n()) %>%
  ungroup()


perm_result_plot2 <- permutation_test(
  data = lncRNA_HCC_2024,
  observations_col = "gc_last_100_cdna",
  groups_col = "first_codon",
  groupA = "ATG",
  groupB = "non_ATG")


plot2 <- ggplot(lncRNA_HCC_2024, aes(x = start_codon, y = gc_last_100_cdna)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  labs(x = "Codón de inicio", y = "Porcentaje de GC") +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_x_discrete(labels = c("ATG", "NTG")) +
  ylim(0, 100) +
  geom_text(data = n, aes(x = start_codon, y = 90, label = paste(n)),
            position = position_dodge(width = 0.75), family = "CMU Serif", size = 4.5) +
  geom_text(data = perm_result_plot2,
            aes(x = 1.5, y = pval_y_position, label = paste("p =", pvalue_mean)),
            family = "CMU Serif", size = 3.5)

# Prepare data for plot 3
lncRNA_HCC_2024 <- lncRNA_HCC_2024 %>%
  rowwise() %>%
  mutate(
    five_prime_utr = gcContentCalc(substr(exp_seq, 1, 100)),
    first_100 = gcContentCalc(substr(open_reading_frame, 1, 100)),
    three_prime_utr = gcContentCalc(substr(exp_seq, nchar(exp_seq) - 99, nchar(exp_seq))),
    last_100_cdna = gcContentCalc(substr(open_reading_frame, nchar(open_reading_frame) - 99, nchar(open_reading_frame)))
  )

q1 <- quantile(lncRNA_HCC_2024$len, 0.33, na.rm = TRUE)
q2 <- quantile(lncRNA_HCC_2024$len, 0.66, na.rm = TRUE)

lncRNA_HCC_2024 <- lncRNA_HCC_2024 %>%
  mutate(size = case_when(
    len > q2 ~ "large",
    len < q1 ~ "small",
    TRUE ~ "medium"
  )) %>%
  pivot_longer(cols = c("five_prime_utr", "first_100", "three_prime_utr"),
               names_to = "region",
               values_to = "gc_content") %>% 
  mutate(region = factor(region, levels = c("five_prime_utr", "first_100", "three_prime_utr")))

region_labels <- c(
  "five_prime_utr" = "5' UTR",
  "first_100" = "Primeros 100",
  "three_prime_utr" = "3' UTR"
)

n <- lncRNA_HCC_2024 %>%
  group_by(size, region) %>%
  summarise(n = n()) %>%
  ungroup()

plot3 <- ggplot(lncRNA_HCC_2024, aes(x = size, y = gc_content)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  labs(x = "Tamaño relativo", y = "Porcentaje de GC") +
  geom_violin() +
  scale_x_discrete(labels = c("Grande", "Mediano", "Pequeño")) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_wrap(~ region, ncol = 3, labeller = labeller(region = region_labels)) +
  ylim(0, 100) +
  geom_text(data = n, aes(x = interaction(size), y = 95, label = paste(n)),
            position = position_dodge(width = 0.75), family = "CMU Serif", size = 4.5)

# Combine all plots
final_plot <- plot_grid(
  plot1, plot2, plot3,
  ncol = 1,
  labels = c("(a)", "(b)", "(c)"),
  label_fontfamily = "CMU Serif",
  label_size = 12
)

final_plot

# Save the combined plot ----

padded_plot1 <- plot_grid(
  NULL, plot1, NULL,
  nrow = 1,
  rel_widths = c(0.2, 0.6, 0.2)
)

print(padded_plot1)


padded_plot2 <- plot_grid(
  NULL, plot2, NULL,
  nrow = 1,
  rel_widths = c(0.3,0.4,0.3)
)

print(padded_plot2)

plot_grid(
  padded_plot1, padded_plot2, plot3,
  ncol = 1,
  labels = c("(c)", "(d)", "(e)"),
  label_fontfamily = "CMU_Serif",
  label_size = 10
)

ggsave(
  plot = last_plot(),
  filename = "./DOCUMENTS/Memoria/Images/gc_content/combined_gc_content_part_2.pdf",
  device = "pdf",
  width = 16,
  height = 20,
  units = "cm"
)
