library(showtext)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)


showtext_auto()


font_add(family = "CMU_Serif", 
         regular = "~/Library/Fonts/computer-modern/cmunrm.ttf", 
         bold = "~/Library/Fonts/computer-modern/cmunbx.ttf", 
         italic = "~/Library/Fonts/computer-modern/cmunti.ttf", 
         bolditalic = "~/Library/Fonts/computer-modern/cmunbi.ttf")


load("DATA/datasets/ORFs_HCC_2024.RData")
load("DATA/datasets/CDSs_HCC_2024.RData")
gcContentCalc <- readRDS("FUNCTIONS/gcContentCalc.rds")

five_prime_analysis <- list(CDSs_HCC_2024, ORFs_HCC_2024) %>%
  lapply(function(x) {
    positive <- x[x$sense == "+", ]
    negative <- x[x$sense == "-", ]
    
    positive_filtrados <- positive %>% 
      group_by(chromosome, posicion_inicial) %>%
      arrange(-len) %>% 
      slice_head()
    
    negative_filtrados <- negative %>% 
      group_by(chromosome, posicion_inicial) %>%
      arrange(-len) %>% 
      slice_head()
    combined <- rbind(positive_filtrados, negative_filtrados)
    return(combined)
  }) %>%
  do.call(rbind, .) %>%
  rowwise() %>% 
  mutate(five_prime_utr = gcContentCalc(substr(exp_seq, 1, 100)),
         first_100 = gcContentCalc(substr(cDNA, 4, 103))) %>% 
  ungroup()

three_prime_analysis <- list(CDSs_HCC_2024, ORFs_HCC_2024) %>% 
  lapply(function(x) {
    
    positive <- x[x$sense == "+", ]
    negative <- x[x$sense == "-", ]
    
    positive_filtrados <- positive %>% 
      group_by(chromosome, posicion_inicial) %>%
      arrange(-len) %>% 
      slice_head()
    
    negative_filtrados <- negative %>% 
      group_by(chromosome, posicion_inicial) %>%
      arrange(-len) %>% 
      slice_head()
    

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

combined_data <- bind_rows(five_prime_analysis, three_prime_analysis, .id = "tabla") %>% 
  pivot_longer(cols = c(five_prime_utr, first_100, three_prime_utr),
               names_to = "region",
               values_to = "gc_content") %>% 
  mutate(region = factor(region, levels = c("five_prime_utr", "first_100", "three_prime_utr")),
         class = factor(class, levels = c("CDS", "uORF", "uoORF", "lncRNA-ORF")),
         start_codon = ifelse(substr(cDNA, 1, 3) == "ATG", "ATG", "non_ATG")) %>% 
  filter(!(class == "CDS" & start_codon != "ATG")) %>% 
  select(orf, class, gene, sense, posicion_inicial, posicion_final,region, start_codon, gc_content, cDNA) %>% 
  drop_na(gc_content)

region_labels <- c(
  "five_prime_utr" = "5' UTR",
  "first_100" = "Primeros 100",
  "three_prime_utr" = "3' UTR"
)

class_labels <- c(
  "CDS" = "CDS",
  "uORF" = "uORF",
  "uoORF" = "uoORF",
  "lncRNA-ORF" = "lncRNA"
)

add_counts <- function(x){
  return(data.frame(y = 99, label = length(x)))
}


small_plot <- ggplot(combined_data, aes(x = class, y = gc_content)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  labs(x = "Clase",
       y = "Porcentaje de GC") +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_grid(~region, labeller = labeller(region = region_labels)) +
  ylim(0, 100) +
  scale_x_discrete(labels = c("CDS", "uORF", "uoORF", "lncRNA"),
                   guide = guide_axis(angle = 25)) +
  stat_summary(fun.data = add_counts, geom = "text", aes(label = ..label.., family = "CMU Serif"), size = 3, vjust = 0)

# ggsave(plot = last_plot(),
#        filename = "DOCUMENTS/Memoria/Images/gc_content_small.pdf",
#        width = 15,
#        height = 7,
#        units = "cm",
#        device = "pdf")



######### Large plot ----
max_gc <- combined_data %>%
  group_by(start_codon, region, class) %>%
  summarize(max_gc_content = max(gc_content, na.rm = TRUE))

counts_large <- combined_data %>%
  group_by(region, class, start_codon) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  left_join(max_gc, by = c("start_codon", "region", "class")) %>%
  mutate(y_position = 100) # Adjust the offset as needed

large_plot <- ggplot(combined_data, aes(x = start_codon, y = gc_content)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  labs(x = "Codón de inicio",
       y = "Porcentaje de GC") +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_grid(region ~ class,
             labeller = labeller(region = region_labels,
                                 class = class_labels)) +
  ylim(0, 100) +
  scale_x_discrete(labels = c("ATG", "NTG")) +
  geom_text(data = counts_large, aes(x = start_codon, y = y_position, label = paste(n)),
            position = position_dodge(width = 0.75), family = "CMU Serif", size = 3)

# ggsave(plot = last_plot(),
#        filename = "DOCUMENTS/Memoria/Images/gc_content_large.pdf",
#        width = 15,
#        height = 12,
#        units = "cm",
#        device = "pdf")

plot_grid(
  small_plot, large_plot,
  ncol = 1,
  labels = c("(a)", "(b)"),
  rel_heights = c(0.4, 0.6),
  label_fontfamily = "CMU_Serif",
  label_size = 10
)

ggsave(
  plot = last_plot(),
  filename = "./DOCUMENTS/Memoria/Images/gc_content/combined_gc_content_part_1.pdf",
  device = "pdf",
  width = 16,
  height = 18,
  units = "cm"
)


# ***************************************************
# ******* STATISTICS OF THE LARGE PLOT ************** ----
# ***************************************************

permutation_test <- function(data, observations_col, groups_col, groupA, groupB, P = 100000, seed = 1234) {
  # Extract the relevant columns 
  observations <- data[[observations_col]]
  groups <- data[[groups_col]]
  
  # Filter the data for the two specified groups
  data <- data[groups %in% c(groupA, groupB), ]
  observations <- data[[observations_col]]
  groups <- data[[groups_col]]
  
  # Calculate means and medians for the groups
  group_means <- with(data, tapply(observations, groups, mean))
  # group_medians <- with(data, tapply(observations, groups, median))
  
  # Calculate test statistics
  test.stat1 <- abs(diff(group_means))
  # test.stat2 <- abs(diff(group_medians))
  
  
  # PERMUTATION TEST
  set.seed(seed)
  n <- length(observations)
  
  PermSamples <- matrix(0, nrow = n, ncol = P)
  
  # Create the permutation samples matrix
  PermSamples <- replicate(P, sample(observations, size = n, replace = FALSE))
  
  
  # Using apply to replace the for loop
  Perm.test.stat1 <- apply(PermSamples, 2, function(column) {
    abs(mean(column[data[[groups_col]] == groupA]) - median(column[data[[groups_col]] == groupB]))
  })
  
  # Perm.test.stat2 <- apply(PermSamples, 2, function(column) {
  #   abs(median(column[data[[groups_col]] == groupA]) - median(column[data[[groups_col]] == groupB]))
  # })
  
  # Calculate the p-values
  p_value1 <- mean(Perm.test.stat1 >= test.stat1)
  # p_value2 <- mean(Perm.test.stat2 >= test.stat2)
  
  # Return the results as a dataframe
  return(data.frame(
    class = unique(data$class),
    region = unique(data$region),
    pvalue_mean = p_value1
  ))
  
}

perm_results_df <- data.frame()

# Define the classes and regions
classes <- c("uORF", "uoORF", "lncRNA-ORF")
regions <- c("five_prime_utr", "first_100", "three_prime_utr")

# Iterate over each combination of class and region
for (cl in classes) {
  for (re in regions) {
    
    # Filter and select the data
    library(dplyr)
    data <- combined_data %>%
      filter(class == cl & region == re) %>%
      dplyr::select(orf, class, region, gc_content, start_codon)
    
    # Perform analysis and store the result
    perm_result <- permutation_test(
      data = data,
      observations_col = "gc_content",
      groups_col = "start_codon",
      groupA = "ATG",
      groupB = "non_ATG"
    )
    
    # Store the result in the list
    perm_results_df <- rbind(perm_results_df, perm_result)
  }
}

p_value_labels <- perm_results_df %>% 
  rowwise() %>% 
  mutate(pvalue = ifelse(pvalue_mean < 0.0001, "< 0.0001", paste("=", as.character( round( pvalue_mean, 4) ) ) ))

# Determine the y-position for the p-value labels
pval_y_position <- 12  

# Create the plot with the p-value text
large_plot <- ggplot(combined_data, aes(x = start_codon, y = gc_content)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif")) +
  labs(x = "Codón de inicio",
       y = "Porcentaje de GC") +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_grid(region ~ class,
             labeller = labeller(region = region_labels,
                                 class = class_labels)) +
  ylim(0, 100) +
  scale_x_discrete(labels = c("ATG", "NTG")) +
  geom_text(data = counts_large, aes(x = start_codon, y = y_position, label = paste(n)),
              position = position_dodge(width = 0.75), family = "CMU Serif", size = 3) +
  geom_text(data = p_value_labels,
            aes(x = 1.5, y = pval_y_position, label = paste("p", pvalue)),
            family = "CMU Serif", size = 3.5)

# Print the plot
print(large_plot)
