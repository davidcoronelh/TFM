library(lme4)
library(dplyr)
library(emmeans)
library(influence.ME)

data <- filtered_data_for_plot %>% 
  mutate(log10_ti = log10(translation_index))

# ajuste del modelo lineal mixto
mixed_model <- lmer(log10_ti ~ patient_id + (1 | orf), data = data)


# calcular la distancia de Cook 
infl <- influence(mixed_model, group = "orf")
cooksd <- cooks.distance(infl)

# graficar la distancia de Cook
plot(cooksd, type = "b",
     pch = 18,
     col = "firebrick2",
     xlab = "Índice del ORF",
     ylab = "Distancia de Cook")

# definir un valor de corte para identificar puntos influyentes
n <- length(unique(data$orf))
cutoff <- 4 / n
abline(h = cutoff, lty = 2)


# filtrar datos basados en la distancia de Cook
below_threshold <- cooksd[cooksd <= cutoff,  ,drop = FALSE]

# filtrar el conjunto de datos para eliminar los puntos influyentes
data_cooked <- data %>% 
  filter(orf %in% rownames(below_threshold))

# histograma del índice de traducción después de eliminar los puntos
hist(data_cooked$translation_index)

# número de filas eliminadas
nrow(data) - nrow(data_cooked)

# volver a ajustar el modelo
cooked_mixed_model <- lmer(log10_ti ~ patient_id + (1 | orf), data = data_cooked)

# hacer el contraste por pares
print(emmeans(cooked_mixed_model, pairwise ~ patient_id, adjust = "tukey"))

