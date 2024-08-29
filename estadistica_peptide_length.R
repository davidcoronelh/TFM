library(emmeans)
library(dplyr)

load("DATA/datasets/all_datasets.RData")

BRE <- bind_rows(CDSs_BRE, ORFs_BRE) %>% mutate(group = "BRE")
HCC <- bind_rows(CDSs_HCC, ORFs_HCC) %>% mutate(group = "HCC")
HCC_2024 <- bind_rows(CDSs_HCC_2024, ORFs_HCC_2024) %>% mutate(group = "HCC_2024")
all_ORFs <- bind_rows(BRE, HCC, HCC_2024) %>% mutate(log10_len = log10(len))

# ajustar el modelo lineal
model <- lm(log10_len ~ group * class, data = all_ORFs)

# extraer los residuos
residuals <- residuals(model, type = "deviance")


par(mfrow = c(2, 2))

# qq plot de los residuos
qqnorm(residuals)
qqline(residuals, col = "red")

# histograma
hist(residuals, breaks = 30)

# residuos vs valores ajustados
plot(model$fitted.values, residuals)
abline(h = 0, col = "red")

# grafico de escala-ubicacion
plot(sqrt(abs(residuals)) ~ model$fitted.values, 
abline(h = 0, col = "red")

# calcular las medias marginales estimadas y graficarlas
emms <- emmeans(model, ~ class * group)
plot(emms, comparisons = TRUE)

# realizar la comparacion pairwise entre todas las combinaciones de clase y grupo con ajuste del p-valor
pairwise_comparisons <- contrast(emms, method = "pairwise", adjust = "bonferroni")

