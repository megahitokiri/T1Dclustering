# Assuming mf is your merged table
library(dplyr)
library(ggpubr)

# Focus on dataset C, full clinical
df_C_full <- mf %>% filter(dataset == "C", feature == "genetic_fullclin")

# 1️⃣ Check PC1–PC5 distributions by cluster
pc_vars <- paste0("PC", 1:5)
p_tests <- sapply(pc_vars, function(pc)
  wilcox.test(df_C_full[[pc]] ~ df_C_full$cluster)$p.value)
p_tests

# 2️⃣ Estimate variance explained by ancestry PCs
summary(lm(as.numeric(cluster) ~ PC1 + PC2 + PC3 + PC4 + PC5, data = df_C_full))