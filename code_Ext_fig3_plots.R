### General information ----
# Title: Establishing a random forest model to predict the maturity of adult, human plasma cells using various cell surface marker
# Author: Tobit D. Steinmetz
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: d.t.steinmetz@umcg.nl
# Collaboration: please ask permission from the author before using this script
# Date created: 19-10-2023
# Date last adjustment: 03-12-2024
# RStudio version: 2023.06.1 Build 524
# References: 

### introduction ----
# this scripts contains the code to generate heatmaps for Extended Figure 3
# correlation matrices are generated for DPI and maturity markers with overlays for pair numbers and p-values
# 

### This R-script requires the following packages:
library(ggplot2)
library(reshape2)
library(dplyr)

## contine with R environment from training the RF model
load(".../ASC_ME/ASC_ME_environment.RData")
setwd(".../ASC_ME")

# Function to calculate the correlation p-value for each pair of columns
cor_pvalue <- function(x, y) {
  cor_test <- cor.test(x, y, method = "pearson")
  return(cor_test$p.value)
}

# Subset the data you want to analyze (columns 2 to 10 in this case)
data_subset <- data[, c(2:10)] %>% filter(dpi>4)
correlation_matrix <- cor(data_subset, method = "pearson", use = "pairwise.complete.obs")
print(correlation_matrix)

# Initialize a matrix to store the p-values
p_value_matrix <- matrix(NA, ncol = ncol(data_subset), nrow = ncol(data_subset))
rownames(p_value_matrix) <- colnames(data_subset)
colnames(p_value_matrix) <- colnames(data_subset)

# Calculate p-values for each pair of variables
for (i in 1:ncol(data_subset)) {
  for (j in 1:ncol(data_subset)) {
    if (i != j) {  # Avoid calculating p-value for a variable with itself
      x <- data_subset[[i]]
      y <- data_subset[[j]]
      p_value_matrix[i, j] <- cor_pvalue(x, y)
    }
  }
}

# Create a matrix to store the number of pairs for each correlation
n_pairs_matrix <- outer(
  1:ncol(data_subset), 1:ncol(data_subset), 
  Vectorize(function(i, j) sum(complete.cases(data_subset[, c(i, j)])))
)
colnames(n_pairs_matrix) <- c("dpi","CD19","CD20","CD28","CD45","CD56","CD138","HLA.DR","Ki67")
rownames(n_pairs_matrix) <- c("dpi","CD19","CD20","CD28","CD45","CD56","CD138","HLA.DR","Ki67")

# Melt the correlation matrix and p-value matrix into long format
cor_melt <- melt(correlation_matrix)
pval_melt <- melt(p_value_matrix)
pair_melt <- melt(n_pairs_matrix)

# Rename the columns for clarity
colnames(cor_melt) <- c("X1", "X2", "Correlation")
colnames(pval_melt) <- c("X1", "X2", "PValue")
colnames(pair_melt) <- c("X1", "X2", "npair")

# Create the heatmap
ggplot(cor_melt, aes(X1, X2)) +
  geom_tile(aes(fill = Correlation), color = "white") +  # Heatmap for correlation matrix
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Pearson\nCorrelation") +  # Color scale
  geom_text(data = pval_melt, aes(X1, X2, label = formatC(PValue, format = "e", digits = 2)), color = "black", size = 4, nudge_y = -0.2) +  # Overlay p-values in scientific notation
  geom_text(data = pair_melt, aes(X1, X2, label = npair), color = "black", size = 4, nudge_y = 0.2, fontface = "bold") +  
  theme_bw() +
  ylim("Ki67","HLA.DR","CD138","CD56","CD45","CD28","CD20","CD19","dpi") +
  theme(axis.text.y = element_text(size = 15,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 15,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  labs(x = "", y = "", title = "correlation heatmap for unimputed data")
ggsave("plot_heatmap_cor_unimp_data.png")

###
###heatmap for imputed data
#
# Subset the data you want to analyze (columns 2 to 10 in this case)
data_subset_imp <- data.imputed_mean[, c(2,4:11)] %>% filter(dpi>4)
correlation_matrix <- cor(data_subset_imp, method = "pearson", use = "pairwise.complete.obs")
# Initialize a matrix to store the p-values
p_value_matrix <- matrix(NA, ncol = ncol(data_subset_imp), nrow = ncol(data_subset_imp))
rownames(p_value_matrix) <- colnames(data_subset_imp)
colnames(p_value_matrix) <- colnames(data_subset_imp)

# Calculate p-values for each pair of variables
for (i in 1:ncol(data_subset_imp)) {
  for (j in 1:ncol(data_subset_imp)) {
    if (i != j) {  # Avoid calculating p-value for a variable with itself
      x <- data_subset_imp[[i]]
      y <- data_subset_imp[[j]]
      p_value_matrix[i, j] <- cor_pvalue(x, y)
    }
  }
}

# Create a matrix to store the number of pairs for each correlation
n_pairs_matrix <- outer(
  1:ncol(data_subset_imp), 1:ncol(data_subset_imp), 
  Vectorize(function(i, j) sum(complete.cases(data_subset_imp[, c(i, j)])))
)
colnames(n_pairs_matrix) <- c("dpi","CD19","CD20","CD28","CD45","CD56","CD138","HLA.DR","Ki67")
rownames(n_pairs_matrix) <- c("dpi","CD19","CD20","CD28","CD45","CD56","CD138","HLA.DR","Ki67")

# Melt the correlation matrix and p-value matrix into long format
cor_melt <- melt(correlation_matrix)
pval_melt <- melt(p_value_matrix)
pair_melt <- melt(n_pairs_matrix)

# Rename the columns for clarity
colnames(cor_melt) <- c("X1", "X2", "Correlation")
colnames(pval_melt) <- c("X1", "X2", "PValue")
colnames(pair_melt) <- c("X1", "X2", "npair")

# Create the heatmap
ggplot(cor_melt, aes(X1, X2)) +
  geom_tile(aes(fill = Correlation), color = "white") +  # Heatmap for correlation matrix
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Pearson\nCorrelation") +  # Color scale
  geom_text(data = pval_melt, aes(X1, X2, label = formatC(PValue, format = "e", digits = 2)), color = "black", size = 4, nudge_y = -0.2) +  # Overlay p-values in scientific notation
  #geom_text(data = pair_melt, aes(X1, X2, label = npair), color = "black", size = 4, nudge_y = 0.2, fontface = "bold") +  
  theme_bw() +
  ylim("Ki67","HLA.DR","CD138","CD56","CD45","CD28","CD20","CD19","dpi") +
  theme(axis.text.y = element_text(size = 15,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 15,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  labs(x = "", y = "", title = "correlation heatmap for imputed data")
ggsave("plot_heatmap_cor_imp_data.png")
