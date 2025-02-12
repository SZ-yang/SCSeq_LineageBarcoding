rm(list=ls())
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

tab_mat <- table(seurat_obj$Lineage, seurat_obj$OG_condition)
tab_mat <- log10(tab_mat+1)

pdf(paste0(plot_folder, "Writeup7_lineage-size-scatterplot.pdf"),
    height = 5, width = 5)
p <- ncol(tab_mat)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    x <- tab_mat[,i]
    x2 <- x + runif(length(x), min = 0, max = 0.2)
    y <- tab_mat[,j]
    y2 <- y + runif(length(y), min = 0, max = 0.2)
    cor_val <- stats::cor(x,y)
    
    plot(x = x2, y = y2,
         col = rgb(0.5,0.5,0.5,0.5),
         xlab = colnames(tab_mat)[i],
         ylab = colnames(tab_mat)[j],
         asp = TRUE,
         pch = 16,
         main = paste0("Correlation: ", round(cor_val,2)))
  }
}
dev.off()

#########

# Load required packages
library(ggplot2)
library(ggcorrplot)
library(corrplot)

# Function to compute correlation matrix and p-values
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

# Compute correlation matrix and p-values
cor_matrix <- cor(tab_mat, use = "pairwise.complete.obs", method = "pearson")
p_matrix <- cor.mtest(tab_mat)

# Create heatmap with formatted correlation values and significance stars
plot1 <- ggcorrplot(cor_matrix, 
                    method = "square",  # Use "square" since "color" is not a valid method
                    type = "full", 
                    lab = TRUE, 
                    lab_col = "black",
                    lab_size = 3,
                    outline.color = "white",
                    colors = c("blue", "white", "red"),
                    tl.cex = 12,
                    digits = 2,
                    p.mat = p_matrix, # Provide p-values for significance testing
                    sig.level = c(0.05, 0.01, 0.001), # Significance levels
                    insig = "blank") + # Hide insignificant values
  theme_minimal() +
  ggtitle("Correlation Heatmap") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup7_lineage-size_correlation-heatmap.png"),
                height = 6, width = 6)
