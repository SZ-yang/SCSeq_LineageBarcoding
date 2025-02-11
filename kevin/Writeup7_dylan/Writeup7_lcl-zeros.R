rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dbscan)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"
source("kevin/Writeup7_dylan/welch_anova.R")

load(paste0(out_folder, "adata_with_lcl.RData"))

embedding <- seurat_obj[["lcl"]]@cell.embeddings

num_nonzeros_cells <- apply(embedding, 1, function(x){length(which(x > 0))})
hist(num_nonzeros_cells)

num_nonzeros_features <- apply(embedding, 2, function(x){length(which(x > 0))})
hist(num_nonzeros_features)
