rm(list=ls())
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

seurat_obj$Lineage <- factor(paste0("Lineage:", seurat_obj$Lineage))

tab_mat <- table(seurat_obj$Lineage, seurat_obj$OG_condition)
rowsum_vec <- rowSums(tab_mat)
tab_mat2 <- tab_mat

# normalize
for(i in 1:nrow(tab_mat2)){
  tab_mat2[i,] <- tab_mat2[i,]/sum(tab_mat2[i,])
}

n <- nrow(tab_mat2)
cor_mat <- matrix(NA, n, n)
rownames(cor_mat) <- rownames(tab_mat)
colnames(cor_mat) <- rownames(tab_mat)
for(i in 1:(n-1)){
  if(i %% floor(n/10) == 0) cat('*')
  idx1 <- which(tab_mat2[i,] != 0)
  
  for(j in (i+1):n){
    idx2 <- which(tab_mat2[j,] != 0)
    if(length(idx1) == length(idx2) & all(idx1 %in% idx2) & length(idx1) > 1){
      tmp <- stats::cor(tab_mat2[i,], tab_mat2[j,])
      cor_mat[i,j] <- tmp
      cor_mat[j,i] <- tmp
    }
  }
}

table(!is.na(cor_mat))
quantile(cor_mat[!is.na(cor_mat)])

idx <- which(cor_mat > 0.9, arr.ind = TRUE)
unique_lineage_idx <- unique(as.numeric(idx))
quantile(rowsum_vec[unique_lineage_idx])
largest_lineage_idx <- unique_lineage_idx[which.max(rowsum_vec[unique_lineage_idx])]
rowsum_vec[largest_lineage_idx]

# zoom in on this lineage
parter_lineage_idx <- which(cor_mat[largest_lineage_idx,] >= 0.9)
rowsum_vec[parter_lineage_idx]

lineage_idx <- c(largest_lineage_idx, parter_lineage_idx)
cor_mat[lineage_idx, lineage_idx]
tab_mat[lineage_idx,]
tab_mat2[lineage_idx,]

lineage_names <- rownames(cor_mat)[lineage_idx]

##############

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    reduction = "lcl",
                                    dims = 1:64)

seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.1)

seurat_obj$RNA_snn_res.0.1 <- factor(paste0("Cluster:", seurat_obj$RNA_snn_res.0.1))

# find all the cells that are in the lineage_names
cell_idx <- which(seurat_obj$Lineage %in% lineage_names)
tab_mat_special <- table(
  droplevels(seurat_obj$Lineage[cell_idx]),
  droplevels(seurat_obj$RNA_snn_res.0.1[cell_idx])
)

# hmmm.....

# how many lienages are in each island?
cluster_names <- levels(seurat_obj$RNA_snn_res.0.1)
num_lineages_per_cluster <- sapply(cluster_names, function(cluster){
  tmp_idx <- which(seurat_obj$RNA_snn_res.0.1 == cluster)
  length(unique(seurat_obj$Lineage[tmp_idx]))
})
names(num_lineages_per_cluster) <- cluster_names
num_lineages_per_cluster[unique(seurat_obj$RNA_snn_res.0.1[cell_idx])]

hist(num_lineages_per_cluster)
