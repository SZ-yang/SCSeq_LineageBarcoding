rm(list=ls())

library(Seurat)
library(GEMLI)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))

# first filter to the top 50 lineages
lineage_tab <- table(seurat_obj$clone_id)
largest_lineages <- names(lineage_tab)[order(lineage_tab, decreasing = TRUE)[1:50]]
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[seurat_obj$clone_id %in% largest_lineages] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# remove any cells that are empty
data_matrix <- SeuratObject::LayerData(seurat_obj,
                                       layer = "data",
                                       assay = "RNA")
sum_vec <- Matrix::colSums(data_matrix)
if(any(sum_vec == 0)){
  keep_vec <- rep(TRUE, length(Seurat::Cells(seurat_obj)))
  keep_vec[sum_vec == 0] <- FALSE
  seurat_obj$keep <- keep_vec
  seurat_obj <- subset(seurat_obj, keep == TRUE)
}

# run GEMLI
data_matrix <- SeuratObject::LayerData(seurat_obj,
                                       layer = "data",
                                       assay = "RNA")
data_matrix <- as.matrix(data_matrix)
rowsum_vec <- Matrix::rowSums(data_matrix)
data_matrix <- data_matrix[rowsum_vec != 0,]

GEMLI_items <- list()
GEMLI_items[['gene_expression']] <- data_matrix
GEMLI_items[['barcodes']] <- seurat_obj$clone_id

set.seed(10)
GEMLI_items <- predict_lineages_custom(GEMLI_items,
                                       desired_cluster_size = c(50, 200),
                                       verbose = 1)

# Error in hclust(as.dist(corr_expr_subset), method = "ward.D2", ) : 
# NA/NaN/Inf in foreign function call (arg 10)
# [[This is from quantify_clusters_iterative... https://github.com/UPSUTER/GEMLI/blob/main/GEMLI_package_v0/R/quantify_clusters_iterative.R#L19]]




