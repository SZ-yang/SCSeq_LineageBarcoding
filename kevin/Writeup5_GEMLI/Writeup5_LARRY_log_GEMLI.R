rm(list=ls())

library(Seurat)
library(GEMLI)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()

source("predict_lineages_custom.R")
source("quantify_clusters_iterative_custom.R")

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
                                       desired_cluster_size = c(50,200),
                                       verbose = 1)

save(date_of_run, session_info,
     GEMLI_items,
     file = paste0(out_folder, "Larry_log_GEMLI.RData"))

zz <- GEMLI_items[['prediction']] # we need to look at how they did it on LARRY...
quantile(zz[zz!=0])

GEMLI_items <- GEMLI::test_lineages(GEMLI_items, 
                                    valid_fam_sizes = 50:200)
GEMLI_items$testing_results

save(date_of_run, session_info,
     GEMLI_items,
     file = paste0(out_folder, "Writeup5_Larry_log_GEMLI.RData"))

print("Done! :)")
