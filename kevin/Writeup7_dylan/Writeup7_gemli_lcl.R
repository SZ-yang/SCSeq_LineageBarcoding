rm(list=ls())

library(Seurat)

source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup5_GEMLI/quantify_clusters_iterative_custom.R")
source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup5_GEMLI/predict_lineages_custom.R")

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/"
plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

seurat_obj <- NormalizeData(seurat_obj)

# first filter to the top 50 lineages
lineage_tab <- table(seurat_obj$Lineage)
largest_lineages <- names(lineage_tab)[order(lineage_tab, decreasing = TRUE)[1:50]]
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[seurat_obj$Lineage %in% largest_lineages] <- TRUE
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
data_matrix <- t(seurat_obj[["lcl"]]@cell.embeddings)

GEMLI_items <- list()
GEMLI_items[['gene_expression']] <- data_matrix
GEMLI_items[['barcodes']] <- seurat_obj$Lineage

set.seed(10)
GEMLI_items <- predict_lineages_custom(GEMLI_items,
                                       bool_find_markers = FALSE,
                                       desired_cluster_size = c(50,200),
                                       verbose = 1)

save(date_of_run, session_info,
     GEMLI_items,
     file = paste0(out_folder, "Writeup7_GEMLI_lcl.RData"))

zz <- GEMLI_items[['prediction']] # we need to look at how they did it on LARRY...
quantile(zz[zz!=0])

GEMLI_items <- test_lineages_custom(GEMLI_items, 
                                    valid_fam_sizes = c(50,200))
GEMLI_items$testing_results

save(date_of_run, session_info,
     GEMLI_items,
     file = paste0(out_folder, "Writeup7_GEMLI_lcl.RData"))

print("Done! :)")

