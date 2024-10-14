rm(list=ls())

library(Seurat)
library(GEMLI)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()

source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup5_GEMLI/memory_gene_calling_custom.R")

#########

data_folder <- "~/kzlinlab/data/larry_hematopoiesis/"
metadata_df <- read.csv(paste0(data_folder, "stateFate_inVitro_metadata.txt"), sep = "\t")
traj_df <- read.csv(paste0(data_folder, "stateFate_inVitro_neutrophil_monocyte_trajectory.txt"))

tmp <- sapply(1:nrow(metadata_df), function(i){
  paste0(metadata_df$Library[i], ":", metadata_df$Cell.barcode[i])
})
metadata_df$linBarcode <- tmp
traj_vec <- rep(FALSE, nrow(metadata_df))
traj_vec[traj_df[,1]+1] <- TRUE
metadata_df$trajectory <- traj_vec

metadata_df <- metadata_df[metadata_df$trajectory == TRUE,]
metadata_df <- metadata_df[metadata_df$Time.point == 2,]
table(metadata_df$Cell.type.annotation)

##########

tmp_df <- seurat_obj@meta.data[,c("Library", "Cell.barcode")]
tmp_df$Library <- as.character(tmp_df$Library)
tmp_df$Cell.barcode <- as.character(tmp_df$Cell.barcode)
tmp <- sapply(1:nrow(tmp_df), function(i){
  paste0(tmp_df[i,], collapse = ":")
})
seurat_obj$linBarcode <- tmp
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[seurat_obj$linBarcode %in% metadata_df$linBarcode] <- TRUE
seurat_obj$keep <- keep_vec

seurat_obj <- subset(seurat_obj, keep == TRUE)

table(table(seurat_obj$clone_id))

##############

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

##################

GEMLI_items <- memory_gene_calling_custom(GEMLI_items = GEMLI_items, 
                                          ground_truth = TRUE,
                                          num_trials = 100,
                                          use_median = TRUE, 
                                          valid_lineage_sizes = c(1:10), 
                                          verbose = 2)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
save(GEMLI_items,
     date_of_run, session_info, 
     file = paste0(out_folder, "Larry_GEMLI-memory-genes_day2.RData"))



