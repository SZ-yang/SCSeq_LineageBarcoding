rm(list=ls())

library(Seurat)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup9b/"
load(paste0(out_folder, "Writeup9b_celltag-data_with_LCL.RData"))

train_metadata <- read.csv(paste0(out_folder, "train_obs.csv"))
rownames(train_metadata) <- train_metadata[,1]

test_metadata <- read.csv(paste0(out_folder, "test_obs.csv"))
rownames(test_metadata) <- test_metadata[,1]

clone_vec <- rep(NA, length(Seurat::Cells(seurat_obj)))
names(clone_vec) <- Seurat::Cells(seurat_obj)
clone_vec[rownames(train_metadata)] <- train_metadata$clone_id
clone_vec[rownames(test_metadata)] <- test_metadata$clone_id

seurat_obj$clone_id <- clone_vec

testing_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
names(testing_vec) <- Seurat::Cells(seurat_obj)
testing_vec[rownames(test_metadata)] <- TRUE
seurat_obj$test <- testing_vec

set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "lcl", 
                              dims = 1:ncol(seurat_obj[["lcl"]]@cell.embeddings),
                              assay = 'RNA', 
                              reduction.name = 'lcl.umap', 
                              reduction.key = 'lclUMAP_')

scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "cell_type")

scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "test")

#############


color_vec <- c(Lineage1 = rgb(203, 50, 39, maxColorValue = 255),
               Lineage2 = rgb(65, 113, 170, maxColorValue = 255),
               Lineage3 = rgb(131, 72, 149, maxColorValue = 255),
               Lineage4 = rgb(92, 163, 77, maxColorValue = 255),
               Lineage5 = rgb(235, 123, 44, maxColorValue = 255),
               other = "gray")

tab_vec <- table(seurat_obj$clone_id)
tab_vec <- sort(tab_vec, decreasing = TRUE)
largest_lineage_vec <- rep("other", length(Seurat::Cells(seurat_obj)))
names(largest_lineage_vec) <- Seurat::Cells(seurat_obj)
for(i in 1:5){
  lineage_name <- names(tab_vec)[i]
  largest_lineage_vec[which(seurat_obj$clone_id == lineage_name)] <- paste0("Lineage", i)
}
seurat_obj$largest_lineage <- largest_lineage_vec

Seurat::Idents(seurat_obj) <- "largest_lineage"
scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "largest_lineage", 
                              colors_use = color_vec,
                              order = paste0("Lineage", 1:5))

Seurat::Idents(seurat_obj) <- "largest_lineage"
scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "python_X_tsne", 
                              group.by = "largest_lineage", 
                              colors_use = color_vec,
                              order = paste0("Lineage", 1:5))


