rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

train_embedding <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup9/train_adata_LCL.csv")
rownames(train_embedding) <- train_embedding[,1]
train_embedding <- as.matrix(train_embedding[,-1])

test_embedding <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup9/test_adata_LCL.csv")
rownames(test_embedding) <- test_embedding[,1]
test_embedding <- as.matrix(test_embedding[,-1])

table(rownames(train_embedding) %in% Seurat::Cells(seurat_obj))
table(rownames(test_embedding) %in% Seurat::Cells(seurat_obj))

LCL_embedding <- matrix(NA, nrow = length(Seurat::Cells(seurat_obj)), ncol = ncol(train_embedding))
rownames(LCL_embedding) <- Seurat::Cells(seurat_obj)
colnames(LCL_embedding) <- paste0("LCL_", 1:ncol(LCL_embedding))
LCL_embedding[rownames(train_embedding),] <- train_embedding
LCL_embedding[rownames(test_embedding),] <- test_embedding

seurat_obj[["lcl"]] <- Seurat::CreateDimReducObject(LCL_embedding)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "lcl", 
                              dims = 1:ncol(LCL_embedding),
                              assay = 'RNA', 
                              reduction.name = 'lcl.umap', 
                              reduction.key = 'lclUMAP_')
testing_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
names(testing_vec) <- Seurat::Cells(seurat_obj)
testing_vec[rownames(test_embedding)] <- TRUE
seurat_obj$test <- testing_vec

scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "OG_condition", 
                              colors_use = seurat_obj@misc[["OG_condition_colors"]])

scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "test")


color_vec <- c(Lineage1 = rgb(203, 50, 39, maxColorValue = 255),
               Lineage2 = rgb(65, 113, 170, maxColorValue = 255),
               Lineage3 = rgb(131, 72, 149, maxColorValue = 255),
               Lineage4 = rgb(92, 163, 77, maxColorValue = 255),
               Lineage5 = rgb(235, 123, 44, maxColorValue = 255),
               other = "gray")

tab_vec <- table(seurat_obj$Lineage)
tab_vec <- sort(tab_vec, decreasing = TRUE)
largest_lineage_vec <- rep("other", length(Seurat::Cells(seurat_obj)))
names(largest_lineage_vec) <- Seurat::Cells(seurat_obj)
for(i in 1:5){
  lineage_name <- names(tab_vec)[i]
  largest_lineage_vec[which(seurat_obj$Lineage == lineage_name)] <- paste0("Lineage", i)
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
                              reduction = "umap", 
                              group.by = "largest_lineage", 
                              colors_use = color_vec,
                              order = paste0("Lineage", 1:5))




