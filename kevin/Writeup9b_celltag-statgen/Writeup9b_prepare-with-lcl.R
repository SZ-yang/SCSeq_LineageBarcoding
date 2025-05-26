rm(list=ls())

library(Seurat)

load("~/kzlinlab/data/biddy_2018_celltag/biddy_seurat.RData")

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup9b/"

train_embedding <- read.csv(paste0(out_folder, "train_adata_LCL.csv"))
rownames(train_embedding) <- train_embedding[,1]
train_embedding <- as.matrix(train_embedding[,-1])

test_embedding <- read.csv(paste0(out_folder, "test_adata_LCL.csv"))
rownames(test_embedding) <- test_embedding[,1]
test_embedding <- as.matrix(test_embedding[,-1])

table(rownames(train_embedding) %in% Seurat::Cells(seurat_obj))
table(rownames(test_embedding) %in% Seurat::Cells(seurat_obj))

LCL_embedding <- matrix(NA, nrow = length(Seurat::Cells(seurat_obj)), ncol = ncol(train_embedding))
rownames(LCL_embedding) <- Seurat::Cells(seurat_obj)
colnames(LCL_embedding) <- paste0("LCL_", 1:ncol(LCL_embedding))
LCL_embedding[rownames(train_embedding),] <- train_embedding
LCL_embedding[rownames(test_embedding),] <- test_embedding

seurat_obj[["lcl"]] <- Seurat::CreateDimReducObject(LCL_embedding,
                                                    assay = "RNA")

keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
names(keep_vec) <- Seurat::Cells(seurat_obj)
keep_vec[which(!is.na(LCL_embedding[,1]))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

Seurat::DefaultAssay(seurat_obj) <- "SCT"
hvg_vec <- Seurat::VariableFeatures(seurat_obj)
Seurat::DefaultAssay(seurat_obj) <- "RNA"
Seurat::VariableFeatures(seurat_obj) <- hvg_vec
seurat_obj <- Seurat::DietSeurat(seurat_obj,
                                 assays = "RNA",
                                 layers = "counts",
                                 features = Seurat::VariableFeatures(seurat_obj),
                                 dimreducs = c("python_X_diff", "python_X_tsne", "pca", "umap", "lcl"))

session_info <- devtools::session_info()
date_of_run <- Sys.time()

save(seurat_obj, 
     session_info, date_of_run,
     file = paste0(out_folder, "Writeup9b_celltag-data_with_LCL.RData"))

