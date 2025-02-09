rm(list=ls())
library(SeuratDisk)
library(Seurat)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"

# load(paste0(out_folder, "Writeup7_shaffer_preprocessed.RData"))

# Convert the H5AD file to an H5Seurat file
Convert(paste0(out_folder, "adata_with_lcl.h5ad"), 
        dest = "h5seurat", 
        overwrite = TRUE)

# Load the converted file as a Seurat object
seurat_obj <- LoadH5Seurat(paste0(out_folder, "adata_with_lcl.h5seurat"))

metadata <- read.csv(paste0(out_folder, "adata_obs_backup.csv"))
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]
seurat_obj@meta.data <- metadata

save(seurat_obj,
     file = paste0(out_folder, "adata_with_lcl.RData"))
