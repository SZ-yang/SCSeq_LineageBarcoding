rm(list=ls())
library(Seurat)
library(SeuratDisk)

# — your folders —
data_folder <- "/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/"
out_folder  <- data_folder

# — load your combined Seurat object —
load(file.path(data_folder, "cell_tag_integration.RData"))


# 1) Save your Seurat object to an H5Seurat file
SaveH5Seurat(
  integrated, 
  filename = "cell_tag_multi_integrated_seurat.h5seura", 
  overwrite = TRUE
)

# 2) Convert that H5Seurat into an AnnData (.h5ad)
Convert(
  "cell_tag_multi_integrated_seurat.h5seura", 
  dest      = "h5ad", 
  overwrite = TRUE
)
