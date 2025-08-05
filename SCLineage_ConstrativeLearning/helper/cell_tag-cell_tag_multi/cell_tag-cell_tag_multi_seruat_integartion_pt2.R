rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(scCustomize)

SeuratDisk::Convert("/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/Seurat_Integrate/adata_simplified.h5ad", 
                    dest = "h5seurat", 
                    overwrite = TRUE)

# Load the converted file as a Seurat object
seurat_obj <- SeuratDisk::LoadH5Seurat("/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/Seurat_Integrate/adata_simplified.h5seurat")
seurat_obj

# Now put the raw counts back in
raw_mat <- Seurat::ReadMtx(
  mtx = "/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/Seurat_Integrate/matrix.mtx",
  features = "/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/Seurat_Integrate/features.tsv",
  cells = "/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/Seurat_Integrate/barcodes.tsv"
)

# Make sure the raw counts has the same features and barcodes
raw_mat <- raw_mat[SeuratObject::Features(seurat_obj),
                   Seurat::Cells(seurat_obj)]

SeuratObject::LayerData(seurat_obj,
                        layer = "counts",
                        assay = "RNA") <- raw_mat

# Put in the metadata
metadata <- read.csv("/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/Seurat_Integrate/adata_obs.csv")
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, metadata)

# Since there's already an scVI and UMAP, we can directly visualize the data
cell_tag_multi <- seurat_obj
save(cell_tag_multi, 
     file = "/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/cell_tag-cell_tag_multi.RData")
