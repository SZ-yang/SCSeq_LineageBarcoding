rm(list=ls())

library(Seurat)
library(mclust)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup12_joshua-celltagmulti/cell_tag_integration.RData")

seurat_obj <- subset(integrated, batch == "tag")
table(seurat_obj$assigned_lineage)
table(seurat_obj$clone_id)

seurat_obj$clone_id <- factor(paste0("Lineage:", seurat_obj$clone_id))

seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1)

DimPlot(seurat_obj,
        reduction = "umap",
        group.by = "seurat_clusters")

head(seurat_obj@meta.data[,c("clone_id", "seurat_clusters")])
summary(seurat_obj@meta.data[,c("clone_id", "seurat_clusters")])

md <- seurat_obj@meta.data
ari <- mclust::adjustedRandIndex(md$clone_id, md$seurat_clusters)
ari
