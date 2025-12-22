rm(list=ls())

library(Seurat)

load("~/kzlinlab/data/celltagging-multi_fibroblast/celltagging-multi_fibroblast.RData")
ctM_metadata <- seurat_obj@meta.data

load("~/kzlinlab/data/biddy_2018_celltag/biddy_seurat.RData")
ct_metadata <- seurat_obj@meta.data

rm(list="seurat_obj")

