rm(list=ls())

library(Seurat)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup12_joshua-celltagmulti/cell_tag_integration.RData")

umap_mat <- integrated[["umap"]]@cell.embeddings