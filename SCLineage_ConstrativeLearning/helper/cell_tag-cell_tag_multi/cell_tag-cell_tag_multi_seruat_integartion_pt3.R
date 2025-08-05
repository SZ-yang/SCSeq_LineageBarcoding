rm(list=ls())
library(Seurat)
library(ggplot2)

# — your folders —
data_folder <- "/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/"
out_folder  <- data_folder

# — load your combined Seurat object —
load(file.path(data_folder, "cell_tag-cell_tag_multi.RData"))

seurat_list <- SplitObject(seurat_obj, split.by = "batch")

# free up memory
rm(seurat_obj); gc()

# — standard per‐dataset preprocessing (only if not already done) —
for (nm in names(seurat_list)) {
  so <- seurat_list[[nm]]
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_list[[nm]] <- so
}

# — pick genes to integrate on & find anchors —
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
anchors  <- FindIntegrationAnchors(
  object.list    = seurat_list,
  anchor.features= features)

# — integrate —
integrated <- IntegrateData(
  anchorset = anchors
)

DefaultAssay(integrated) <- "integrated"

# 2) Scale, PCA, and UMAP on the integrated values
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

# UMAP:
p <- DimPlot(
  integrated,
  reduction = "umap",
  group.by  = "batch",  # or your actual meta column
  pt.size   = 0.4
)

# Add a title with ggplot2
p + ggtitle("Integrated UMAP — Colored by Dataset")

# — save your results —
save(anchors, integrated,
     file = file.path(out_folder, "cell_tag_integration.RData"))