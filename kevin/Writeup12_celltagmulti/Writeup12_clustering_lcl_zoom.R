rm(list=ls())

library(Seurat)
library(mclust)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup12_joshua-celltagmulti/cell_tag_integration.RData")

seurat_obj <- subset(integrated, batch == "tag")

lcl_csv <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-lcl/test_proj_embe_w_clone_id.csv")
rownames(lcl_csv) <- lcl_csv[,"X"]
lcl_csv <- lcl_csv[,grep("embed_", colnames(lcl_csv))]
lcl_csv <- lcl_csv[Seurat::Cells(seurat_obj),]
colnames(lcl_csv) <- paste0("LCL_", 1:ncol(lcl_csv))
lcl_csv <- as.matrix(lcl_csv)

umap_csv <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-lcl/cellTag-cellTagMulti_UMAP_test.csv")
rownames(umap_csv) <- umap_csv$X
umap_csv <- umap_csv[,c("UMAP0", "UMAP1")]
colnames(umap_csv) <- paste0("lclumap_", 1:2)
umap_csv <- umap_csv[Seurat::Cells(seurat_obj),]
umap_csv <- as.matrix(umap_csv)

seurat_obj[["LCL"]] <- Seurat::CreateDimReducObject(lcl_csv,
                                                    assay = "integrated")
seurat_obj[["LCLUMAP"]] <- Seurat::CreateDimReducObject(umap_csv,
                                                        assay = "integrated")

DimPlot(seurat_obj,
        reduction = "LCLUMAP")

seurat_obj$clone_id <- factor(paste0("Lineage:", seurat_obj$clone_id))

########

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    dims = 1:2, 
                                    reduction = "LCLUMAP",
                                    k.param = 5)
seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.01,
                                   group.singletons = FALSE)

DimPlot(seurat_obj,
        reduction = "LCLUMAP",
        group.by = "seurat_clusters",
        label = TRUE)

###################################

seurat_obj2 <- subset(seurat_obj, seurat_clusters %in% c(2,13,6,9))
DimPlot(seurat_obj2,
        reduction = "LCLUMAP",
        group.by = "seurat_clusters")

vec <- table(seurat_obj$clone_id)
vec2 <- table(seurat_obj2$clone_id)

plot(as.numeric(vec),
     as.numeric(vec2), asp = TRUE)
lines(c(0,1e5),
      c(0,1e5),
      col = 2)

###################################

seurat_obj2 <- subset(seurat_obj, seurat_clusters %in% c(1))
DimPlot(seurat_obj2,
        reduction = "LCLUMAP",
        group.by = "seurat_clusters")

vec <- table(seurat_obj$clone_id)
vec2 <- table(seurat_obj2$clone_id)

plot(as.numeric(vec),
     as.numeric(vec2), asp = TRUE)
lines(c(0,1e5),
      c(0,1e5),
      col = 2)
