rm(list=ls())

library(Seurat)
library(scCustomize)


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

##########

ig_df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-IG/IG_analysis_res.csv")
rownames(ig_df) <- ig_df$gene
ig_df <- ig_df[order(ig_df$score_001_15, decreasing = TRUE),]

#########

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    dims = 1:31, 
                                    reduction = "LCL",
                                    k.param = 5)
seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.1,
                                   group.singletons = FALSE)

Seurat::FeaturePlot(seurat_obj,
                    features = "Apoa1",
                    reduction = "LCLUMAP")

DimPlot(seurat_obj,
        reduction = "LCLUMAP",
        group.by = "seurat_clusters",
        label = TRUE)

de_res <- Seurat::FindAllMarkers(seurat_obj, group.by = "seurat_clusters")
de_res2 <- de_res
de_res2 <- de_res2[which(de_res2$avg_log2FC >= 3),]
de_res2 <- de_res2[which(de_res2$pct.1 >= 0.3),]
de_res2 <- de_res2[which(de_res2$p_val_adj <= 0.05),]
de_res2 <- de_res2[order(abs(de_res2$avg_log2FC), decreasing = TRUE),]
rownames(de_res2) <- de_res2$gene

gene_vec <- c("Ttr", "Acox2", "Bhlha15", "Shh", "Sfrp1", "Ptn", "Dcn", "Cyp2f2")
for(gene in gene_vec){
  print(scCustomize::FeaturePlot_scCustom(seurat_obj,
                                    features = gene,
                                    reduction = "LCLUMAP"))
}

# Dcn,Ptn,Sfrp1, maybe Shh
