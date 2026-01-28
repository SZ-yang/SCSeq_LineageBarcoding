rm(list=ls())

library(Seurat)
library(scCustomize)


load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup12_joshua-celltagmulti/cell_tag_integration.RData")

seurat_obj <- subset(integrated, batch == "tag")

gene_vec <- c("Sfrp1", "Ptn", "Dcn", "Shh")
for(gene in gene_vec){
  print(scCustomize::FeaturePlot_scCustom(seurat_obj,
                                    features = gene,
                                    reduction = "umap"))
}

# Dcn,Ptn,Sfrp1

##########

ig_df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-IG/IG_analysis_res.csv")
rownames(ig_df) <- ig_df$gene
ig_df <- ig_df[order(ig_df$score_001_15, decreasing = TRUE),]

ig_df[gene_vec,]
#########