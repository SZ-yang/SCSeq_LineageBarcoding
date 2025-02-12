rm(list=ls())

library(Seurat)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

seurat_obj <- NormalizeData(seurat_obj)
zz <- SeuratObject::LayerData(seurat_obj,
                              assay = "RNA",
                              layer = "data")
quantile(zz["FTL",])

idx1 <- which(seurat_obj$OG_condition == "cis")
idx2 <- which(seurat_obj$OG_condition == "cistocis")

gene_mat <- sapply(1:nrow(zz), function(j){
  if(j %% floor(nrow(zz)/10) == 0) cat('*')
  
  vec1 <- as.numeric(zz[j,idx1])
  vec2 <- as.numeric(zz[j,idx2])
  
  logfc <- mean(vec1) - mean(vec2)
  pvalue <- stats::wilcox.test(x = vec1, y = vec2)$p.value
  
  c(logfc = logfc, pvalue = pvalue)
})
gene_mat <- data.frame(t(gene_mat))
rownames(gene_mat) <- rownames(zz)

x_lim <- quantile(abs(gene_mat$logfc), probs = 0.9)

plot1 <- EnhancedVolcano::EnhancedVolcano(gene_mat,
                                          lab = rownames(gene_mat),
                                          x = 'logfc',
                                          y = 'pvalue',
                                          FCcutoff = x_lim)
ggplot2::ggsave(plot1,
                file = paste0(plot_folder, "Writeup7_de_cis-vs-cistocis.png"),
                height = 8, width = 8)
