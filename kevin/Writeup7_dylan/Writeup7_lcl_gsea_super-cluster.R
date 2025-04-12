rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"
source("kevin/Writeup7_dylan/welch_anova.R")

load(paste0(out_folder, "adata_with_lcl.RData"))

seurat_obj <- subset(seurat_obj, in_lcl_cluster == TRUE)

dat <- SeuratObject::LayerData(
  seurat_obj,
  layer = "data",
  assay = "RNA"
)

gene_mat <- sapply(1:nrow(dat), function(i){
  if(i %% floor(nrow(dat)/10) == 0) cat('*')
  
  x <- as.numeric(dat[i,])
  y <- droplevels(seurat_obj$fine_cluster)
  df <- data.frame(response_col = x,
                   group_col = y)
  
  df <- remove_zero_variance(data = df, 
                             response_col = "response_col", 
                             group_col = "group_col")
  
  res <- tryCatch({welch_anova(data = df, 
                               response_col = "response_col", 
                               group_col = "group_col")},
                  error = function(e){
                    list(p.value = NA, R2_welch = NA)
                  })
  
  c(res$p.value, res$R2_welch)
})
gene_mat <- t(gene_mat)
rownames(gene_mat) <- rownames(dat)
colnames(gene_mat) <- c("p.value", "R2_welch")

library(org.Mm.eg.db)
library(clusterProfiler)

tmp <- gene_mat[,2]
tmp <- tmp[order(tmp, decreasing = TRUE)]
all(!is.na(tmp))

set.seed(10)
gse <- clusterProfiler::gseGO(
  tmp,
  ont = "MF", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500,            # maximum gene set size
  scoreType = "pos"
)

gse_df <- as.data.frame(gse)

gse_df[ c("Description","p.adjust")]



rm_idx <- unique(c(which(is.na(gene_mat[,1])), which(is.na(gene_mat[,2]))))
if(length(rm_idx) > 0){
  gene_mat <- gene_mat[-rm_idx,]
}