rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"
source("kevin/Writeup7_dylan/welch_anova.R")

load(paste0(out_folder, "adata_with_lcl.RData"))

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    reduction = "lcl",
                                    dims = 1:64)

seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.01)

seurat_obj$Lineage <- factor(paste0("Lineage:", seurat_obj$Lineage))
seurat_obj$RNA_snn_res.0.01 <- factor(paste0("Cluster:", seurat_obj$RNA_snn_res.0.01))

seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = 'lcl', 
                              dims = 1:64, 
                              assay = 'RNA', 
                              reduction.name = 'lcl.umap', 
                              reduction.key = 'lclUMAP_')

######

tab_mat <- table(droplevels(seurat_obj$Lineage), 
                 droplevels(seurat_obj$RNA_snn_res.0.01))
cluster_sizes <- colSums(tab_mat)

# keep only the cells in clusters with more than 100 cells
cluster_names <- colnames(tab_mat)[which(cluster_sizes > 200)]
keep_vec <- (seurat_obj$RNA_snn_res.0.01 %in% cluster_names)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

Seurat::DimPlot(seurat_obj,
                group.by = "RNA_snn_res.0.01",
                reduction = "lcl.umap")

seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.01)
tab_vec <- table(droplevels(seurat_obj$RNA_snn_res.0.01))
cluster_names <- names(tab_vec)[which(tab_vec > 200)]
keep_vec <- (seurat_obj$RNA_snn_res.0.01 %in% cluster_names)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

Seurat::DimPlot(seurat_obj,
                group.by = "RNA_snn_res.0.01",
                reduction = "lcl.umap")

seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.01)
tab_vec <- table(droplevels(seurat_obj$RNA_snn_res.0.01))
cluster_names <- names(tab_vec)[which(tab_vec > 200)]
keep_vec <- (seurat_obj$RNA_snn_res.0.01 %in% cluster_names)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

Seurat::DimPlot(seurat_obj,
                group.by = "RNA_snn_res.0.01",
                reduction = "lcl.umap")

###########################


dat <- SeuratObject::LayerData(
  seurat_obj,
  layer = "data",
  assay = "RNA"
)

gene_mat <- sapply(1:nrow(dat), function(i){
  if(i %% floor(nrow(dat)/10) == 0) cat('*')
  
  x <- as.numeric(dat[i,])
  y <- droplevels(seurat_obj$RNA_snn_res.0.1)
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

rm_idx <- unique(c(which(is.na(gene_mat[,1])), which(is.na(gene_mat[,2]))))
if(length(rm_idx) > 0){
  gene_mat <- gene_mat[-rm_idx,]
}


hallmark_csv <- readxl::read_excel(paste0(out_folder, "41586_2023_6130_MOESM6_ESM.xlsx"),
                                   sheet = "Cancer MPs")
hallmark_csv <- as.data.frame(hallmark_csv)
for(j in 1:ncol(hallmark_csv)){
  vec <- hallmark_csv[,j]
  vec <- intersect(vec, rownames(gene_mat))
  hist(gene_mat[,2],
       main = paste0("Column ", j, " (", colnames(hallmark_csv)[j], "),",
                     "\n# genes: ", length(vec)))
  idx <- which(rownames(gene_mat) %in% vec)
  rug(gene_mat[idx,2], 
      col = 2,
      lwd = 2)
}

############################

library(org.Mm.eg.db)
library(clusterProfiler)

tmp <- gene_mat[,2]
tmp <- tmp[order(tmp, decreasing = TRUE)]

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

actin_genes <- gse_df["GO:0051015", "core_enrichment"]
actin_genes <- strsplit(actin_genes, split = "/")[[1]]

seurat_obj <- Seurat::AddModuleScore(
  seurat_obj,
  features = list(actin_genes),
  ctrl = 5,
  name = "Actin"
)

seurat_obj$Actin2 <- pmax(seurat_obj$Actin1, 0)

Seurat::FeaturePlot(seurat_obj,
                    features = "Actin2",
                    reduction = "lcl.umap")


