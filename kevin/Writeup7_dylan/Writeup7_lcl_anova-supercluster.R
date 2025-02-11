rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dbscan)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"
source("kevin/Writeup7_dylan/welch_anova.R")

load(paste0(out_folder, "adata_with_lcl.RData"))
seurat_obj$OG_condition <- factor(seurat_obj$OG_condition)
seurat_obj$Lineage <- factor(seurat_obj$Lineage)

color_vec <- c("cis" = rgb(202,111,55, maxColorValue = 255),
               "cistocis" = rgb(248,210,152, maxColorValue = 255),
               "cistococl2" = rgb(240,148,71, maxColorValue = 255),
               "cistodabtram" = rgb(160,61,38, maxColorValue = 255),
               "cocl2" = rgb(69,132,69, maxColorValue = 255),
               "cocl2tocis" = rgb(131,202,163, maxColorValue = 255),
               "cocl2tococl2" = rgb(126,191,90, maxColorValue = 255),
               "cocl2todabtram" = rgb(35,63,58, maxColorValue = 255),
               "dabtram" = rgb(68,49,147, maxColorValue = 255),
               "dabtramtocis" = rgb(147,137,193, maxColorValue = 255),
               "dabtramtococl2" = rgb(145,54,147, maxColorValue = 255),
               "dabtramtodabtram" = rgb(68,32,85, maxColorValue = 255))

embedding <- seurat_obj[["lcl"]]@cell.embeddings
apply(embedding, 2, quantile)

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    reduction = "lcl",
                                    dims = 1:64)

seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.01)

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

# keep only the cells in clusters with more than 200 cells
cluster_names <- colnames(tab_mat)[which(cluster_sizes > 200)]
keep_vec <- (seurat_obj$RNA_snn_res.0.01 %in% cluster_names)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)
plot1 <- Seurat::DimPlot(seurat_obj,
                group.by = "RNA_snn_res.0.01",
                reduction = "lcl.umap")
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_anova-refinement1.png"),
                height = 5, width = 8)

embedding <- seurat_obj[["lcl"]]@cell.embeddings
dbscan_res <- dbscan::dbscan(
  x = embedding,
  eps = 2,
  minPts = 50
)
dbscan_cluster <- dbscan_res$cluster
dbscan_cluster[dbscan_cluster == 0] <- NA
names(dbscan_cluster) <- Seurat::Cells(seurat_obj)
seurat_obj$dbscan <- dbscan_cluster

plot1 <- Seurat::DimPlot(seurat_obj,
                group.by = "dbscan",
                reduction = "lcl.umap")
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_anova-refinement2.png"),
                height = 5, width = 8)

seurat_obj$keep <- !is.na(seurat_obj$dbscan)
seurat_obj <- subset(seurat_obj, keep == TRUE)

seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.01)
tab_vec <- table(seurat_obj$RNA_snn_res.0.01)
cluster_names <- names(tab_vec)[which(tab_vec > 200)]
keep_vec <- (seurat_obj$RNA_snn_res.0.01 %in% cluster_names)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)
plot1 <- Seurat::DimPlot(seurat_obj,
                group.by = "RNA_snn_res.0.01",
                reduction = "lcl.umap")
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_anova-refinement3.png"),
                height = 5, width = 8)

seurat_obj$Lineage <- factor(paste0("Lineage:", seurat_obj$Lineage))
seurat_obj$Lineage <- droplevels(seurat_obj$Lineage)
seurat_obj$RNA_snn_res.0.01 <- factor(paste0("Cluster:", seurat_obj$RNA_snn_res.0.01))

plot1 <- scCustomize::DimPlot_scCustom(seurat_obj,
                              group.by = "Lineage",
                              reduction = "lcl.umap") + Seurat::NoLegend()
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_anova-refinement3_by-lineage.png"),
                height = 5, width = 8)

plot1 <- scCustomize::DimPlot_scCustom(seurat_obj,
                                       group.by = "OG_condition",
                                       reduction = "lcl.umap",
                                       colors_use = color_vec) + Seurat::NoLegend()
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_anova-refinement3_by-OG_condition.png"),
                height = 5, width = 8)

###########################

dat <- SeuratObject::LayerData(
  seurat_obj,
  layer = "data",
  assay = "RNA"
)

gene_mat <- sapply(1:nrow(dat), function(i){
  if(i %% floor(nrow(dat)/10) == 0) cat('*')
  
  x <- as.numeric(dat[i,])
  y <- droplevels(seurat_obj$RNA_snn_res.0.01)
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

pdf(paste0(plot_folder, "Writeup7_lcl_anova_hallmark.pdf"),
    height = 5, width = 8, onefile = TRUE)
for(j in 1:ncol(hallmark_csv)){
  vec <- hallmark_csv[,j]
  vec <- intersect(vec, rownames(gene_mat))
  hist(gene_mat[,2],
       main = paste0("Column ", j, " (", colnames(hallmark_csv)[j], "),",
                     "\n# genes: ", length(vec)),
       xlab = "Percent explained by LCL cluster (per gene)")
  idx <- which(rownames(gene_mat) %in% vec)
  rug(gene_mat[idx,2], 
      col = 2,
      lwd = 2)
}
graphics.off()

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
  ctrl = 50,
  name = "Actin"
)

seurat_obj$Actin2 <- pmax(seurat_obj$Actin1, 0)

plot1 <- scCustomize::FeaturePlot_scCustom(seurat_obj,
                                  features = "Actin2",
                                  reduction = "lcl.umap")
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_anova_Actin_by-lcl.png"),
                height = 5, width = 8)

plot1 <- scCustomize::FeaturePlot_scCustom(seurat_obj,
                                  features = "Actin2",
                                  reduction = "umap")
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_anova_Actin_by-umap.png"),
                height = 5, width = 8)

##########

# it really comes down to Cluster 0...
tab_mat <- table(seurat_obj$Lineage, 
                 seurat_obj$RNA_snn_res.0.01)

# number of lineages per cluster
num_lineage_per_cluster <- apply(tab_mat, 2, function(x){length(which(x!=0))})

# let's look at one with 3 lineages
cluster_name <- names(num_lineage_per_cluster)[which(num_lineage_per_cluster == 3)[2]]
lineage_names <- rownames(tab_mat)[which(tab_mat[,cluster_name] != 0)]

idx <- which(seurat_obj$Lineage %in% lineage_names)
tab_mat2 <- table(droplevels(seurat_obj$Lineage[idx]),
                  droplevels(seurat_obj$OG_condition[idx]))
tab_mat2

##############

# among the remaining lineages, see if any are similar
tab_mat <- table(droplevels(seurat_obj$Lineage), 
                 droplevels(seurat_obj$OG_condition))
rowsum_vec <- rowSums(tab_mat)
tab_mat2 <- tab_mat
# normalize
for(i in 1:nrow(tab_mat2)){
  tab_mat2[i,] <- tab_mat2[i,]/sum(tab_mat2[i,])
}
n <- nrow(tab_mat2)
cor_mat <- stats::cor(t(tab_mat2))
diag(cor_mat) <- NA
hist(as.numeric(cor_mat))

idx <- which(cor_mat > 0.9, arr.ind = TRUE)
unique_lineage_idx <- unique(as.numeric(idx))
largest_lineage_idx <- unique_lineage_idx[which.max(rowsum_vec[unique_lineage_idx])]
rowsum_vec[largest_lineage_idx]

# zoom in on this lineage
partner_lineage_idx <- which(cor_mat[largest_lineage_idx,] >= 0.9)
rowsum_vec[partner_lineage_idx]

lineage_idx <- c(largest_lineage_idx, partner_lineage_idx)
cor_mat[lineage_idx, lineage_idx]
tab_mat[lineage_idx,]
tab_mat2[lineage_idx,]

lineage_names <- rownames(cor_mat)[lineage_idx]
cell_idx <- which(seurat_obj$Lineage %in% lineage_names)
tab_mat <- table(droplevels(seurat_obj$Lineage[cell_idx]),
                 droplevels(seurat_obj$RNA_snn_res.0.01[cell_idx]))
