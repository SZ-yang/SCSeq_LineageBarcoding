rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

#################

tab_mat <- table(seurat_obj$Lineage, seurat_obj$RNA_snn_res.0.1)
rowsum_vec <- rowSums(tab_mat)
tab_mat <- tab_mat[order(rowsum_vec, decreasing = TRUE),]
colsum_vec <- colSums(tab_mat)
tab_mat <- tab_mat[,order(colsum_vec, decreasing = TRUE)]

rowsum_vec <- rowSums(tab_mat)
colsum_vec <- colSums(tab_mat)
num_lineages_per_cluster <- apply(tab_mat, 2, function(x){length(which(x!=0))})

df <- data.frame(
  cluster = colnames(tab_mat),
  num_cells = colsum_vec,
  num_lineages = num_lineages_per_cluster
)


threshold_cells <- 400
threshold_lineage <- 10

# Example scatter plot
plot1 <- ggplot(df, aes(x = num_lineages, y = num_cells,label = cluster)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter points
  geom_text_repel(data = df[df$num_cells > threshold_cells | df$num_lineages > threshold_lineage, ], # Adjust threshold
                  aes(label = cluster, color = "red"),
                  box.padding = 0.5, 
                  point.padding = 0.3) +  # Prevent label overlapping
  labs(x = "Number of Lineages", y = "Number of Cells", title = "Cluster Scatter Plot") +
  Seurat::NoLegend()
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lineage-by-cells_per-lcl-cluster.png"),
                height = 5, width = 5)

##################

# find all the cluster 0's
cluster_name <- "Cluster:11"
idx <- which(seurat_obj$RNA_snn_res.0.1 == cluster_name)
lineage_names <- unique(seurat_obj$Lineage[idx])
percentage_mat <- sapply(lineage_names, function(lineage){
  tmp_idx <- which(seurat_obj$Lineage == lineage)
  in_number <- which(seurat_obj$RNA_snn_res.0.1[tmp_idx] == cluster_name)
  c(in_ratio = length(in_number)/length(tmp_idx), 
    total_size = length(tmp_idx))
})
percentage_mat <- t(percentage_mat)
rownames(percentage_mat) <- lineage_names
percentage_mat <- percentage_mat[order(percentage_mat[,"total_size"], decreasing = TRUE),]
head(percentage_mat)
tail(percentage_mat)
plot(percentage_mat[,1],
     log10(percentage_mat[,2]),
     xlab = "In ratio",
     ylab = "Total size of lineage",
     pch = 16)

# now, for all the cells in the "pure" lineages, see what conditions they are in
lineage_names <- rownames(percentage_mat)[which(percentage_mat[,"in_ratio"] > 0.8)]
tmp_df <- seurat_obj@meta.data[,c("OG_condition", "Lineage")]
tmp_df <- tmp_df[which(tmp_df$Lineage %in% lineage_names),]
tab_mat <- table(tmp_df$Lineage, tmp_df$OG_condition)
tab_mat <- tab_mat[order(rowSums(tab_mat), decreasing = TRUE),]
tmp <- colSums(tab_mat)

original_percentage <- table(seurat_obj$OG_condition)
cluster_percentage <- rep(0, length(original_percentage))
names(cluster_percentage) <- names(original_percentage)
cluster_percentage[names(tmp)] <- tmp

round(100*cluster_percentage/original_percentage,2)

#################

hotspot_csv <- read.csv(paste0(out_folder, "Writeup7_LCL_hotspot_autocorrelations.csv"))
rownames(hotspot_csv) <- hotspot_csv$Gene
min_pval <- min(hotspot_csv$Pval[hotspot_csv$Pval > 0])
hotspot_csv$Pval <- pmax(hotspot_csv$Pval, min_pval)
hist(-log10(hotspot_csv$Pval))

idx <- which(hotspot_csv$C >= 0.2)

writeLines(hotspot_csv[idx, "Gene"], paste0(out_folder, "hotspot_gene_list.txt"))

hallmark_csv <- readxl::read_excel(paste0(out_folder, "41586_2023_6130_MOESM6_ESM.xlsx"),
                                   sheet = "Cancer MPs")
hallmark_csv <- as.data.frame(hallmark_csv)
for(j in 1:ncol(hallmark_csv)){
  vec <- hallmark_csv[,j]
  vec <- intersect(vec, rownames(hotspot_csv))
  hist(-log10(hotspot_csv$Pval),
       main = paste0("Column ", j, " (", colnames(hallmark_csv)[j], "),",
                     "\n# genes: ", length(vec)))
  idx <- which(rownames(hotspot_csv) %in% vec)
  rug(-log10(hotspot_csv$Pval[idx]), 
      col = 2,
      lwd = 2)
}

