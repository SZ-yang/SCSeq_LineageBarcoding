rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

tab_mat <- table(seurat_obj$Lineage, seurat_obj$fine_cluster)
# filter out lineages too small
rowsum_vec <- rowSums(tab_mat)
idx <- which(rowsum_vec >= 20)
tab_mat <- tab_mat[idx,]
colsum_vec <- colSums(tab_mat)
tab_mat <- tab_mat[,colsum_vec > 0]

# find a cluster with a ton of lineages
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
                filename = paste0(plot_folder, "Writeup7_lineage-by-cells_per-lcl-finecluster.png"),
                height = 4, width = 4)

# zoom in
# finecluster <- df$cluster[which.max(df$num_lineages)]
finecluster <- "LCL:285"
lineages <- rownames(tab_mat)[which(tab_mat[,finecluster]!=0)]

cell_idx <- which(seurat_obj$Lineage %in% lineages)
tab_mat_small <- table(droplevels(seurat_obj$Lineage[cell_idx]),
                       droplevels(seurat_obj$OG_condition[cell_idx]))
colnames(tab_mat_small) <- c("Cis", "Cis-to-CoCl2", "Cis-to-DabTram")
tab_mat_small

#######




# Load necessary library
library(ggplot2)
library(reshape2)

# Reshape the data to long format for ggplot2
df_long <- melt(tab_mat_small, 
                id.vars = "Lineage", 
                variable.name = "Treatment", 
                value.name = "Count")
colnames(df_long) <- c("Lineage", "Treatment", "Count")

# Create the bar plot
plot1 <- ggplot(df_long, aes(x = Treatment, y = Count, fill = Lineage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sequential treatment", y = "Number of surviving cells in lineage") +
  theme_minimal() +
  scale_fill_manual(values = c("#1b9e77", "#d95f02"))  # Custom colors

ggplot2::ggsave(plot1,
                file = paste0(plot_folder, "Writeup7_lcl_super-cluster_", finecluster, "_by-OG_condition.png"),
                height = 2.5, width = 5)

##################

tab_mat <- table(seurat_obj$Lineage, seurat_obj$OG_condition)
other_vec <- c("cistocis", "cocl2", "cocl2tocis",
               "cocl2tococl2", "cocl2todabtram", "dabtram",
               "dabtramtocis", "dabtramtococl2", "dabtramtodabtram")
bool_idx <- which(apply(tab_mat, 1, function(x){
  bool1 <- all(x[other_vec] == 0)
  bool2 <- x["cis"] > 1
  bool3 <- x["cistococl2"] > 2*x["cis"]
  bool4 <- x["cistodabtram"] > 0
  
  all(c(bool1, bool2, bool3, bool4))
}))


color_vec <- c(main = rgb(194, 155, 70, maxColorValue = 255),
               secondary = "black",
               other = "gray")
embedding <- seurat_obj[["lcl.umap"]]@cell.embeddings

lineage_vec <- as.character(seurat_obj$Lineage)
lineage_tmp <- lineage_vec
lineage_tmp[!lineage_vec %in% names(bool_idx)] <- "other"
lineage_tmp[lineage_vec %in% names(bool_idx)] <- "secondary"
lineage_tmp[lineage_vec %in% lineages] <- "main"

df <- cbind(data.frame(embedding), 
            factor(lineage_tmp))
colnames(df) <- c("x", "y", "lineage")

# shuffle indicies
cell_idx <- which(df$lineage == "other")
df <- df[c(cell_idx, setdiff(1:nrow(df), cell_idx)),]

plot1 <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y, ))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(
  color = lineage,
  size  = ifelse(lineage != "other", 2, 1)
))
plot1 <- plot1 + ggplot2::scale_colour_manual(values = color_vec)
plot1 <- plot1 + ggplot2::scale_size_identity() 
plot1 <- plot1 + cowplot::theme_cowplot()
plot1 <- plot1 + ggplot2::labs(x = "", y = "")
plot1 <- plot1 + Seurat::NoLegend() 

ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_super-cluster_", finecluster, "_UMAP.png"),
                height = 3, width = 5)

