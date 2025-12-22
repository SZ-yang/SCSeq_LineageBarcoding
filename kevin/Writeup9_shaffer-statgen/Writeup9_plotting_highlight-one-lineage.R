rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup9b/"

load(paste0(out_folder, "adata_with_lcl.RData"))

tab_mat <- table(seurat_obj$clone_id, seurat_obj$OG_condition)
dabtram_idx <- grep("^dabtram", colnames(tab_mat))
dabtram_count <- rowSums(tab_mat[,dabtram_idx])
chosen_lineage <- rownames(tab_mat)[which.max(dabtram_count)]

highlight_cells <- which(seurat_obj$clone_id == as.numeric(chosen_lineage))

########

# Plot just: dabtramtocis and cistodabtram
kk <- 1
embedding <- embedding_list[[kk]]
embedding_name <- names(embedding_list)[kk]
conditions <- c("dabtramtocis", "cistodabtram")

order_vec <- levels(seurat_obj$OG_condition)
order_vec <- c(setdiff(order_vec, conditions), conditions)

color_vec <- seurat_obj@misc[["OG_condition_colors"]]
color_vec[!names(color_vec) %in% conditions] <- "gray"

df <- cbind(data.frame(embedding), 
            seurat_obj$OG_condition)
colnames(df) <- c("x", "y", "OG_condition")

# shuffle indicies
cell_idx <- unlist(lapply(order_vec, function(grouping){
  which(seurat_obj$OG_condition == grouping)
}))
df <- df[cell_idx,]

# shuffle the conditions
cell_idx <- which(df$OG_condition %in% conditions)
df[cell_idx,] <- df[sample(cell_idx),]

plot1 <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y, ))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(
  color = OG_condition,
  size  = ifelse(OG_condition %in% conditions, 2, 1)
))
plot1 <- plot1 + ggplot2::scale_colour_manual(values = color_vec)
plot1 <- plot1 + ggplot2::scale_size_identity() 
plot1 <- plot1 + cowplot::theme_cowplot()
plot1 <- plot1 + ggplot2::labs(x = "", y = "")
plot1 <- plot1 + Seurat::NoLegend() 

ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl_embedding-", embedding_name, "_dabtram-cis.png"),
                height = 5, width = 5)
