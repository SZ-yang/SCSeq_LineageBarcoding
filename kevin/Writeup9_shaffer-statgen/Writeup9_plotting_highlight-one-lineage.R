rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup9/"

load(paste0(out_folder, "adata_with_lcl.RData"))

tab_mat <- table(seurat_obj$clone_id, seurat_obj$OG_condition)
dabtram_idx <- grep("^dabtram", colnames(tab_mat))
dabtram_count <- rowSums(tab_mat[,dabtram_idx])
chosen_lineage <- rownames(tab_mat)[which.max(dabtram_count)]
tab_mat[chosen_lineage,]

highlight_cells <- which(seurat_obj$clone_id == as.numeric(chosen_lineage))

########

embedding <- seurat_obj[["umap"]]@cell.embeddings
condition_vec <- as.character(seurat_obj$OG_condition)
condition_vec[which(seurat_obj$clone_id != as.numeric(chosen_lineage))] <- "other"
condition_vec <- factor(condition_vec, levels = c("other",  colnames(tab_mat)[grep("^dabtram", colnames(tab_mat))]))

color_vec <- seurat_obj@misc[["OG_condition_colors"]]
color_vec[!names(color_vec) %in% conditions] <- "gray"
color_vec$other <- "gray"

df <- cbind(data.frame(embedding), 
            condition_vec)
colnames(df) <- c("x", "y", "OG_condition")

# shuffle indicies
cell_idx <- unlist(lapply(levels(condition_vec), function(grouping){
  which(condition_vec == grouping)
}))
df <- df[cell_idx,]

plot1 <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y))
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
                filename = paste0(plot_folder, "Writeup9_umap-highlighting-one-lineage.png"),
                height = 5, width = 5)
