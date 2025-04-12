rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

color_vec <- c(Lineage1 = rgb(203, 50, 39, maxColorValue = 255),
               Lineage2 = rgb(65, 113, 170, maxColorValue = 255),
               Lineage3 = rgb(131, 72, 149, maxColorValue = 255),
               Lineage4 = rgb(92, 163, 77, maxColorValue = 255),
               Lineage5 = rgb(235, 123, 44, maxColorValue = 255),
               other = "gray")

embedding_list <- list(
  umap = seurat_obj[["umap"]]@cell.embeddings,
  lcl = seurat_obj[["lcl.umap"]]@cell.embeddings
)

lineage_vec <- seurat_obj$Lineage
tab_vec <- table(lineage_vec)
tab_vec <- sort(tab_vec, decreasing = TRUE)
top_lineages <- names(tab_vec)[1:5]

lineage_tmp <- as.character(lineage_vec)
lineage_tmp[!lineage_tmp %in% top_lineages] <- "other"
for(i in 1:length(top_lineages)){
  lineage <- top_lineages[i]
  lineage_tmp[lineage_tmp == lineage] <- paste0("Lineage", i)
}

seurat_obj$lineage_tmp <- lineage_tmp

for(kk in 1:length(embedding_list)){
  embedding <- embedding_list[[kk]]
  embedding_name <- names(embedding_list)[kk]
  
  df <- cbind(data.frame(embedding), 
              factor(seurat_obj$lineage_tmp))
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
                  filename = paste0(plot_folder, "Writeup7_lcl_embedding-", embedding_name, "_top-lineages.png"),
                  height = 5, width = 5)
}

