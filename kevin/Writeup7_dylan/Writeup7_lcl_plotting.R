rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

############

plot1 <- scCustomize::DimPlot_scCustom(seurat_obj,
                                       reduction = "lcl.umap",
                                       group.by = "OG_condition",
                                       colors_use = seurat_obj@misc[["OG_condition_colors"]])
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl-umap_OG_condition.png"),
                height = 5, width = 8)

################

plot1 <- scCustomize::DimPlot_scCustom(seurat_obj,
                                       reduction = "lcl.umap",
                                       group.by = "RNA_snn_res.0.01")
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup7_lcl-umap_RNA_snn_res.0.01.png"),
                height = 5, width = 5)

################

Seurat::Idents(seurat_obj) <- "OG_condition"
embedding_list <- list(
  umap = seurat_obj[["umap"]]@cell.embeddings,
  lcl = seurat_obj[["lcl.umap"]]@cell.embeddings
)

# now plot by first timepoint
treatment_vec <- c("cis", "cocl2", "dabtram")
for(treatment in treatment_vec){
  for(kk in 1:length(embedding_list)){
    embedding <- embedding_list[[kk]]
    embedding_name <- names(embedding_list)[kk]
    
    tmp <- grep(paste0("^", treatment), levels(seurat_obj$OG_condition))
    conditions <- levels(seurat_obj$OG_condition)[tmp]
    
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
                    filename = paste0(plot_folder, "Writeup7_lcl_embedding-", embedding_name, "_first-", treatment, ".png"),
                    height = 5, width = 5)
  }
}

####

# now plot by second timepoint
treatment_vec <- c("cis", "cocl2", "dabtram")
for(treatment in treatment_vec){
  for(kk in 1:length(embedding_list)){
    embedding <- embedding_list[[kk]]
    embedding_name <- names(embedding_list)[kk]
    tmp <- grep(paste0(treatment, "$"), levels(seurat_obj$OG_condition))
    conditions <- levels(seurat_obj$OG_condition)[tmp]
    
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
                    filename = paste0(plot_folder, "Writeup7_lcl_embedding-", embedding_name, "_last-", treatment, ".png"),
                    height = 5, width = 5)
  }
}



