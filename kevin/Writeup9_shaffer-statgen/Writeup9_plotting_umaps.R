rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup9/"

load(paste0(out_folder, "adata_with_lcl.RData"))

train_embedding <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup9/train_adata_LCL.csv")
rownames(train_embedding) <- train_embedding[,1]
train_embedding <- as.matrix(train_embedding[,-1])

test_embedding <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup9/test_adata_LCL.csv")
rownames(test_embedding) <- test_embedding[,1]
test_embedding <- as.matrix(test_embedding[,-1])

table(rownames(train_embedding) %in% Seurat::Cells(seurat_obj))
table(rownames(test_embedding) %in% Seurat::Cells(seurat_obj))

LCL_embedding <- matrix(NA, nrow = length(Seurat::Cells(seurat_obj)), ncol = ncol(train_embedding))
rownames(LCL_embedding) <- Seurat::Cells(seurat_obj)
colnames(LCL_embedding) <- paste0("LCL_", 1:ncol(LCL_embedding))
LCL_embedding[rownames(train_embedding),] <- train_embedding
LCL_embedding[rownames(test_embedding),] <- test_embedding

seurat_obj[["lcl"]] <- Seurat::CreateDimReducObject(LCL_embedding)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "lcl", 
                              dims = 1:ncol(LCL_embedding),
                              assay = 'RNA', 
                              reduction.name = 'lcl.umap', 
                              reduction.key = 'lclUMAP_')

########

Seurat::Idents(seurat_obj) <- "OG_condition"
embedding_list <- list(
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
    cell_idx1 <- which(!seurat_obj$OG_condition %in% conditions)
    cell_idx2 <- sample(which(seurat_obj$OG_condition %in% conditions))
      
    df <- df[c(cell_idx1, cell_idx2),]
    
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
                    filename = paste0(plot_folder, "Writeup9_lcl_embedding-", embedding_name, "_first-", treatment, ".png"),
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
    cell_idx1 <- which(!seurat_obj$OG_condition %in% conditions)
    cell_idx2 <- sample(which(seurat_obj$OG_condition %in% conditions))
    
    df <- df[c(cell_idx1, cell_idx2),]
    
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
                    filename = paste0(plot_folder, "Writeup9_lcl_embedding-", embedding_name, "_last-", treatment, ".png"),
                    height = 5, width = 5)
  }
}
