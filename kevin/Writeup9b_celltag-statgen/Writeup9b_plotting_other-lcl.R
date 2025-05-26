rm(list=ls())

library(Seurat)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup9b/"
load(paste0(out_folder, "Writeup9b_celltag-data_with_LCL.RData"))

train_embedding <- read.csv(paste0(out_folder, "train_adata_LCL_other-lcl.csv"))
rownames(train_embedding) <- train_embedding[,1]
train_embedding <- as.matrix(train_embedding[,-1])

test_embedding <- read.csv(paste0(out_folder, "test_adata_LCL_other-lcl.csv"))
rownames(test_embedding) <- test_embedding[,1]
test_embedding <- as.matrix(test_embedding[,-1])

table(rownames(train_embedding) %in% Seurat::Cells(seurat_obj))
table(rownames(test_embedding) %in% Seurat::Cells(seurat_obj))

LCL_embedding <- matrix(NA, nrow = length(Seurat::Cells(seurat_obj)), ncol = ncol(train_embedding))
rownames(LCL_embedding) <- Seurat::Cells(seurat_obj)
colnames(LCL_embedding) <- paste0("LCL_", 1:ncol(LCL_embedding))
LCL_embedding[rownames(train_embedding),] <- train_embedding
LCL_embedding[rownames(test_embedding),] <- test_embedding

seurat_obj[["lcl"]] <- Seurat::CreateDimReducObject(LCL_embedding,
                                                    assay = "RNA")

train_metadata <- read.csv(paste0(out_folder, "train_obs.csv"))
rownames(train_metadata) <- train_metadata[,1]

test_metadata <- read.csv(paste0(out_folder, "test_obs.csv"))
rownames(test_metadata) <- test_metadata[,1]

clone_vec <- rep(NA, length(Seurat::Cells(seurat_obj)))
names(clone_vec) <- Seurat::Cells(seurat_obj)
clone_vec[rownames(train_metadata)] <- train_metadata$clone_id
clone_vec[rownames(test_metadata)] <- test_metadata$clone_id

seurat_obj$clone_id <- clone_vec

testing_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
names(testing_vec) <- Seurat::Cells(seurat_obj)
testing_vec[rownames(test_metadata)] <- TRUE
seurat_obj$test <- testing_vec

set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "lcl", 
                              dims = 1:ncol(seurat_obj[["lcl"]]@cell.embeddings),
                              assay = 'RNA', 
                              reduction.name = 'lcl.umap', 
                              reduction.key = 'lclUMAP_')

scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "cell_type")

scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "test")


#############

color_vec <- c(Lineage1 = rgb(203, 50, 39, maxColorValue = 255),
               Lineage2 = rgb(65, 113, 170, maxColorValue = 255),
               Lineage3 = rgb(92, 163, 77, maxColorValue = 255),
               Lineage4 = rgb(131, 72, 149, maxColorValue = 255),
               Lineage5 = rgb(235, 123, 44, maxColorValue = 255),
               other = "gray80")

tab_vec <- table(seurat_obj$clone_id)
tab_vec <- sort(tab_vec, decreasing = TRUE)
largest_lineage_vec <- rep("other", length(Seurat::Cells(seurat_obj)))
names(largest_lineage_vec) <- Seurat::Cells(seurat_obj)
for(i in 1:5){
  lineage_name <- names(tab_vec)[i]
  largest_lineage_vec[which(seurat_obj$clone_id == lineage_name)] <- paste0("Lineage", i)
}
seurat_obj$largest_lineage <- largest_lineage_vec

Seurat::Idents(seurat_obj) <- "largest_lineage"
scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "lcl.umap", 
                              group.by = "largest_lineage", 
                              colors_use = color_vec,
                              order = paste0("Lineage", 1:5),
                              shuffle = TRUE)

Seurat::Idents(seurat_obj) <- "largest_lineage"
scCustomize::DimPlot_scCustom(seurat_obj, 
                              reduction = "python_X_tsne", 
                              group.by = "largest_lineage", 
                              colors_use = color_vec,
                              order = paste0("Lineage", 1:5))

#######################

embedding <- seurat_obj[["lcl.umap"]]@cell.embeddings
condition_vec <- factor(seurat_obj$largest_lineage)

df <- cbind(data.frame(embedding), 
            condition_vec)
colnames(df) <- c("x", "y", "Lineage")

# shuffle indicies
cell_idx <- which(condition_vec == "other")
df <- df[c(cell_idx, setdiff(1:nrow(df), cell_idx)),]

plot1 <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(
  color = Lineage,
  size  = ifelse(Lineage %in% paste0("Lineage", 1:5), 2, 0.5)
))
plot1 <- plot1 + ggplot2::scale_colour_manual(values = color_vec)
plot1 <- plot1 + ggplot2::scale_size_identity() 
plot1 <- plot1 + cowplot::theme_cowplot()
plot1 <- plot1 + ggplot2::labs(x = "", y = "")
plot1 <- plot1 + Seurat::NoLegend() 

ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup9_umap-highlighting-one-lineage.png"),
                height = 5, width = 5)






