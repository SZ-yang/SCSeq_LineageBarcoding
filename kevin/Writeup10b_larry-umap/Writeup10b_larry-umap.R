rm(list=ls())

library(Seurat)

source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup10c_celltagmulti-anova/anova_functions.R")
load("~/kzlinlab/data/larry_hematopoiesis/larry-dataset_KZL.RData")
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/"
fig_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup10b/"

new_names <- with(seurat_obj@meta.data, paste(Library, Cell.barcode, sep = ":"))
new_names <- gsub("-", replacement = "", x = new_names)
seurat_obj <- RenameCells(seurat_obj, new.names = new_names)

# load in the LCL embedding
data_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10b/"
embedding <- read.csv(paste0(data_folder, "Larry_full_train_base_embedding.csv"),
                      row.names = 1)
rownames(embedding) <- embedding[,1]
embedding <- embedding[,-1]

keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[Seurat::Cells(seurat_obj) %in% rownames(embedding)] <- TRUE
table(keep_vec)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)
embedding <- embedding[Seurat::Cells(seurat_obj),]

embedding <- as.matrix(embedding)
colnames(embedding) <- paste0("LCL_", 1:ncol(embedding))
seurat_obj[["LCL"]] <- Seurat::CreateDimReducObject(embedding)

set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              reduction = "LCL",
                              dims = 1:ncol(embedding),
                              reduction.name = "LCLUMAP")


umap_mat <- seurat_obj[["LCLUMAP"]]@cell.embeddings

# assigned colors
lineage_vec <- seurat_obj$assigned_lineage
tab_vec <- table(lineage_vec)
top_lineages <- names(sort(tab_vec, decreasing = TRUE))[1:5]
col_vec <- rep("gray", nrow(umap_mat))
for(i in 1:length(top_lineages)){
  lineage <- top_lineages[i]
  idx <- which(lineage_vec == lineage)
  col_vec[idx] <- i+1
}

png(paste0(fig_folder, "Writeup10b_LARRY_UMAP.png"),
    height = 8, width = 8, units = "in", res = 300)
plot(umap_mat[,1],
       umap_mat[,2],
       col = col_vec,
       pch = 16,
       asp = TRUE)
graphics.off()






