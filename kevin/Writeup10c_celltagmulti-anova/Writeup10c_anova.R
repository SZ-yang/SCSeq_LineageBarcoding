rm(list=ls())

library(Seurat)

source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup10c_celltagmulti-anova/anova_functions.R")
load("~/kzlinlab/data/celltagging-multi_fibroblast/celltagging-multi_fibroblast.RData")
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/"

head(seurat_obj@meta.data)

keep_vec <- !is.na(seurat_obj$assigned_lineage)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

tab_vec <- table(seurat_obj$assigned_lineage)
lineage_names <- names(tab_vec)[which(tab_vec >= 50)]

keep_vec <- seurat_obj$assigned_lineage %in% lineage_names
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

data_mat <- SeuratObject::LayerData(
  seurat_obj,
  layer = "data",
  assay = "RNA",
  features = Seurat::VariableFeatures(seurat_obj)
)
depth_vec <- Matrix::rowSums(data_mat)
if(any(depth_vec <= 1e-6)){
  idx <- which(depth_vec >= 1e-6)
  data_mat <- data_mat[idx,]
}
data_mat <- as.matrix(data_mat)

covariates <- c("predicted.id_cca_co", "sample", "replicate", "assigned_lineage", "Phase", "md_fate_rev1", "RNA_snn_res.0.2")
metadata <- seurat_obj@meta.data[,covariates]
for(j in 1:ncol(metadata)){
  metadata[,j] <- droplevels(factor(metadata[,j]))
}

##############

anova_mat <- matrix(NA, nrow = nrow(data_mat), ncol = ncol(metadata))
rownames(anova_mat) <- rownames(data_mat)
colnames(anova_mat) <- colnames(metadata)

p <- nrow(data_mat)
k <- ncol(metadata)

for(j in 1:p){
  if(p > 10 && j %% floor(p/10) == 0) {
    cat('*')
    
    save(anova_mat, 
         file = paste0(out_folder, "Writeup10c_celltagmulti_anova_tmp.RData"))
  }
  
  for(i in 1:k){
    gene <- rownames(data_mat)[j]
    variable <- colnames(metadata)[i]
    
    x <- metadata[,variable]
    y <- data_mat[j,]
    
    anova_mat[gene,variable] <- .anova_percentage(x,y)
  }
}

save(anova_mat, 
     file = paste0(out_folder, "Writeup10c_celltagmulti_anova.RData"))
