rm(list=ls())

library(Seurat)
library(Matrix)
library(variancePartition)
library(lme4)

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

#############################

seurat_obj <- Seurat::DietSeurat(seurat_obj, 
                                 layers = "counts")
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
scaled_data <- SeuratObject::LayerData(seurat_obj,
                                       assay = "RNA",
                                       layer = "scale.data")
var_info <- seurat_obj@meta.data

# fit
var_data <- scaled_data
form <- ~ (1 | predicted.id_cca_co) + (1 | assigned_lineage)  + (1 | sample)
varPart <- variancePartition::fitExtractVarPartModel(var_data, 
                                                     form, 
                                                     var_info, 
                                                     showWarnings = FALSE)

save(varPart,
     file = "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/Writeup10c_celltagmulti_variancePartition.RData")

# variance_explained <- Matrix::rowSums(varPart[,c("assigned_lineage", "Cell.type.annotation", "Time.point")])
# varPart <- varPart[order(variance_explained, decreasing = TRUE),]
# variance_explained <- variance_explained[order(variance_explained, decreasing = TRUE),]
# 
# head(varPart[1:50,])
# quantile(varPart[1:50,"assigned_lineage"]/variance_explained[1:50])

print("Done! :)")


