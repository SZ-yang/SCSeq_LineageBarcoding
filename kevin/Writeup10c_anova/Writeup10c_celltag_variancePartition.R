rm(list=ls())

library(Seurat)
library(Matrix)
library(variancePartition)
library(lme4)

load("~/kzlinlab/data/biddy_2018_celltag/biddy_seurat.RData")
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/"

head(seurat_obj@meta.data)

keep_vec <- !is.na(seurat_obj$CellTagD0_85k)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

tab_vec <- table(seurat_obj$CellTagD0_85k)
lineage_names <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- seurat_obj$CellTagD0_85k %in% lineage_names
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

#############################
Seurat::DefaultAssay(seurat_obj) <- "SCT"
gene_vec <- Seurat::VariableFeatures(seurat_obj)

Seurat::DefaultAssay(seurat_obj) <- "RNA"
Seurat::VariableFeatures(seurat_obj) <- gene_vec

seurat_obj <- Seurat::DietSeurat(seurat_obj, 
                                 layers = "counts",
                                 features = Seurat::VariableFeatures(seurat_obj))
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
scaled_data <- SeuratObject::LayerData(seurat_obj,
                                       assay = "RNA",
                                       layer = "data")
seurat_obj$reprogramming_day <- factor(paste0("T", seurat_obj$reprogramming_day))
seurat_obj$CellTagD0_85k <- factor(paste0("Lin", seurat_obj$CellTagD0_85k))
seurat_obj$cell_type <- factor(seurat_obj$cell_type) 
var_info <- seurat_obj@meta.data

# fit
var_data <- scaled_data
form <- ~ (1 | cell_type) + (1 | CellTagD0_85k)  + (1 | reprogramming_day)
varPart <- variancePartition::fitExtractVarPartModel(var_data, 
                                                     form, 
                                                     var_info, 
                                                     showWarnings = FALSE)

save(varPart,
     file = "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/Writeup10c_celltag_variancePartition.RData")

# variance_explained <- Matrix::rowSums(varPart[,c("assigned_lineage", "Cell.type.annotation", "Time.point")])
# varPart <- varPart[order(variance_explained, decreasing = TRUE),]
# variance_explained <- variance_explained[order(variance_explained, decreasing = TRUE),]
# 
# head(varPart[1:50,])
# quantile(varPart[1:50,"assigned_lineage"]/variance_explained[1:50])

print("Done! :)")


