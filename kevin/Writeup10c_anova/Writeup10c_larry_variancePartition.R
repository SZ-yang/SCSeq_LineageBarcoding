rm(list=ls())

library(Seurat)
library(Matrix)
library(variancePartition)
library(lme4)


source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup10c_celltagmulti-anova/anova_functions.R")
load("~/kzlinlab/data/larry_hematopoiesis/larry-dataset_KZL.RData")
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/"

# for the lineages w/ more than 50 cells, stratify by lineage AND cell-type. 

keep_vec <- !is.na(seurat_obj$assigned_lineage)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

keep_vec <- seurat_obj$Cell.type.annotation %in% c("Undifferentiated", "Neutrophil", "Monocyte")
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# for any lineage x cell-type with less than 5 cells, remove those cells
tab_mat <- table(seurat_obj$assigned_lineage, seurat_obj$Cell.type.annotation)

# ============================================================
# 1) Remove cells in any (lineage x cell-type) with < 3 cells
# ============================================================
md <- seurat_obj@meta.data
combo_key <- paste(md$assigned_lineage, md$Cell.type.annotation, sep = "||")
combo_tab <- table(combo_key)

keep_step1 <- combo_tab[combo_key] >= 3
cells_step1 <- rownames(md)[keep_step1]
seurat_obj <- subset(seurat_obj, cells = cells_step1)

# ============================================================
# 2) Remove any lineage whose remaining cells are only 1 cell-type
# ============================================================
md1 <- seurat_obj@meta.data

# number of distinct cell types per lineage (after step 1)
n_types_per_lineage <- tapply(
  md1$Cell.type.annotation,
  md1$assigned_lineage,
  function(x) length(unique(x[!is.na(x)]))
)

keep_lineages <- names(n_types_per_lineage)[n_types_per_lineage > 2]
cells_step2 <- rownames(md1)[md1$assigned_lineage %in% keep_lineages]
seurat_obj <- subset(seurat_obj, cells = cells_step2)

# look at the final meta.data
tab_mat <- table(seurat_obj$assigned_lineage, seurat_obj$Cell.type.annotation)

seurat_obj <- Seurat::DietSeurat(seurat_obj,
                                 features = Seurat::VariableFeatures(seurat_obj))

#############################

seurat_obj <- Seurat::DietSeurat(seurat_obj, 
                                         layers = "counts")
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
scaled_data <- SeuratObject::LayerData(seurat_obj,
                                       assay = "RNA",
                                       layer = "scale.data")
seurat_obj$Time.point <- factor(paste0("T", seurat_obj$Time.point))
var_info <- seurat_obj@meta.data

# fit
var_data <- scaled_data
form <- ~ (1 | Cell.type.annotation) + (1 | assigned_lineage)  + (1 | Time.point)
varPart <- variancePartition::fitExtractVarPartModel(var_data, 
                                                     form, 
                                                     var_info, 
                                                     showWarnings = FALSE)

save(varPart,
     file = "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/Writeup10c_larry_variancePartition.RData")

variance_explained <- Matrix::rowSums(varPart[,c("assigned_lineage", "Cell.type.annotation", "Time.point")])
varPart <- varPart[order(variance_explained, decreasing = TRUE),]
variance_explained <- variance_explained[order(variance_explained, decreasing = TRUE),]

head(varPart[1:50,])
quantile(varPart[1:50,"assigned_lineage"]/variance_explained[1:50])

#######################


