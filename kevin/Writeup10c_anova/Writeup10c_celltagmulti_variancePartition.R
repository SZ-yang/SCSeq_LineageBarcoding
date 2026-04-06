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

tmp <- as.character(seurat_obj$predicted.id_cca_co)
tmp2 <- sapply(tmp, function(x){
  strsplit(x, split = "_")[[1]][1]
})
names(tmp2) <- NULL
seurat_obj$celltype <- as.factor(tmp2)
seurat_obj$predicted.id_cca_co <- as.factor(seurat_obj$predicted.id_cca_co)


#############################

seurat_obj <- Seurat::DietSeurat(seurat_obj, 
                                 layers = "counts")
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
scaled_data <- SeuratObject::LayerData(seurat_obj,
                                       assay = "RNA",
                                       layer = "scale.data")
var_data <- scaled_data
var_info <- seurat_obj@meta.data

# remove genes that are too sparse
## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
.nonzero_col <- function(mat, # is a cell-by-cell sparse matrix
                         col_idx, # this is an integer for the cell you're interested in
                         bool_value, # this is TRUE if you want to know the values, or FALSE if you just want to know the indices
                         bool_exclude_self = FALSE){
  stopifnot(col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  if(bool_exclude_self) stopifnot(nrow(mat) == ncol(mat))
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  
  # sometimes the value of "0" is accidentally stored in count_mat
  idx <- mat@i[(val1+1):val2]+1
  values <- mat@x[(val1+1):val2]
  idx <- idx[values != 0]
  values <- values[values != 0]
  if(length(values) == 0) return(numeric(0))
  
  if(bool_exclude_self){
    if(col_idx %in% idx){
      rm_idx <- which(col_idx %in% idx)
      idx <- idx[-rm_idx]
      values <- values[-rm_idx]
    }
  }
  
  if(bool_value){
    # return the value
    return(values)
  } else {
    # return the row index
    return(idx)
  }
}
count_data <- SeuratObject::LayerData(seurat_obj,
                                      assay = "RNA",
                                      layer = "counts",
                                      features = Seurat::VariableFeatures(seurat_obj))
count_data_t <- Matrix::t(count_data)
num_nonzero <- sapply(1:ncol(var_data_t), function(j){
  length(.nonzero_col(count_data_t,
               col_idx = j,
               bool_value = FALSE,
               bool_exclude_self = FALSE))/nrow(var_data_t)
})
names(num_nonzero) <- colnames(count_data_t)
var_data <- var_data[names(num_nonzero)[num_nonzero > 0.1],]

# fit
form <- ~ (1 | predicted.id_cca_co) + (1 | assigned_lineage)  + (1 | sample)
varPart <- variancePartition::fitExtractVarPartModel(var_data, 
                                                     form, 
                                                     var_info, 
                                                     showWarnings = FALSE)

save(varPart,
     file = "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/Writeup10c_celltagmulti_variancePartition.RData")

variance_explained <- Matrix::rowSums(varPart[,c("assigned_lineage", "predicted.id_cca_co", "sample")])
varPart <- varPart[order(variance_explained, decreasing = TRUE),]

round(varPart[1:50,],2)

print("Done! :)")


