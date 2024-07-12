rm(list=ls()); gc(TRUE)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v5")

# from https://github.com/theislab/zellkonverter/issues/38
adata <- zellkonverter::readH5AD(
  file = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup3/Writeup3_Larry_scBaseEncoderFeat_Z_bs30_tau05.h5ad"
)

names(adata@assays)

tmp <- Seurat::as.Seurat(adata, counts = "X", data = "X")
# converted most of the things. It creates an assay by default called "originalexp"
# see also https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html

seurat_obj <- Seurat::CreateSeuratObject(
  counts = tmp[["originalexp"]],
  data = tmp[["originalexp"]],
  meta.data = tmp@meta.data
)

seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

# put in the assays
name_vec <- names(adata@assays)
gene_vec <- SeuratObject::Features(seurat_obj[["RNA"]])

# put in the gene metafeatures
gene_metadata <- SingleCellExperiment::rowData(adata)
seurat_obj[["RNA"]]@misc <- as.data.frame(gene_metadata)

# put in the dimension reductions
name_vec <- SingleCellExperiment::reducedDimNames(adata)
for(name_val in name_vec){
  mat <- SingleCellExperiment::reducedDim(adata, name_val)
  name_val2 <- paste0("python_", name_val)
  colnames(mat) <- paste0(name_val2, "_", 1:ncol(mat))
  
  seurat_obj[[name_val2]] <- Seurat::CreateDimReducObject(embeddings = mat,
                                                          assay = "RNA")
}

# put in the metadata
metadata_list <- adata@metadata
idx <- which(sapply(1:length(metadata_list), function(i){
  class(metadata_list[[i]]) %in% c("dgCMatrix", "dgRMatrix")
}))
if(length(idx) > 0){
  graph_list <- metadata_list[idx]
  metadata_list <- metadata_list[-idx]
} else {
  graph_list <- numeric(0)
}
seurat_obj@misc <- metadata_list

if(length(graph_list) > 0){
  for(name_val in names(graph_list)){
    print(paste0("Putting in graph ", name_val))
    
    seurat_obj@graphs[[name_val]] <- graph_list[[name_val]]
    rownames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)
    colnames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)
    
  }
}

Seurat::DefaultAssay(seurat_obj) <- "RNA"

# for the cluster color specifically, add the names for convenience
names(seurat_obj@misc[["Cell.type.annotation_colors"]]) <- sort(unique(seurat_obj$Cell.type.annotation))

# clearing environment
ls_vec <- ls()
ls_vec <- ls_vec[ls_vec != "seurat_obj"]
rm(list = ls_vec)
gc(TRUE)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Conversion of Joshua's Larry analysis based on the scBaseEncoderFeat_Z_bs30_tau0.5.npy embedding (July 11, 2024)")


###############

# let's put in the usual Seurat pipeline
Seurat::VariableFeatures(seurat_obj) <- SeuratObject::Features(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
set.seed(10)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             features = Seurat::VariableFeatures(seurat_obj),
                             verbose = FALSE)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              dims = 1:30)

# fixing annoying names

seurat_obj$clone_id <- paste0("lin_", seurat_obj$clone_id)
seurat_obj$clone_id <- factor(seurat_obj$clone_id)
seurat_obj$Time.point <- paste0("t_", seurat_obj$Time.point)
seurat_obj$Time.point <- factor(seurat_obj$Time.point)
seurat_obj$Well <- paste0("w_", seurat_obj$Well)
seurat_obj$Well <- factor(seurat_obj$Well)

seurat_obj <- SeuratObject::RenameCells(seurat_obj,
                                        new.names =paste0("c:", Seurat::Cells(seurat_obj)))

save(seurat_obj,
     date_of_run, session_info, note,
     file = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup3/Writeup3_Larry_scBaseEncoderFeat_Z_bs30_tau05.RData")

######

Seurat::DimPlot(seurat_obj, 
                reduction = "python_X_umap", 
                group.by = "Cell.type.annotation",
                cols = seurat_obj@misc$Cell.type.annotation_colors)

Seurat::DimPlot(seurat_obj, 
                reduction = "umap", 
                group.by = "Cell.type.annotation",
                cols = seurat_obj@misc$Cell.type.annotation_colors)











