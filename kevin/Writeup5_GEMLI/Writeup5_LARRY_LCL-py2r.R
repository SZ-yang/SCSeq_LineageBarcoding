rm(list=ls()); gc(TRUE)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v5")

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"

# from https://github.com/theislab/zellkonverter/issues/38
adata <- zellkonverter::readH5AD(
  file = paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.h5ad")
)

names(adata@assays)

tmp <- Seurat::as.Seurat(adata, counts = "raw_counts", data = "X")
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

# check
mat <- SeuratObject::LayerData(seurat_obj,
                               layer = "counts",
                               assay = "RNA")
mat[1:10,1:10]
mat <- SeuratObject::LayerData(seurat_obj,
                               layer = "data",
                               assay = "RNA")
mat[1:10,1:10]

# loading the LCL embedding
lcl_embedding <- read.csv("/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/scBaseEncoderFeat_Z_bs260_tau0.5.csv")
colnames(lcl_embedding) <- paste0("LCL_", 1:ncol(lcl_embedding))
rownames(lcl_embedding) <- Seurat::Cells(seurat_obj)
lcl_embedding <- as.matrix(lcl_embedding)
seurat_obj[["LCL"]] <- Seurat::CreateDimReducObject(lcl_embedding)

# clearing environment
ls_vec <- ls()
ls_vec <- ls_vec[ls_vec != "seurat_obj"]
rm(list = ls_vec)
gc(TRUE)

# some minor adjustments
## renaming it SPRING
seurat_obj[["SPRING"]] <- seurat_obj[["python_X_emb"]]
seurat_obj[["python_X_emb"]] <- NULL
seurat_obj[["python_X_clone"]] <- NULL

# compute the LCL UMAP
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "LCL",
                              reduction.name = "LCLumap",
                              dims = 1:64)

date_of_run <- Sys.time()
session_info <- devtools::session_info()

seurat_obj$clone_id <- paste0("Lineage_", seurat_obj$clone_id)

save(seurat_obj,
     date_of_run, session_info, 
     file = paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))

#############

plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup5/"

pdf(paste0(plot_folder, "Writeup5_Larry_41093_2000_norm_log_cleaned-covariates.pdf"),
    onefile = T, width = 5, height = 5)

plot1 <- Seurat::DimPlot(seurat_obj, 
                         reduction = "LCLumap",
                         group.by = "clone_id",
                         raster = TRUE)
plot1 <- plot1 + Seurat::NoLegend()
print(plot1)

tab_vec <- table(seurat_obj$clone_id)
largest_lineage <- names(tab_vec)[which.max(tab_vec)]
idx <- which(seurat_obj$clone_id == largest_lineage)

cell_names <- Seurat::Cells(seurat_obj)[idx]
plot1 <- Seurat::DimPlot(object = seurat_obj, 
                         reduction = "LCLumap",
                         cells.highlight = cell_names, 
                         cols.highlight = "red", 
                         cols = "gray", 
                         order = TRUE,
                         raster = TRUE)
plot1 <- plot1 + Seurat::NoLegend()
print(plot1)


plot1 <- Seurat::DimPlot(seurat_obj2, 
                         reduction = "LCLumap",
                         group.by = "clone_id",
                         raster = TRUE)
plot1 <- plot1 + Seurat::NoLegend()
print(plot1)

plot1 <- Seurat::DimPlot(seurat_obj, 
                         reduction = "LCLumap",
                         group.by = "Time.point",
                         raster = TRUE)
print(plot1)

plot1 <- Seurat::DimPlot(seurat_obj, 
                         reduction = "LCLumap",
                         group.by = "state_info",
                         raster = TRUE)
print(plot1)

plot1 <- Seurat::DimPlot(seurat_obj, 
                         reduction = "SPRING",
                         group.by = "Time.point",
                         raster = TRUE)
print(plot1)

plot1 <- Seurat::DimPlot(seurat_obj, 
                         reduction = "SPRING",
                         group.by = "state_info",
                         raster = TRUE)
print(plot1)

dev.off()

###########################

# plotting the top 20 lineages

pdf(paste0(plot_folder, "Writeup5_Larry_41093_2000_norm_log_cleaned_top-lineages.pdf"),
    onefile = T, width = 5, height = 5)

tab_vec <- table(seurat_obj$clone_id)
order_idx <- order(tab_vec, decreasing = TRUE)[1:20]
lineage_vec <- names(tab_vec)[order_idx]

for(i in 1:20){
  lineage_name <- lineage_vec[i]
  idx <- which(seurat_obj$clone_id == lineage_name)
  
  cell_names <- Seurat::Cells(seurat_obj)[idx]
  plot1 <- Seurat::DimPlot(object = seurat_obj, 
                           reduction = "LCLumap",
                           cells.highlight = cell_names, 
                           cols.highlight = "red", 
                           cols = "gray", 
                           order = TRUE)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0(lineage_name, " with ", length(idx), " cells"))
  
  print(plot1)
}

dev.off()

########


