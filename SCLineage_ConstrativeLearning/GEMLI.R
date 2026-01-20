############################################
## h5ad -> Seurat (store normalized X in RNA@data)
## - NO scVI embedding
## - NO NormalizeData() (avoid double-normalization)
## - Uses layers["counts"] as counts if available; otherwise makes 0-count placeholder
############################################

rm(list = ls())

############################################
## 0. Setup Python (reticulate + conda env)
############################################
library(reticulate)
use_condaenv("r-h5ad", required = TRUE)

############################################
## 1. Read h5ad via anndata (NO basilisk)
############################################
ad <- import("anndata", delay_load = FALSE)
adata <- ad$read_h5ad(
  "/Users/apple/Desktop/KB/data/Cell_tag-Cell_tag_multi_integrated/Seurat_method/cellTag_test_multi_Seurat_clone_id.h5ad"
)

############################################
## 2. Extract clone IDs + names
############################################
obs <- py_to_r(adata$obs)
clone_id <- obs$clone_id

cell_names <- as.character(py_to_r(adata$obs_names$to_list()))
gene_names <- as.character(py_to_r(adata$var_names$to_list()))
names(clone_id) <- cell_names

############################################
## 3. Extract normalized expression from adata.X (cells x genes)
############################################
X_norm_cxg <- py_to_r(adata$X)
X_norm_cxg <- as.matrix(X_norm_cxg)  # may be large

rownames(X_norm_cxg) <- cell_names
colnames(X_norm_cxg) <- gene_names

# Seurat expects genes x cells
X_norm_gxc <- t(X_norm_cxg)

############################################
## 4. Optionally extract raw counts from layers["counts"] (for Seurat counts layer)
############################################
has_layers <- !is.null(adata$layers)

layer_keys <- if (has_layers) {
  tryCatch(py_to_r(adata$layers$keys()), error = function(e) NULL)
} else {
  NULL
}

has_counts_layer <- !is.null(layer_keys) && ("counts" %in% as.character(layer_keys))

if (has_counts_layer) {
  counts_cxg <- py_to_r(adata$layers[["counts"]])
  counts_cxg <- as.matrix(counts_cxg)  # may be large
  rownames(counts_cxg) <- cell_names
  colnames(counts_cxg) <- gene_names
  counts_gxc <- t(counts_cxg)
} else {
  # Placeholder counts required by CreateSeuratObject; not used downstream
  counts_gxc <- matrix(
    0,
    nrow = length(gene_names),
    ncol = length(cell_names),
    dimnames = list(gene_names, cell_names)
  )
}

############################################
## 5. Create Seurat object and set normalized data layer directly
############################################
library(Seurat)
library(Matrix)

counts_gxc <- Matrix(counts_gxc, sparse = TRUE)
X_norm_gxc <- Matrix(X_norm_gxc, sparse = TRUE)

seurat_obj <- CreateSeuratObject(counts = counts_gxc)
seurat_obj$clone_id <- clone_id[Cells(seurat_obj)]

# IMPORTANT: put normalized expression into RNA "data" layer (NO NormalizeData())
SeuratObject::LayerData(seurat_obj, assay = "RNA", layer = "data") <- X_norm_gxc

############################################
## 6. Sanity checks
############################################
cat("RNA dim (genes x cells): ",
    paste(dim(seurat_obj[["RNA"]]), collapse = " x "), "\n")

cat("Layers in RNA: ",
    paste(SeuratObject::Layers(seurat_obj[["RNA"]]), collapse = ", "), "\n")

data_layer <- SeuratObject::LayerData(seurat_obj, assay = "RNA", layer = "data")
cat("Nonzero entries in data layer: ", Matrix::nnzero(data_layer), "\n")

# Optional: confirm a small slice prints (won't error)
print(data_layer[1:5, 1:5])

############################################
## 7. Save
############################################
saveRDS(seurat_obj, "cellTagMulti_normX_seurat.rds")

cat("Saved: cellTagMulti_normX_seurat.rds\n")
