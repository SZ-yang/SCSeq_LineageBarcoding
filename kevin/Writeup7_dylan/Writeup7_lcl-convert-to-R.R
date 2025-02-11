rm(list=ls())
library(SeuratDisk)
library(Seurat)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"

# load(paste0(out_folder, "Writeup7_shaffer_preprocessed.RData"))

# Convert the H5AD file to an H5Seurat file
Convert(paste0(out_folder, "adata_with_lcl.h5ad"), 
        dest = "h5seurat", 
        overwrite = TRUE)

# Load the converted file as a Seurat object
seurat_obj <- LoadH5Seurat(paste0(out_folder, "adata_with_lcl.h5seurat"))

metadata <- read.csv(paste0(out_folder, "adata_obs_backup.csv"))
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]
seurat_obj@meta.data <- metadata

###########

# put in the colors

color_vec <- c("cis" = rgb(202,111,55, maxColorValue = 255),
               "cistocis" = rgb(248,210,152, maxColorValue = 255),
               "cistococl2" = rgb(240,148,71, maxColorValue = 255),
               "cistodabtram" = rgb(160,61,38, maxColorValue = 255),
               "cocl2" = rgb(69,132,69, maxColorValue = 255),
               "cocl2tocis" = rgb(131,202,163, maxColorValue = 255),
               "cocl2tococl2" = rgb(126,191,90, maxColorValue = 255),
               "cocl2todabtram" = rgb(35,63,58, maxColorValue = 255),
               "dabtram" = rgb(68,49,147, maxColorValue = 255),
               "dabtramtocis" = rgb(147,137,193, maxColorValue = 255),
               "dabtramtococl2" = rgb(145,54,147, maxColorValue = 255),
               "dabtramtodabtram" = rgb(68,32,85, maxColorValue = 255))

seurat_obj@misc[["OG_condition_colors"]] <- color_vec

#########

# label the cells in the main LCL clusters

seurat_obj$OG_condition <- factor(seurat_obj$OG_condition)
seurat_obj$Lineage <- factor(seurat_obj$Lineage)

# first round of clustering
set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    reduction = "lcl",
                                    dims = 1:64)
seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.01)

tab_mat <- table(droplevels(seurat_obj$Lineage), 
                 droplevels(seurat_obj$RNA_snn_res.0.01))
cluster_sizes <- colSums(tab_mat)
cluster_names <- colnames(tab_mat)[which(cluster_sizes > 200)]
keep_vec <- (seurat_obj$RNA_snn_res.0.01 %in% cluster_names)
seurat_obj$keep <- keep_vec
seurat_obj2 <- subset(seurat_obj, keep == TRUE)

# dbscan
embedding <- seurat_obj2[["lcl"]]@cell.embeddings
dbscan_res <- dbscan::dbscan(
  x = embedding,
  eps = 2,
  minPts = 50
)
dbscan_cluster <- dbscan_res$cluster
dbscan_cluster[dbscan_cluster == 0] <- NA
names(dbscan_cluster) <- Seurat::Cells(seurat_obj2)
seurat_obj2$dbscan <- dbscan_cluster
seurat_obj2$keep <- !is.na(seurat_obj2$dbscan)
seurat_obj2 <- subset(seurat_obj2, keep == TRUE)

# second round of clustering
seurat_obj2 <- Seurat::FindClusters(seurat_obj2, 
                                   resolution = 0.01)
tab_vec <- table(seurat_obj2$RNA_snn_res.0.01)
cluster_names <- names(tab_vec)[which(tab_vec > 200)]
keep_vec <- (seurat_obj2$RNA_snn_res.0.01 %in% cluster_names)
seurat_obj2$keep <- keep_vec
seurat_obj2 <- subset(seurat_obj2, keep == TRUE)

full_cells <- Seurat::Cells(seurat_obj)
selected_cells <- Seurat::Cells(seurat_obj2)
keep_vec <- rep(FALSE, length(full_cells))
keep_vec[which(full_cells %in% selected_cells)] <- TRUE

seurat_obj$in_lcl_cluster <- keep_vec

seurat_obj$Lineage <- droplevels(factor(paste0("Lineage:", seurat_obj$Lineage)))
seurat_obj$RNA_snn_res.0.01 <- factor(paste0("Cluster:", seurat_obj$RNA_snn_res.0.01))

###############

seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = 'lcl', 
                              dims = 1:64, 
                              assay = 'RNA', 
                              reduction.name = 'lcl.umap', 
                              reduction.key = 'lclUMAP_')

save(seurat_obj,
     file = paste0(out_folder, "adata_with_lcl.RData"))


