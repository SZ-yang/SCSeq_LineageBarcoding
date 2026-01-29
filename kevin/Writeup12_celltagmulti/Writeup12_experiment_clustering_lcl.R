rm(list=ls())

library(Seurat)
library(mclust)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup12_joshua-celltagmulti/cell_tag_integration.RData")

seurat_obj <- subset(integrated, batch == "tag")

lcl_csv <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-lcl/test_proj_embe_w_clone_id.csv")
rownames(lcl_csv) <- lcl_csv[,"X"]
lcl_csv <- lcl_csv[,grep("embed_", colnames(lcl_csv))]
lcl_csv <- lcl_csv[Seurat::Cells(seurat_obj),]
colnames(lcl_csv) <- paste0("LCL_", 1:ncol(lcl_csv))
lcl_csv <- as.matrix(lcl_csv)

umap_csv <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-lcl/cellTag-cellTagMulti_UMAP_test.csv")
rownames(umap_csv) <- umap_csv$X
umap_csv <- umap_csv[,c("UMAP0", "UMAP1")]
colnames(umap_csv) <- paste0("lclumap_", 1:2)
umap_csv <- umap_csv[Seurat::Cells(seurat_obj),]
umap_csv <- as.matrix(umap_csv)

seurat_obj[["LCL"]] <- Seurat::CreateDimReducObject(lcl_csv,
                                                    assay = "integrated")
seurat_obj[["LCLUMAP"]] <- Seurat::CreateDimReducObject(umap_csv,
                                                        assay = "integrated")

DimPlot(seurat_obj,
        reduction = "LCLUMAP")

seurat_obj$clone_id <- factor(paste0("Lineage:", seurat_obj$clone_id))

########

seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    dims = 1:32, 
                                    reduction = "LCL")
seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                   resolution = 0.1,
                                   group.singletons = FALSE)

DimPlot(seurat_obj,
        reduction = "LCLUMAP",
        group.by = "seurat_clusters")

DimPlot(seurat_obj,
        reduction = "LCLUMAP",
        group.by = "cell_type")


tab_mat <- table(seurat_obj$clone_id,
                 seurat_obj$cell_type)
lineage_assignment <- sapply(1:nrow(tab_mat), function(i){
  if(tab_mat[i,"Fibroblast"] >= 2*(tab_mat[i,"iEP"]+1)){
    "commit_fibroblast"
  } else if(tab_mat[i,"iEP"] >= 2*(tab_mat[i,"Fibroblast"]+1)){
    "commit_iep"
  } else {
    "na"
  }
})
names(lineage_assignment) <- rownames(tab_mat)
vec <- lineage_assignment[seurat_obj$clone_id]
names(vec) <- Seurat::Cells(seurat_obj)
seurat_obj$fate <- vec


DimPlot(seurat_obj,
        reduction = "LCLUMAP",
        group.by = "fate")


################

md <- seurat_obj@meta.data

# ---- helper: sample up to mmax UNIQUE unordered pairs and compute agreement ----
pair_agreement <- function(clust, mmax = 500) {
  n <- length(clust)
  if (n < 2) return(NA_real_)
  
  M <- n * (n - 1) / 2
  m <- min(mmax, M)
  
  # If total unique pairs is small (<= 500), just enumerate all pairs safely
  if (M <= mmax) {
    ij <- utils::combn(n, 2)
    return(mean(clust[ij[1, ]] == clust[ij[2, ]]))
  }
  
  # Otherwise: sample m unique unordered pairs without enumerating all M
  keys_seen <- integer(0)
  a_keep <- integer(0)
  b_keep <- integer(0)
  
  while (length(keys_seen) < m) {
    p <- max(2000L, (m - length(keys_seen)) * 10L)  # batch size
    i <- sample.int(n, p, replace = TRUE)
    j <- sample.int(n, p, replace = TRUE)
    ok <- i != j
    i <- i[ok]; j <- j[ok]
    
    a <- pmin(i, j)
    b <- pmax(i, j)
    key <- a + (b - 1L) * n  # unique encoding for unordered (a<b)
    
    # keep new keys only
    new <- !key %in% keys_seen
    if (any(new)) {
      key <- key[new]; a <- a[new]; b <- b[new]
      keys_seen <- c(keys_seen, key)
      a_keep   <- c(a_keep, a)
      b_keep   <- c(b_keep, b)
    }
  }
  
  a_keep <- a_keep[seq_len(m)]
  b_keep <- b_keep[seq_len(m)]
  mean(clust[a_keep] == clust[b_keep])
}

set.seed(1)  # for reproducibility of subsampling

# ---- compute per-clone and overall average ----
clone_ids <- md$clone_id
clust_all <- md$seurat_clusters

split_idx <- split(seq_len(nrow(md)), clone_ids)

per_clone <- data.frame(
  clone_id = names(split_idx),
  n_cells  = vapply(split_idx, length, integer(1)),
  same_cluster_rate = vapply(split_idx, function(ii) {
    pair_agreement(clust_all[ii], mmax = 500)
  }, numeric(1))
)
plot(per_clone$n_cells, per_clone$same_cluster_rate,
     main = paste0("LCL: ", length(unique(md$seurat_clusters))))

overall_mean <- mean(per_clone$same_cluster_rate, na.rm = TRUE)
overall_mean
