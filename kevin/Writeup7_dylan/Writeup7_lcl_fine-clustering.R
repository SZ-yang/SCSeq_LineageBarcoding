rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dbscan)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

# binarize the lcl embedding
embedding <- seurat_obj[["lcl"]]@cell.embeddings
embedding[embedding > 0] <- 1
character_str <- sapply(1:nrow(embedding), function(i){
  tmp <- which(embedding[i,] != 0)
  paste0(tmp, collapse = "-")
})
length(unique(character_str))
character_str <- factor(character_str)

seurat_obj$cluster <- character_str

plot1 <- scCustomize::DimPlot_scCustom(seurat_obj,
                                       group.by = "cluster",
                                       reduction = "lcl.umap") + Seurat::NoLegend()

##################

embedding <- seurat_obj[["lcl"]]@cell.embeddings

set.seed(10)
knn_res <- RANN::nn2(data = embedding, 
                     query = embedding,
                     k = 10)
nn_mat <- knn_res$nn.idx
nn_mat <- nn_mat[,-1]
i_vec <- as.numeric(t(nn_mat))
j_vec <- rep(1:nrow(nn_mat), each = ncol(nn_mat))
sparse_mat <- Matrix::sparseMatrix(i = i_vec, 
                                   j = j_vec, 
                                   x = rep(1, length(i_vec)),
                                   dims = c(nrow(embedding), nrow(embedding)))
sparse_mat <- sparse_mat * Matrix::t(sparse_mat)


# Function to find connected components in a large sparse adjacency matrix
find_connected_components_rows <- function(sparse_mat_r) {
  n <- nrow(sparse_mat_r)
  visited <- rep(FALSE, n)
  component_labels <- integer(n)  # or rep(0, n)
  component_id <- 0
  
  for (row_i in seq_len(n)) {
    if(row_i %% floor(n/10) == 0) cat('*')
    
    if (!visited[row_i]) {
      component_id <- component_id + 1
      # Start BFS at row_i
      queue <- row_i
      
      visited[row_i] <- TRUE
      component_labels[row_i] <- component_id
      
      while (length(queue) > 0) {
        current <- queue[1]
        queue   <- queue[-1]
        
        # In row-compressed format, neighbors of 'current' row
        start_idx <- sparse_mat_r@p[current]
        end_idx   <- sparse_mat_r@p[current + 1] - 1
        
        # If start_idx <= end_idx, it means we have neighbors
        if (start_idx <= end_idx) {
          # @j holds the column indices of nonzeros in row `current`
          neighbors <- sparse_mat_r@j[(start_idx + 1):(end_idx + 1)] + 1
          for (nb in neighbors) {
            if (!visited[nb]) {
              visited[nb] <- TRUE
              component_labels[nb] <- component_id
              queue <- c(queue, nb)
            }
          }
        }
      }
    }
  }
  
  return(component_labels)
}


# Example usage
set.seed(42)
sparse_mat <- as(sparse_mat, "RsparseMatrix")
components <- find_connected_components_rows(sparse_mat)
components_safe <- components

# Check number of unique components
length(unique(components))

seurat_obj$cluster <- factor(components)

plot1 <- scCustomize::DimPlot_scCustom(seurat_obj,
                                       group.by = "cluster",
                                       reduction = "lcl.umap") + Seurat::NoLegend()
plot1

#########

# merge all the lineages
iter <- 1
components_updated <- components
while(any(components_updated != components) | iter == 1){
  print(iter)
  components <- components_updated
  
  lineage_vec <- seurat_obj$Lineage
  lineage_names <- unique(lineage_vec)
  for(lineage in lineage_names){
    idx <- which(lineage_vec == lineage)
    min_val <- min(components_updated[idx])
    components_updated[idx] <- min_val
  }
  
  iter <- iter+1
}
components <- components_updated

seurat_obj$cluster <- factor(components)
plot1 <- scCustomize::DimPlot_scCustom(seurat_obj,
                                       group.by = "cluster",
                                       reduction = "lcl.umap") + Seurat::NoLegend()
plot1