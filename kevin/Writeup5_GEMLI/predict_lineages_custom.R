predict_lineages_custom <- function(GEMLI_items, 
                                    bool_set_seed = TRUE,
                                    desired_cluster_size = c(2,3), 
                                    fast = TRUE,
                                    repetitions = 100, 
                                    sample_size = (2/3), 
                                    verbose = 0) # check
{
  data_matrix <- GEMLI_items[['gene_expression']]
  
  # check
  colsum_vec <- Matrix::colSums(data_matrix)
  rowsum_vec <- Matrix::rowSums(data_matrix)
  stopifnot(all(colsum_vec != 0), all(rowsum_vec != 0))
  
  if(verbose > 0) print("Compute potential markers")
  marker_genes <- GEMLI::potential_markers(data_matrix)
  results <- data.matrix(matrix(0, nrow=ncol(data_matrix), ncol=ncol(data_matrix)))
  rownames(results) <- colnames(data_matrix)
  colnames(results) <- colnames(data_matrix)
  
  for (i in 1:repetitions){
    if(verbose > 0) print(paste0("On iteration ", i, " in ", repetitions))
    
    if(bool_set_seed) set.seed(10*i)
    marker_genes_sample <- sample(intersect(marker_genes, rownames(data_matrix)), 
                                  round(length(intersect(marker_genes, rownames(data_matrix)))*sample_size,0))
    # the main workhorse
    cell_clusters <- quantify_clusters_iterative_custom(data_matrix, 
                                                        marker_genes_sample, 
                                                        N = 2, 
                                                        fast)
    
    cell_clusters_unique_name <- cell_clusters
    for (colname in 1:ncol(cell_clusters)){
      tmp <- paste0(colname,'_',cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname])
      cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname] <- tmp
    }
    clustersize_dict <- table(cell_clusters_unique_name)
    
    smallest_clusters <- names(clustersize_dict)[clustersize_dict %in% desired_cluster_size]
    best_prediction <- data.matrix(matrix(FALSE, nrow=ncol(data_matrix), ncol=ncol(data_matrix)))
    rownames(best_prediction) <- colnames(data_matrix)
    colnames(best_prediction) <- colnames(data_matrix)
    
    # only set "best_prediction" to clusters of the appropriate size
    for (cluster in smallest_clusters){
      cells_in_cluster <- rownames(best_prediction)[rowSums(cell_clusters_unique_name==cluster, na.rm=T)>0]
      best_prediction[cells_in_cluster,cells_in_cluster] <- TRUE
    }
    
    diag(best_prediction) <- FALSE
    results <- results + best_prediction
  }
  
  GEMLI_items[["prediction"]] <- results
  return(GEMLI_items)
}