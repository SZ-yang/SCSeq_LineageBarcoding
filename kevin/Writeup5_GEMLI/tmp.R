# zoom in
set.seed(10)
bool_set_seed = TRUE
desired_cluster_size = c(75, 120)
fast = TRUE
repetitions = 100
sample_size = (2/3)
verbose = 1

data_matrix <- GEMLI_items[['gene_expression']]
if(verbose > 0) print("Compute potential markers")
marker_genes <- GEMLI::potential_markers(data_matrix)
results <- data.matrix(matrix(0, nrow=ncol(data_matrix), ncol=ncol(data_matrix)))
rownames(results) <- colnames(data_matrix)
colnames(results) <- colnames(data_matrix)

i <- 6
if(verbose > 0) print(paste0("On iteration ", i, " in ", repetitions))
if(bool_set_seed) set.seed(10*i)
marker_genes_sample <- sample(intersect(marker_genes, rownames(data_matrix)), 
                              round(length(intersect(marker_genes, rownames(data_matrix)))*sample_size,0))
# the main workhorse
cell_clusters <- quantify_clusters_iterative_custom(data_matrix, 
                                                    marker_genes_sample, 
                                                    N = 2, 
                                                    fast)

##########

# zoom in
marker_genes <- marker_genes_sample
N <- 2
iterate <- TRUE
i <- 2
genes <- intersect(marker_genes, rownames(data_matrix)[rowMeans(data_matrix)>0])
data_matrix <- data_matrix[genes,]

# this probably can be changed to handle sparse matrices...
corr_expr_raw <- GEMLI:::calculate_correlations(t(data_matrix), fast=FALSE)
corr_expr_raw[is.na(corr_expr_raw)] <- 0
corr_expr <- (1 - corr_expr_raw)/2

cell_clusters <- data.matrix(matrix(0, nrow=ncol(data_matrix), ncol=1))
rownames(cell_clusters) <- colnames(data_matrix)
cell_clusters[,1] <- rep(1, ncol(data_matrix))

while (iterate) {
  cell_clusters <- cbind(cell_clusters, rep(0,nrow(cell_clusters)))
  
  for (cluster in setdiff(unique(cell_clusters[,(i-1)]),0)) {
    cells_in_cluster <- rownames(cell_clusters)[cell_clusters[,(i-1)]==cluster]
    if (length(cells_in_cluster) >= 4) { # this line ends the sub clustering # min of desired cluster size
      correlation <- mean((corr_expr_raw[cells_in_cluster,cells_in_cluster])[lower.tri(corr_expr_raw[cells_in_cluster,cells_in_cluster], diag=F)])
      corr_expr_subset <- corr_expr[cells_in_cluster,cells_in_cluster]
      
      # main function
      clustering <- stats::cutree(stats::hclust(stats::as.dist(corr_expr_subset), method = "ward.D2", ), k=N)
      
      cell_clusters[names(clustering),i] <- as.vector(clustering) + max(c(0, cell_clusters[,i]), na.rm=TRUE)
    }
    else {cell_clusters[cells_in_cluster,i] = 0}
  }
  
  if (sum(cell_clusters[,i], na.rm=TRUE)==0) {iterate = FALSE}
  i <- i+1
}

# verify this doesn't work
clustering <- stats::cutree(stats::hclust(stats::as.dist(corr_expr_subset), method = "ward.D2"), k=N)

table(is.na(corr_expr_subset))
table(is.nan(corr_expr_subset))
table(is.infinite(corr_expr_subset))

tmp <- which(is.na(corr_expr_subset), arr.ind = TRUE)

# index 1694
# LK_d6_2_2:CCGACTTTGCTTACCT

corr_expr_raw[1694,1:10]
zz <- t(data_matrix)[1694,]
