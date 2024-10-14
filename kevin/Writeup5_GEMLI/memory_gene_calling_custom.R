memory_gene_calling_custom <- function(GEMLI_items, 
                                       ground_truth = TRUE,
                                       num_trials = 20,
                                       use_median = TRUE, 
                                       valid_lineage_sizes = c(50:200), 
                                       verbose = 0){
  data_matrix <- GEMLI_items[['gene_expression']]
  
  if (ground_truth) {
    lineage_dict <- GEMLI_items[['barcodes']]
  } else {
    if (length(GEMLI_items[['predicted_lineages']])>0){
      lineage_dict <- GEMLI_items[['predicted_lineages']]
    } else {
      lineage_dict <- GEMLI_items[['predicted_lineage_table']]$clone.ID
      names(lineage_dict) <- GEMLI_items[['predicted_lineage_table']]$cell.ID
    }
  }
  
  if (verbose > 0) print("Compute value on observed data")
  lineage_center_variation <- .markers_by_cvsq_of_lineage_means(data_matrix, 
                                                                lineage_dict = lineage_dict, 
                                                                valid_lineage_sizes = valid_lineage_sizes, 
                                                                use_median = use_median,
                                                                verbose = verbose - 1)
  
  data_matrix_control <- matrix(NA, ncol=num_trials, nrow=nrow(data_matrix))
  rownames(data_matrix_control) <- rownames(data_matrix)
  
  if (verbose > 0) print("Compute null distribution")
  for (i in c(1:num_trials)) {
    if(verbose == 1 && num_trials > 10 && i %% floor(num_trials/10) == 0) cat('*')
    if(verbose > 2) print(paste0("On trial", i))
    
    lineage_dict_sampled <- lineage_dict
    names(lineage_dict_sampled) <- sample(names(lineage_dict))
    tmp <- .markers_by_cvsq_of_lineage_means(as.matrix(data_matrix), 
                                             lineage_dict = lineage_dict_sampled, 
                                             valid_lineage_sizes = valid_lineage_sizes, 
                                             use_median = use_median,
                                             verbose = verbose - 1)
    data_matrix_control[names(tmp),i] <- tmp
  }
  
  if (verbose > 0) print("Compute p-value")
  df_1 <- data_matrix_control[intersect(rownames(data_matrix_control), names(lineage_center_variation)),,drop = FALSE]
  df_2 <- lineage_center_variation[intersect(rownames(data_matrix_control), names(lineage_center_variation)),drop = FALSE]
  markers_pvalue <- Matrix::rowSums(df_1 > df_2, na.rm=TRUE)/num_trials
  
  shared_genes <- intersect(names(lineage_center_variation), names(markers_pvalue))
  marker_table <- data.frame(cbind(lineage_center_variation[shared_genes], markers_pvalue[shared_genes]))
  rownames(marker_table) <- shared_genes
  colnames(marker_table) <- c("var","p")
  marker_table <- marker_table[with(marker_table, order(p, -var)),]
  
  GEMLI_items[["memory_genes"]] <- marker_table
  
  return(GEMLI_items)
}

##################

.cv_sq <- function(mat) {
  sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
  mean <- apply(mat, 1, mean, na.rm = TRUE)
  noise <- (sd/mean)**2
  return(noise)
}

.markers_by_cvsq_of_lineage_means <- function(data_matrix, 
                                              lineage_dict, 
                                              valid_lineage_sizes = c(50:200), 
                                              use_median = TRUE,
                                              verbose = 0){
  lineage_dict_filt <- lineage_dict[intersect(colnames(data_matrix), names(lineage_dict))]
  valid_lineage_dict <- lineage_dict_filt[as.character(lineage_dict_filt) %in% names(table(lineage_dict_filt))[table(lineage_dict_filt) %in% valid_lineage_sizes]]
  lineage_center <- matrix(NA, ncol=length(unique(valid_lineage_dict)), nrow=nrow(data_matrix))
  colnames(lineage_center) <- unique(valid_lineage_dict)
  rownames(lineage_center) <- rownames(data_matrix)
  
  if (verbose > 0) print("Compute lineage centers")
  for (lineage in as.character(unique(valid_lineage_dict))) {
    tmp_mat <- data_matrix[,names(valid_lineage_dict)[valid_lineage_dict==lineage],drop = FALSE]
    if (use_median){
      lineage_center[,lineage] <- apply(tmp_mat, 1, stats::median, na.rm = TRUE)
    }
    else {
      lineage_center[,lineage] <- Matrix::rowMeans(tmp_mat)
    }
  }
  
  x <- Matrix::rowMeans(lineage_center)
  y <- .cv_sq(lineage_center)
  x <- log2(x)
  y <- log2(y)
  filter <- (is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y) | (x==0))
  x <- x[!filter]
  y <- y[!filter]
  if(verbose > 0) print("Compute loess")
  loess_means <- stats::loess(y ~ x, span = 0.75, control = stats::loess.control(surface="direct"))
  filter <- names(which(!filter))
  lineage_center_variation <- loess_means$residuals
  
  return(sort(lineage_center_variation[filter], decreasing = TRUE))
}