
plot_anova <- function(seurat_object,
                       cell_imputed_score,
                       assigned_lineage_variable,
                       time_celltype_variable,
                       day_later,
                       bool_mark_mean = TRUE,
                       bool_mark_median = TRUE,
                       min_lineage_size = 2,
                       num_lineages = 20,
                       ylab = "",
                       ylim = NA){
  
  # grab the vector of which celltype-time each cell is
  assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
  names(assigned_lineage) <- Seurat::Cells(seurat_object)
  
  time_celltype <- seurat_object@meta.data[,time_celltype_variable]
  names(time_celltype) <- Seurat::Cells(seurat_object)
  stopifnot(day_later %in% time_celltype)
  
  # determine which lineages qualify to be in the plot
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  tab_mat <- table(assigned_lineage, time_celltype)
  lineage_future_size <- tab_mat[, day_later]
  names(lineage_future_size) <- rownames(tab_mat)
  
  .plot_anova_helper(seurat_object = seurat_object,
                     cell_imputed_score = cell_imputed_score,
                     assigned_lineage_variable = assigned_lineage_variable,
                     lineage_future_size = lineage_future_size,
                     bool_mark_mean = bool_mark_mean,
                     bool_mark_median = bool_mark_median,
                     min_lineage_size = min_lineage_size,
                     num_lineages = num_lineages,
                     ylab = ylab,
                     ylim = ylim)
}

#########################

.plot_anova_helper <- function(seurat_object,
                               cell_imputed_score,
                               assigned_lineage_variable,
                               lineage_future_size,
                               bool_mark_mean = TRUE,
                               bool_mark_median = TRUE,
                               min_lineage_size = 2,
                               num_lineages = 20,
                               ylab = "",
                               ylim = NA){
  stopifnot(length(names(cell_imputed_score)) == length(cell_imputed_score))
  
  if(any(is.na(cell_imputed_score))){
    cell_imputed_score <- cell_imputed_score[!is.na(cell_imputed_score)]
  }
  
  # grab the vector of which celltype-time each cell is
  assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
  names(assigned_lineage) <- Seurat::Cells(seurat_object)
  
  # determine which lineages qualify to be in the plot
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  tab_vec <- table(assigned_lineage)
  tab_vec <- tab_vec[tab_vec >= min_lineage_size] # current size needs to be big enough
  lineage_names <- names(lineage_future_size)[order(lineage_future_size, decreasing = TRUE)[1:num_lineages]]
  idx <- which(lineage_vec %in% lineage_names)
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[idx],
                   imputed_count = cell_imputed_score[idx])
  df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
  anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    imputed_count = cell_imputed_score)
  df <- rbind(df, df2)
  
  # compute percentage
  lineage_effect <- .anova_percentage(
    df = df_tmp,
    lineage_variable = "lineage",
    value_variable = "imputed_count"
  )
  
  col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
  names(col_vec) <- c(lineage_names, "All")
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=imputed_count))
  plot1 <- plot1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=lineage))
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = col_vec) 
  plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                        position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::geom_boxplot(width=0.05)
  plot1 <- plot1 + ggplot2::scale_x_discrete(limits = c(lineage_names, "All"),
                                             guide = ggplot2::guide_axis(angle = 45))
  plot1 <- plot1 + ggplot2::ylab(ylab)
  
  if(!all(is.na(ylim))){
    plot1 <- plot1 + ggplot2::ylim(ylim[1], ylim[2])
  }
  
  if(bool_mark_mean) 
    plot1 <- plot1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  
  if(bool_mark_median) 
    plot1 <- plot1 + ggplot2::stat_summary(fun=max, geom="point", shape=10, size=5, color="blue")
  
  plot1 <- plot1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2), 
                                           ", Lineage effect = ", lineage_effect, "%"))
  
  plot1
}

##################

.anova_percentage <- function(df,
                              lineage_variable,
                              value_variable){
  stopifnot(is.factor(df[,lineage_variable]))
  imputed_count <- df[,value_variable]
  lineage <- df[,lineage_variable]
  
  total_std <- sum((imputed_count - mean(imputed_count))^2)
  
  within_lineage_std <- sum(sapply(levels(lineage), function(lineage_name){
    idx <- which(lineage == lineage_name)
    sum((imputed_count[idx] - mean(imputed_count[idx]))^2)
  }))
  
  across_lineage_std <- sum(sapply(levels(lineage), function(lineage_name){
    idx <- which(lineage == lineage_name)
    mean_val <- mean(imputed_count[idx])
    length(idx) * (mean_val - mean(imputed_count))^2 
  }))
  
  lineage_effect <- round(across_lineage_std/total_std*100,1)
  lineage_effect
}