rm(list=ls())

library(Seurat)

data_folder <- "~/kzlinlab/data/shaffer_clonal-treatment/"
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/"
fig_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup7/"

load(paste0(out_folder, "Writeup7_shaffer_preprocessed.RData"))

seurat_object <- all_data
assigned_lineage_variable <- "Lineage"
num_lineages <- 50
ylim <- NA
bool_mark_mean <- TRUE
bool_mark_median <- TRUE
min_lineage_size <- 2

scaled_data <- SeuratObject::LayerData(all_data,
                                       assay = "RNA",
                                       layer = "scale.data",
                                       features = Seurat::VariableFeatures(all_data))

# grab the vector of which celltype-time each cell is
assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
names(assigned_lineage) <- Seurat::Cells(seurat_object)

# determine which lineages qualify to be in the plot
lineage_vec <- assigned_lineage
tab_vec <- table(assigned_lineage)
lineage_names <- names(tab_vec)[order(tab_vec, decreasing = TRUE)[1:num_lineages]]
cell_idx <- which(lineage_vec %in% lineage_names)

genes_of_interest <- c("FMR1NB", "LRRC52", "SLRC22A3", "NPY2R", "S100B", "IGHG3")
for(gene in genes_of_interest){
  
  if(!gene %in% rownames(scaled_data)) next()
  gene_vec <- scaled_data[gene,]
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[cell_idx],
                   expression = gene_vec[cell_idx])
  df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
  anova_res <- stats::oneway.test(expression ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    expression = gene_vec)
  df <- rbind(df, df2)
  
  # compute percentage
  lineage_effect <- .anova_percentage(
    df = df_tmp,
    lineage_variable = "lineage",
    value_variable = "expression"
  )
  
  col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
  names(col_vec) <- c(lineage_names, "All")
  
  ylab <- paste0("Scaled expression of ", gene)
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=expression))
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
  
  ggplot2::ggsave(plot1, 
                 file = paste0(fig_folder, "violin_high-lineage_", gene, ".png"),
                 height = 4, width = 10)
}


genes_of_interest <- c("FMR1NB", "LRRC52", "SLRC22A3", "NPY2R", "S100B", "IGHG3")
for(gene in genes_of_interest){
  
  if(!gene %in% rownames(scaled_data)) next()
  gene_vec <- scaled_data[gene,]
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[cell_idx],
                   expression = gene_vec[cell_idx])
  df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
  anova_res <- stats::oneway.test(expression ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    expression = gene_vec)
  df <- rbind(df, df2)
  
  # compute percentage
  lineage_effect <- .anova_percentage(
    df = df_tmp,
    lineage_variable = "lineage",
    value_variable = "expression"
  )
  
  col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
  names(col_vec) <- c(lineage_names, "All")
  
  ylab <- paste0("Scaled expression of ", gene)
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=expression))
  plot1 <- plot1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=lineage))
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = col_vec) 
  plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                        position=ggplot2::position_jitter(width = 0.2, height = 0), alpha = 0.3, size = 0.5)
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
  
  ggplot2::ggsave(plot1, 
                  file = paste0(fig_folder, "violin_high-lineage_", gene, ".png"),
                  height = 4, width = 10)
}

#########


genes_of_interest <- c("FMR1NB", "LRRC52", "SLRC22A3", "NPY2R", "S100B", "IGHG3")
for(gene in genes_of_interest){
  
  if(!gene %in% rownames(scaled_data)) next()
  gene_vec <- scaled_data[gene,]
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[cell_idx],
                   expression = gene_vec[cell_idx])
  df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
  anova_res <- stats::oneway.test(expression ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    expression = gene_vec)
  df <- rbind(df, df2)
  
  # compute percentage
  lineage_effect <- .anova_percentage(
    df = df_tmp,
    lineage_variable = "lineage",
    value_variable = "expression"
  )
  
  col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
  names(col_vec) <- c(lineage_names, "All")
  
  ylab <- paste0("Scaled expression of ", gene)
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=expression))
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
  
  ggplot2::ggsave(plot1, 
                  file = paste0(fig_folder, "violin_high-lineage_", gene, ".png"),
                  height = 4, width = 10)
}


genes_of_interest <- c("FN1", "PMEL", "ANXA2", "CALD1", "ALCAM", "MYOF")
for(gene in genes_of_interest){
  
  if(!gene %in% rownames(scaled_data)) next()
  gene_vec <- scaled_data[gene,]
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[cell_idx],
                   expression = gene_vec[cell_idx])
  df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
  anova_res <- stats::oneway.test(expression ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    expression = gene_vec)
  df <- rbind(df, df2)
  
  # compute percentage
  lineage_effect <- .anova_percentage(
    df = df_tmp,
    lineage_variable = "lineage",
    value_variable = "expression"
  )
  
  col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
  names(col_vec) <- c(lineage_names, "All")
  
  ylab <- paste0("Scaled expression of ", gene)
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=expression))
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
  
  ggplot2::ggsave(plot1, 
                  file = paste0(fig_folder, "violin_low-lineage_", gene, ".png"),
                  height = 4, width = 10)
}

