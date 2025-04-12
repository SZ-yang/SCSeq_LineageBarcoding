rm(list=ls())
library(Seurat)
library(ggplot2)

load("~/kzlinlab/data/larry_hematopoiesis/larry-dataset_KZL.RData")
source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup8_larry-anova/welch_anova.R")
plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup8/"

tab_vec <- table(seurat_obj$assigned_lineage)
kept_lineages <- names(tab_vec)[order(tab_vec, decreasing = TRUE)[1:50]]

seurat_obj <- subset(seurat_obj, assigned_lineage %in% kept_lineages)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj,
                             verbose = FALSE)

gene_embedding <- seurat_obj[["pca"]]@feature.loadings
gene_vec <- intersect(Seurat::VariableFeatures(seurat_obj),
                      rownames(gene_embedding))
data <- tcrossprod(seurat_obj[["pca"]]@cell.embeddings[,1:30],
                   gene_embedding[gene_vec,1:30])

anova_res <- lapply(1:ncol(data), function(j){
  if(j %% floor(ncol(data)/10) == 0) cat('*')
  df <- data.frame(gene_exp = data[,j],
                   lineage = seurat_obj$assigned_lineage)
  df$lineage <- factor(df$lineage)

  df <- .remove_zero_variance(data = df,
                              response_col = "gene_exp",
                              group_col = "lineage")

  if(nrow(df) == 0) return(list(
    statistic    = NA,
    parameter    = c(df_num = NA, df_den = NA),
    p.value      = NA,
    R2_welch     = NA
  ))

  res <- .welch_anova(data = df,
                      response_col = "gene_exp",
                      group_col = "lineage")
})
names(anova_res) <- colnames(data)

pvalue_vec <- sapply(anova_res, function(x){x$p.value})
quantile(pvalue_vec, na.rm=TRUE)

R2_vec <- sapply(anova_res, function(x){x$R2_welch})
quantile(R2_vec, na.rm=TRUE)

results_df <- data.frame(
  gene = names(anova_res),
  pvalue = pvalue_vec,
  R2 = R2_vec
)
results_df <- results_df[order(results_df$R2, decreasing = TRUE),]
results_df$gene <- NULL
results_df <- results_df[!is.na(results_df$R2),,drop = FALSE]

write.csv(results_df,
          file = "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup8/larry-anova-genes.csv",
          quote = FALSE)

#############

pdf(paste0(plot_folder, "Writeup8_larry-anova_highest-R2.pdf"),
    height = 4, width = 10, onefile = TRUE)

for(gene in rownames(results_df)[1:25]){
  df <- data.frame(gene_exp = data[,gene],
                   lineage = seurat_obj$assigned_lineage)
  df$lineage <- factor(df$lineage)

  median_vec <- sapply(unique(df$lineage), function(lineage){
    idx <- which(df$lineage == lineage)
    stats::median(df$gene_exp[idx])
  })
  names(median_vec) <- unique(df$lineage)
  lineage_names <- names(median_vec)[order(median_vec, decreasing = TRUE)]

  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage, y = gene_exp)) +
    geom_violin(trim = TRUE, scale = "width", color = "black") +  # Black border for clarity
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "red") +  # Median line
    scale_x_discrete(limits = lineage_names, guide = guide_axis(angle = 45)) +  # Rotated x-axis labels
    labs(title = paste0("Gene: ", gene, ", R2=", round(results_df[gene,"R2"], 2)))

  print(plot1)
}

graphics.off()

###


pdf(paste0(plot_folder, "Writeup8_larry-anova_lowest-R2.pdf"),
    height = 4, width = 10, onefile = TRUE)

for(gene in rev(rownames(results_df))[1:25]){
  df <- data.frame(gene_exp = data[,gene],
                   lineage = seurat_obj$assigned_lineage)
  df$lineage <- factor(df$lineage)

  median_vec <- sapply(unique(df$lineage), function(lineage){
    idx <- which(df$lineage == lineage)
    stats::median(df$gene_exp[idx])
  })
  names(median_vec) <- unique(df$lineage)
  lineage_names <- names(median_vec)[order(median_vec, decreasing = TRUE)]

  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage, y = gene_exp)) +
    geom_violin(trim = TRUE, scale = "width", color = "black") +  # Black border for clarity
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "red") +  # Median line
    scale_x_discrete(limits = lineage_names, guide = guide_axis(angle = 45)) +  # Rotated x-axis labels
    labs(title = paste0("Gene: ", gene, ", R2=", round(results_df[gene,"R2"], 2)))

  print(plot1)
}

graphics.off()








