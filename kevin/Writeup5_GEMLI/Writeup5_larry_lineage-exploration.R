rm(list=ls())

library(Seurat)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))

# which lineages share the same island as the top lineage?

tab_vec <- table(seurat_obj$clone_id)
largest_lineage <- names(tab_vec)[which.max(tab_vec)]
idx <- which(seurat_obj$clone_id == largest_lineage)

umap_res <- seurat_obj[["LCLumap"]]@cell.embeddings

lineage_names <- unique(seurat_obj$clone_id)
lineage_idx_list <- lapply(lineage_names, function(lineage){
  which(seurat_obj$clone_id == lineage)
})
lineage_means <- sapply(lineage_idx_list, function(idx){
  Matrix::colMeans(umap_res[idx,,drop=FALSE])
})
lineage_means <- t(lineage_means)
rownames(lineage_means) <- lineage_names

range_vec <- apply(lineage_means, 2, function(x){diff(range(x))})
tolerance <- mean(range_vec)*0.05

# find all lienages within tolerance
dist_mat <- as.matrix(stats::dist(lineage_means))

idx <- which(dist_mat[largest_lineage,] <= tolerance)
close_lineages <- colnames(dist_mat)[idx]
  
plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup5/"
pdf(paste0(plot_folder, "Writeup5_Larry_41093_2000_norm_log_cleaned_close-lineages.pdf"),
    onefile = T, width = 5, height = 5)

for(i in 1:length(close_lineages)){
  lineage_name <- close_lineages[i]
  idx <- which(seurat_obj$clone_id == lineage_name)
  
  cell_names <- Seurat::Cells(seurat_obj)[idx]
  plot1 <- Seurat::DimPlot(object = seurat_obj, 
                           reduction = "LCLumap",
                           cells.highlight = cell_names, 
                           cols.highlight = "red", 
                           cols = "gray", 
                           order = TRUE,
                           raster = TRUE)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0(lineage_name, " with ", length(idx), " cells"))
  
  print(plot1)
}

dev.off()
  

##############################

unique_celltypes <- sort(unique(seurat_obj$state_info))

composition_mat <- sapply(close_lineages, function(lineage){
  idx <- which(seurat_obj$clone_id == lineage)
  tab_vec <- table(seurat_obj$state_info[idx])
  tab_vec <- tab_vec/sum(tab_vec)
  tab_vec
})
  
round(100*composition_mat)


  