rm(list=ls())

library(Seurat)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))
length(levels(seurat_obj$fine_cluster))

tab_mat_condition <- table(seurat_obj$Lineage,
                           seurat_obj$OG_condition)
tab_mat_cluster <- table(seurat_obj$Lineage,
                         seurat_obj$fine_cluster)

# find all lineages with cis and dabtramtocis
col_idx <- which(colnames(tab_mat_condition) %in% c("cis", "dabtramtocis"))
lineage_idx <- which(apply(tab_mat_condition, 1, function(x){
  all(x[col_idx] > 0)
}))
fine_clusters <- unique(seurat_obj$fine_cluster[which(seurat_obj$Lineage %in% names(lineage_idx))])
unique_lineages <- unique(seurat_obj$Lineage[which(seurat_obj$fine_cluster %in% fine_clusters)])

tab_mat_condition[unique_lineages,]
paste0("Original lineages: ", length(lineage_idx), " vs. new lineages: ", length(unique_lineages))
