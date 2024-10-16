rm(list=ls())

library(Seurat)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))

seurat_obj_day6 <- subset(seurat_obj, Time.point == 6)

# unique lineages
tab_mat <- table(seurat_obj_day6$clone_id, seurat_obj_day6$state_info)
vec <- rowSums(tab_mat[,c("Baso", "Ccr7_DC", "Eos", "Erythroid", "Lymphoid", "Mast", "Meg", "pDC")])
tab_mat <- cbind(tab_mat[,c("Monocyte", "Neutrophil", "Undifferentiated")], vec)

# monocyte_lineage is more than 80% monocyte
# neutrophil_lineage is more than 80% neutrophil
# undecided_lineage is more than 90% in {monocyte, neutrophil, undecided} but not in the first 3

sum_vec <- rowSums(tab_mat)
monocyte_idx <- which(tab_mat[,"Monocyte"] > 0.95*sum_vec)
neutrophil_idx <- which(tab_mat[,"Neutrophil"] > 0.95*sum_vec)
monocyte_undiff_idx <- which(tab_mat[,"Monocyte"] < 0.7*sum_vec)
neutrophil_undiff_idx <- which(tab_mat[,"Neutrophil"] < 0.7*sum_vec)
undiff_idx <- which(tab_mat[,"vec"] < 0.05*sum_vec)
undiff_idx <- intersect(undiff_idx, which(tab_mat[,"Undifferentiated"] < 0.2*sum_vec))
undiff_idx <- intersect(undiff_idx, monocyte_undiff_idx)
undiff_idx <- intersect(undiff_idx, neutrophil_undiff_idx)
undiff_idx <- setdiff(undiff_idx, c(monocyte_idx, neutrophil_idx))

monocyte_lineages <- rownames(tab_mat)[monocyte_idx]
neutrophil_lineages <- rownames(tab_mat)[neutrophil_idx]
undiff_lineages <- rownames(tab_mat)[undiff_idx]

metadata_df <- seurat_obj@meta.data
metadata_df$assignment <- "Other"
metadata_df$assignment[seurat_obj_day6$clone_id %in% monocyte_lineages] <- "Monocyte"
metadata_df$assignment[seurat_obj_day6$clone_id %in% neutrophil_lineages] <- "Neutrophil"
metadata_df$assignment[seurat_obj_day6$clone_id %in% undiff_lineages] <- "Undifferentiated"

table(metadata_df$time_info, metadata_df$assignment)

out_folder2 <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup6/"
write.csv(metadata_df,
          file = paste0(out_folder2, "LARRY_tabulate-lineages.csv"))



