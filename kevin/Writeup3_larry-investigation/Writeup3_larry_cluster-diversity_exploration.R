rm(list=ls())

library(Seurat)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup3/Writeup3_Larry_scBaseEncoderFeat_Z_bs30_tau05.RData")

plot1 <- Seurat::DimPlot(seurat_obj, 
                reduction = "python_X_umap", 
                group.by = "Cell.type.annotation",
                cols = seurat_obj@misc$Cell.type.annotation_colors)
plot1 <- plot1 + Seurat::NoLegend()
plot2 <- Seurat::DimPlot(seurat_obj, 
                reduction = "umap", 
                group.by = "Cell.type.annotation",
                label = TRUE,
                label.size = 2,
                cols = seurat_obj@misc$Cell.type.annotation_colors)
plot2 <- plot2 + Seurat::NoLegend()
plot_all <- cowplot::plot_grid(plotlist = list(plot1, plot2), ncol = 2)

ggplot2::ggsave(plot_all, filename = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup3/Writeup3_umap_celltype.png",
                height = 1200, width = 2000, dpi = 300, units = "px")

# compute a lineage statistic based on a Shannon entropy 

shannon_entropy <- function(vec){
  vec <- vec[vec > 0]
  vec <- vec/sum(vec)
  -sum(vec * log(vec))
}

#################

lineage_names <- sort(unique(seurat_obj$clone_id))
lineage_entropy <- sapply(lineage_names, function(lineage){
  idx <- intersect(which(seurat_obj$clone_id == lineage),
                   which(seurat_obj$Time.point == "t_6"))
  if(length(idx) == 0) return(NA)
  tab_vec <- table(seurat_obj$Cell.type.annotation[idx])
  if("Undifferentiated" %in% names(tab_vec)){
    idx <- which(names(tab_vec) == "Undifferentiated")
    tab_vec <- tab_vec[-idx]
  }
  if(sum(tab_vec) <= 4) return(NA)
  shannon_entropy(tab_vec)
})
names(lineage_entropy) <- lineage_names
table(is.na(lineage_entropy))

png("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup3/Writeup3_entropy_hist.png",
    width = 1500, height = 1200, res = 300, units = "px")
hist(lineage_entropy,
     xlab = "Shannon entropy (one per lineage)",
     main = paste0("Shannon entropy,\nonly for lineages with more than\n4 non-undiff cells at D6\n",
                   length(which(!is.na(lineage_entropy))), " of ", length(lineage_entropy), " non-NA"))
graphics.off()

# compute the hypothetical maximum entropy
num_celltypes <- length(unique(seurat_obj$Cell.type.annotation))-1
max_entropy <- shannon_entropy(rep(1/num_celltypes, length = num_celltypes))
stopifnot(all(max_entropy >= max(lineage_entropy, na.rm = TRUE)))

# map the entropy back to the cells
cell_entropy <- lineage_entropy[seurat_obj$clone_id]
names(cell_entropy) <- Seurat::Cells(seurat_obj)
seurat_obj$entropy <- cell_entropy

set.seed(10)
cells_vec <- sample(Seurat::Cells(seurat_obj))

plot1 <- scCustomize::FeaturePlot_scCustom(seurat_obj,
                                           features = "entropy",
                                           reduction = "python_X_umap",
                                           na_cutoff = 0,
                                           order = FALSE,
                                           cells = cells_vec)
plot2 <- scCustomize::FeaturePlot_scCustom(seurat_obj, 
                                           features = "entropy",
                                           reduction = "umap",
                                           na_cutoff = 0,
                                           order = FALSE,
                                           cells = cells_vec)
plot_all <- cowplot::plot_grid(plotlist = list(plot1, plot2), ncol = 2)

ggplot2::ggsave(plot_all, filename = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup3/Writeup3_umap_entropy.png",
                height = 1200, width = 2000, dpi = 300, units = "px")

# only monocytes and their entropy

entropy_monocyte <- seurat_obj$entropy
entropy_monocyte[seurat_obj$Cell.type.annotation != "Monocyte"] <- NA
seurat_obj$entropy_monocyte <- entropy_monocyte

set.seed(10)
na_idx <- which(is.na(entropy_monocyte))
non_na_idx <- which(!is.na(entropy_monocyte))
ordering <- c(na_idx, sample(non_na_idx))
cells_vec <- Seurat::Cells(seurat_obj)[ordering]

plot1 <- scCustomize::FeaturePlot_scCustom(seurat_obj,
                                           features = "entropy_monocyte",
                                           reduction = "python_X_umap",
                                           na_cutoff = 0,
                                           order = FALSE,
                                           cells = cells_vec)
plot2 <- scCustomize::FeaturePlot_scCustom(seurat_obj, 
                                           features = "entropy_monocyte",
                                           reduction = "umap",
                                           na_cutoff = 0,
                                           order = FALSE,
                                           cells = cells_vec)
plot_all <- cowplot::plot_grid(plotlist = list(plot1, plot2), ncol = 2)

ggplot2::ggsave(plot_all, filename = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup3/Writeup3_umap_entropy-monocyte.png",
                height = 1200, width = 2000, dpi = 300, units = "px")

# only neutrophil and their entropy

entropy_neutrophil <- seurat_obj$entropy
entropy_neutrophil[seurat_obj$Cell.type.annotation != "Neutrophil"] <- NA
seurat_obj$entropy_neutrophil <- entropy_neutrophil

set.seed(10)
na_idx <- which(is.na(entropy_neutrophil))
non_na_idx <- which(!is.na(entropy_neutrophil))
ordering <- c(na_idx, sample(non_na_idx))
cells_vec <- Seurat::Cells(seurat_obj)[ordering]

plot1 <- scCustomize::FeaturePlot_scCustom(seurat_obj,
                                           features = "entropy_neutrophil",
                                           reduction = "python_X_umap",
                                           na_cutoff = 0,
                                           order = FALSE,
                                           cells = cells_vec)
plot2 <- scCustomize::FeaturePlot_scCustom(seurat_obj, 
                                           features = "entropy_neutrophil",
                                           reduction = "umap",
                                           na_cutoff = 0,
                                           order = FALSE,
                                           cells = cells_vec)
plot_all <- cowplot::plot_grid(plotlist = list(plot1, plot2), ncol = 2)

ggplot2::ggsave(plot_all, filename = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup3/Writeup3_umap_entropy-neutrophil.png",
                height = 1200, width = 2000, dpi = 300, units = "px")