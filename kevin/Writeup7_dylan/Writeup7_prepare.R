rm(list=ls())
library(Seurat)

data_folder <- "~/kzlinlab/data/shaffer_clonal-treatment/"
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/"
fig_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup7/"

load(paste0(data_folder, "all_data_final_lineages.RData"))

keep_vec <- rep(FALSE, length(Seurat::Cells(all_data)))
keep_vec[grep("^Lin*", all_data$Lineage)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

tab_vec <- table(all_data$Lineage)
lineage_names <- names(tab_vec)[tab_vec >= 5]
lineage_names <- setdiff(lineage_names, "Lin221115")
keep_vec <- rep(FALSE, length(Seurat::Cells(all_data)))
keep_vec[all_data$Lineage %in% lineage_names] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)


all_data

variable_vec <- c("OG_condition", "Phase", "RNA_snn_res.0.5", "Lineage")
for(variable in variable_vec){
  all_data@meta.data[,variable] <- factor(all_data@meta.data[,variable])
}

summary(all_data@meta.data)

################

var_vec <- c(Lineage = FALSE,
             OG_condition = TRUE)

pdf(paste0(fig_folder, "Writeup7_umap-covariates_original.pdf"),
    onefile = T, width = 7, height = 5)

for(kk in 1:length(var_vec)){
  plot1 <- Seurat::DimPlot(all_data, 
                           reduction = "umap",
                           group.by = names(var_vec)[kk],
                           raster = TRUE)
  if(!var_vec[kk]) plot1 <- plot1 + Seurat::NoLegend()
  print(plot1)
}

dev.off()

tab_vec <- table(all_data$Lineage)
mean(tab_vec)
quantile(tab_vec)
length(tab_vec)

table(all_data$OG_condition)

tab_mat <- table(all_data$Lineage,all_data$OG_condition)
row_sum <- rowSums(tab_mat)
tab_mat <- tab_mat[order(row_sum, decreasing = TRUE),]
tab_mat[1:25,]

######################

all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, 
                           features = Seurat::VariableFeatures(object = all_data), verbose = FALSE)
all_data <- Seurat::RunUMAP(all_data, dims = 1:30, verbose = FALSE)

color_vec <- c("cis" = rgb(202,111,55, maxColorValue = 255),
               "cistocis" = rgb(248,210,152, maxColorValue = 255),
               "cistococl2" = rgb(240,148,71, maxColorValue = 255),
               "cistodabtram" = rgb(160,61,38, maxColorValue = 255),
               "cocl2" = rgb(69,132,69, maxColorValue = 255),
               "cocl2tocis" = rgb(131,202,163, maxColorValue = 255),
               "cocl2tococl2" = rgb(126,191,90, maxColorValue = 255),
               "cocl2todabtram" = rgb(35,63,58, maxColorValue = 255),
               "dabtram" = rgb(68,49,147, maxColorValue = 255),
               "dabtramtocis" = rgb(147,137,193, maxColorValue = 255),
               "dabtramtococl2" = rgb(145,54,147, maxColorValue = 255),
               "dabtramtodabtram" = rgb(68,32,85, maxColorValue = 255))

plot1 <- Seurat::DimPlot(all_data, 
                         cols = color_vec, 
                         group.by = "OG_condition")
ggplot2::ggsave(plot1, file = paste0(fig_folder, "Writeup7_umap_OG_condition.png"),
                width = 7, height = 5)

plot1 <- Seurat::DimPlot(all_data, 
                         group.by = "Lineage") + Seurat::NoLegend()
ggplot2::ggsave(plot1, file = paste0(fig_folder, "Writeup7_umap_Lineage.png"),
                width = 7, height = 5)

######################

save(all_data,
     file = paste0(out_folder, "Writeup7_shaffer_preprocessed.RData"))

####

# export the count matrix

library(Seurat)
library(SeuratObject)
library(Matrix)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/"

# 1. Read the obs (cell) and var (gene) names from Python
obs_python <- read.csv(paste0(out_folder, "obs_names_python.csv"), header=FALSE, stringsAsFactors=FALSE)[,1]
var_python <- read.csv(paste0(out_folder, "var_names_python.csv"), header=FALSE, stringsAsFactors=FALSE)[,1]

# 2. Subset your Seurat object to match these cells and genes
#    In Seurat, cells are columns and features (genes) are rows.
#    Make sure you only keep the intersection to avoid missing data errors.
common_cells <- intersect(Seurat::Cells(all_data), obs_python)
common_genes <- intersect(SeuratObject::Features(all_data), var_python)

# Subset the Seurat object to keep only the common cells and genes
all_data_subset <- subset(all_data, cells = common_cells, features = common_genes)

# 3. Extract the raw count matrix (layer = "counts") 
#    If you only stored counts in the default assay slot, you could do GetAssayData(..., slot="counts").
#    But since you specifically mention LayerData:
zz <- SeuratObject::LayerData(all_data_subset,
                              layer = "counts", 
                              assay = "RNA",
                              features = common_genes)

# 4. Reorder rows and columns to match EXACTLY the Python var/cell order
#    By default, LayerData might return features x cells. Check dimension names with rownames() / colnames().
#    Confirm which dimension is genes vs cells, then reorder:

# Genes in rows
zz <- zz[match(var_python, rownames(zz)), ]   # Reorder rows to match var_python
# Cells in columns
zz <- zz[, match(obs_python, colnames(zz))]   # Reorder columns to match obs_python

# 5. Finally, write out the matrix in .mtx format
writeMM(zz, file = paste0(out_folder, "counts_matrix.mtx"))

# Also optionally export row and column names if you want them in separate files:
writeLines(rownames(zz), paste0(out_folder, "genes_ordered.txt"))
writeLines(colnames(zz), paste0(out_folder, "cells_ordered.txt"))

