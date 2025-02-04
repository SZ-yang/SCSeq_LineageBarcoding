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

