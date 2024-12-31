rm(list=ls())
library(Seurat)

load("/Users/kevinlin/Downloads/all_data_final_lineages.RData")

all_data

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

Seurat::DimPlot(all_data, group.by = "OG_condition")

tab_vec <- table(all_data$Lineage)
mean(tab_vec)
length(tab_vec)
