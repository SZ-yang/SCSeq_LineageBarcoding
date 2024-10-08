# from https://github.com/UPSUTER/GEMLI
rm(list=ls())
# library(devtools); devtools::install_github("UPSUTER/GEMLI", subdir="GEMLI_package_v0")
library(GEMLI)
library(igraph)
library(tidyverse)
library(UpSetR)
library(Seurat)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup5/"
load(paste0(out_folder, 'GEMLI_example_data_matrix.RData'))
load(paste0(out_folder, 'GEMLI_example_barcode_information.RData'))

dim(data_matrix)
table(lineage_dict_bc)
length(unique(lineage_dict_bc))

GEMLI_items = list()
GEMLI_items[['gene_expression']] = data_matrix
GEMLI_items[['barcodes']] = lineage_dict_bc

GEMLI_items[['gene_expression']][9:14,1:5]
GEMLI_items[['barcodes']][1:5]

GEMLI_items = GEMLI::predict_lineages(GEMLI_items)
GEMLI_items[['prediction']][1:5,15:19]

GEMLI_items = GEMLI::test_lineages(GEMLI_items)
GEMLI_items$testing_results
GEMLI_items = GEMLI::test_lineages(GEMLI_items, plot_results=T)

GEMLI::visualize_as_network(GEMLI_items, cutoff=90) # top image
GEMLI::visualize_as_network(GEMLI_items, cutoff=50) # lower image

GEMLI::visualize_as_network(GEMLI_items, cutoff=90, ground_truth=T, highlight_FPs=T) # top image
GEMLI::visualize_as_network(GEMLI_items, cutoff=50, ground_truth=T, highlight_FPs=T) # lower image

GEMLI_items = GEMLI::prediction_to_lineage_information(GEMLI_items, cutoff=50)
GEMLI_items$predicted_lineage_table[1:5,]
GEMLI_items$predicted_lineages[1:5]

GEMLI::suggest_network_trimming_to_size(GEMLI_items, max_size=2, cutoff=50) # left image
GEMLI_items_post_processed = GEMLI::trim_network_to_size(GEMLI_items, max_size=2, cutoff=50)
GEMLI::visualize_as_network(GEMLI_items_post_processed, cutoff=50) # right image

##################
##################
##################

load(paste0(out_folder, 'GEMLI_crypts_example_data_matrix.RData'))
load(paste0(out_folder, 'GEMLI_crypts_example_barcode_information.RData'))

GEMLI_items_crypts = list()
GEMLI_items_crypts[['prediction']] = Crypts
GEMLI_items_crypts[['barcodes']] = Crypts_bc_dict

GEMLI::visualize_as_network(GEMLI_items_crypts, cutoff=70, display_orphan=F, max_edge_width=1, ground_truth=T, include_labels=F, layout_style="fr") # first image
GEMLI::visualize_as_network(GEMLI_items_crypts, cutoff=70, display_orphan=F, max_edge_width=1, ground_truth=T, include_labels=F, layout_style="kk") # second image
GEMLI::visualize_as_network(GEMLI_items_crypts, cutoff=70, display_orphan=F, max_edge_width=1, ground_truth=T, include_labels=F, layout_style="grid") # third image

load(paste0(out_folder, 'GEMLI_crypts_example_cell_type_annotation.RData'))
GEMLI_items_crypts[['cell_type']] = Crypts_annotation

GEMLI::visualize_as_network(GEMLI_items_crypts, cutoff=70, max_edge_width=5, display_orphan=F, include_labels=F, ground_truth=T, highlight_FPs=T, layout_style="kk", cell_type_colors=T)

cell.type <- unique(GEMLI_items_crypts[['cell_type']]$cell.type)
color <- c("#5386BD", "skyblue1", "darkgreen", "gold", "red", "darkred", "black")
Cell_type_color <- data.frame(cell.type, color)
GEMLI_items_crypts[['cell_type_color']] = Cell_type_color

GEMLI::visualize_as_network(GEMLI_items_crypts, cutoff=70, max_edge_width=5, display_orphan=F, include_labels=F, ground_truth=T, highlight_FPs=T, layout_style="kk", cell_type_colors=T)

GEMLI_items_crypts = GEMLI::prediction_to_lineage_information(GEMLI_items_crypts, cutoff=50)
GEMLI::cell_type_composition_plot(GEMLI_items_crypts, cell_type_colors=T, type=c("bubble"))

GEMLI::cell_type_composition_plot(GEMLI_items_crypts, ground_truth=F, cell_type_colors=T, type=c("upsetR")) 
GEMLI::cell_type_composition_plot(GEMLI_items_crypts, ground_truth=F, cell_type_colors=T, type=c("plain"))

############################
############################
############################

load(paste0(out_folder, 'GEMLI_cancer_example_norm_count.RData'))
load(paste0(out_folder, 'GEMLI_cancer_example_predicted_lineages.RData'))
load(paste0(out_folder, 'GEMLI_cancer_example_cell_type_annotation.RData'))

GEMLI_items = list()
GEMLI_items[['gene_expression']] = Cancer_norm_count
GEMLI_items[['predicted_lineage_table']] = Cancer_predicted_lineages
GEMLI_items[['cell_type']] = Cancer_annotation

GEMLI::cell_type_composition_plot(GEMLI_items, type=c("plain"))
GEMLI::cell_type_composition_plot(GEMLI_items, type=c("upsetR"))

GEMLI_items <- GEMLI::extract_cell_fate_lineages(GEMLI_items, selection=c("inv_tumor", "DCIS"), unique=FALSE, threshold=c(10,10))
GEMLI_items[['cell_fate_analysis']][1:10,] 
table(GEMLI_items[['cell_fate_analysis']]$cell.fate)

GEMLI_items <- GEMLI::cell_fate_DEG_calling(GEMLI_items, 
                                            ident1="sym_DCIS", 
                                            ident2="asym_DCIS", 
                                            min.pct=0.05, 
                                            logfc.threshold=0.1)

###

GEMLI_Seurat <- CreateSeuratObject(GEMLI_items[['gene_expression']], project = "SeuratProject", assay = "RNA")
Metadata <- GEMLI_items[['cell_fate_analysis']]; 
Metadata$ident <- "blank"
ident1="sym_DCIS"
ident2="asym_DCIS"
Metadata$ident[Metadata$cell.fate %in% ident1]<-"ident1"
Metadata$ident[Metadata$cell.fate %in% ident2]<-"ident2"
Meta <- as.data.frame(Metadata[,c(5)])
rownames(Meta) <- Metadata$cell.ID; colnames(Meta)<-c("cell.fate")
GEMLI_Seurat <- AddMetaData(GEMLI_Seurat, Meta, col.name = NULL)
DefaultAssay(object = GEMLI_Seurat) <- "RNA"
Idents(GEMLI_Seurat) <- GEMLI_Seurat$cell.fate
# the following is the missing line
SeuratObject::LayerData(GEMLI_Seurat, 
                        layer = "data",
                        assay = "RNA") <- SeuratObject::LayerData(GEMLI_Seurat, 
                                                                  layer = "counts",
                                                                  assay = "RNA")
DEG <- FindMarkers(object = GEMLI_Seurat, 
                   ident.1 = "ident1", 
                   ident.2 = "ident2", 
                   min.pct=0.05, 
                   logfc.threshold=0.1)



GEMLI::DEG_volcano_plot(GEMLI_items, name1="Sym_DCIS", name2="Asym_DCIS")