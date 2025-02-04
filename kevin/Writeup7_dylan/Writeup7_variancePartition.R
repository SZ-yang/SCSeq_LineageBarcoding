rm(list=ls())

library(Seurat)
library(variancePartition)

data_folder <- "~/kzlinlab/data/shaffer_clonal-treatment/"
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/"
fig_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup7/"

load(paste0(out_folder, "Writeup7_shaffer_preprocessed.RData"))

gene_vec <- Seurat::VariableFeatures(all_data)

form <- ~ (1 | Lineage) + (1 | OG_condition) + percent.mt + G2M.Score

scaled_data <- SeuratObject::LayerData(all_data,
                                       assay = "RNA",
                                       layer = "scale.data",
                                       features = gene_vec)

# subset for top variable genes
var_data <- scaled_data
var_info <- seurat_obj@meta.data

varPart <- variancePartition::fitExtractVarPartModel(var_data,
                                                     form,
                                                     var_info,
                                                     showWarnings = FALSE)

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(varPart, date_of_run, session_info,
     file = paste0(out_folder, "Writeup7_shaffer_variancePartition.RData"))

print("Done! :)")






