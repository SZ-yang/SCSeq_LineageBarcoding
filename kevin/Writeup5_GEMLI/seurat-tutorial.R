# from https://satijalab.org/seurat/articles/de_vignette

# devtools::install_github('satijalab/seurat-data')

library(Seurat)
library(SeuratData)
library(ggplot2)
SeuratData::AvailableData()
SeuratData::InstallData("cbmc")
cbmc <- SeuratData::LoadData("cbmc")

# Normalize the data
cbmc <- NormalizeData(cbmc)

# Find DE features between CD16 Mono and CD1 Mono
Idents(cbmc) <- "rna_annotations"
monocyte.de.markers <- FindMarkers(cbmc, 
                                   ident.1 = "Memory CD4 T", 
                                   ident.2 = "CD14+ Mono")
# view results
head(monocyte.de.markers)