rm(list=ls())

library(Seurat)
library(scCustomize)

fig_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup10/"


# load("~/kzlinlab/data/celltagging-multi_fibroblast/celltagging-multi_fibroblast.RData")
# ctM_metadata <- seurat_obj@meta.data
# write.csv(ctM_metadata, 
#           file = "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup10/train_celltagMulti_metadata.csv")
ctM_metadata <- read.csv("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup10/train_celltagMulti_metadata.csv",
                         row.names = 1)

# load("~/kzlinlab/data/biddy_2018_celltag/biddy_seurat.RData")
# ct_metadata <- seurat_obj@meta.data
# write.csv(ct_metadata, 
#           file = "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup10/test_celltag_metadata.csv")
ct_metadata <- read.csv("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup10/test_celltag_metadata.csv",
                        row.names = 1)

rm(list="seurat_obj")

ctM_lcl <- read.csv("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup10/train_proj_embe_w_clone_id.csv",
                    row.names = 1)
ct_lcl <- read.csv("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup10/test_proj_embe_w_clone_id.csv",
                   row.names = 1)

### checks

table(rownames(ctM_lcl) %in% rownames(ctM_metadata))
table(rownames(ct_lcl) %in% rownames(ct_metadata))

#########

colnames(ctM_metadata)
# [1] "orig.ident"                         "nCount_RNA"                        
# [3] "nFeature_RNA"                       "percent.mt"                        
# [5] "S.Score"                            "G2M.Score"                         
# [7] "Phase"                              "old.ident"                         
# [9] "CC.Difference"                      "RNA_snn_res.0.8"                   
# [11] "seurat_clusters"                    "sample"                            
# [13] "replicate"                          "predicted.id_cca_co"               
# [15] "prediction.score.Fib_1_cca_co"      "prediction.score.Fib_0_cca_co"     
# [17] "prediction.score.Fib_2_cca_co"      "prediction.score.Early_0_cca_co"   
# [19] "prediction.score.Tran_0_cca_co"     "prediction.score.Tran_1_cca_co"    
# [21] "prediction.score.Early_1_cca_co"    "prediction.score.Early_2_cca_co"   
# [23] "prediction.score.iEP_1_cca_co"      "prediction.score.Tran_2_cca_co"    
# [25] "prediction.score.iEP_2_cca_co"      "prediction.score.Dead.end_1_cca_co"
# [27] "prediction.score.Dead.end_0_cca_co" "prediction.score.iEP_0_cca_co"     
# [29] "prediction.score.Dead.end_2_cca_co" "prediction.score.max_cca_co"       
# [31] "RNA_snn_res.0.2"                    "cellranger_ident"                  
# [33] "md_fate_rev1"                       "md_fate_coarse_rev1"               
# [35] "cell_name"                          "assigned_lineage"                  
# [37] "has_lineage"                       

colnames(ct_metadata)
# [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"        
# [4] "nCount_originalexp"   "nFeature_originalexp" "timecourse"          
# [7] "reprogramming_day"    "reprogramming"        "cell_type"           
# [10] "cell_cycle"           "cluster"              "monocle_state"       
# [13] "pseudotime"           "CellTagD0_85k"        "CellTagD3_85k"       
# [16] "CellTagD13_85k"       "CellTagD0_48k"        "CellTagD3_48k"       
# [19] "CellTagD13_48k"       "keep"                 "percent.mt"          
# [22] "S.Score"              "G2M.Score"            "Phase"               
# [25] "nCount_SCT"           "nFeature_SCT" 

#################

ctM_clone <- ctM_lcl[,ncol(ctM_lcl)]
ct_clone <- ct_lcl[,ncol(ct_lcl)]
ctM_lcl <- ctM_lcl[,-ncol(ctM_lcl)]
ct_lcl <- ct_lcl[,-ncol(ct_lcl)]
ctM_clone <- paste0("ctM_", ctM_clone)
ct_clone <- paste0("ct_", ct_clone)

ctM_metadata <- ctM_metadata[rownames(ctM_lcl),]
ct_metadata <- ct_metadata[rownames(ct_lcl),]

# create Seurat
df1 <- cbind(ctM_metadata[,c("sample", "predicted.id_cca_co")], ctM_clone, "celltag-Multi")
df2 <- cbind(ct_metadata[,c("reprogramming_day", "cell_type")], ct_clone, "celltag")
colnames(df1) <- c("time", "cell_type", "clone", "dataset")
colnames(df2) <- c("time", "cell_type", "clone", "dataset")

metadata <- rbind(df1, df2)
for(j in 1:ncol(metadata)){
  metadata[,j] <- factor(metadata[,j])
}
lcl_embedding <- rbind(ctM_lcl, ct_lcl)
lcl_embedding <- as.matrix(lcl_embedding)
colnames(lcl_embedding) <- paste0("gene_", 1:ncol(lcl_embedding))
seurat_obj <- Seurat::CreateSeuratObject(counts = Matrix::t(lcl_embedding),
                                         meta.data = metadata)
colnames(lcl_embedding) <- paste0("LCL_", 1:ncol(lcl_embedding))
seurat_obj[["lcl"]] <- Seurat::CreateDimReducObject(lcl_embedding)

# Run UMAPs
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "lcl", 
                              dims = 1:ncol(lcl_embedding))

plot1 <- scCustomize::DimPlot_scCustom(
  seurat_obj, 
  group.by = "dataset"
)

ggplot2::ggsave(plot1,
                filename = paste0(fig_folder, "Writeup10_lcl_umap.png"),
                width = 12, height = 12)

##########################

unique(seurat_obj$cell_type)

celltype_cols <- c(
  "Ambiguous"   = "gray1",
  Fib_0       = "#2A7ABF",
  Fib_1       = "#5AA6E6",
  Fib_2       = "#0B4F8A",
  Fibroblast  = "#1F77B4",
  Early_0     = "#2FBF9B",
  Early_1     = "#66D2B8",
  Early_2     = "#0F8F73",
  Tran_0      = "#FF8C2A",
  Tran_1      = "#FFB066",
  Tran_2      = "#D96B00",
  iEP_0       = "#8E63C7",
  iEP_1       = "#B08AE0",
  iEP_2       = "#6A3D9A",
  iEP         = "#9467BD",
  "Dead-end_0" = "#D62728",
  "Dead-end_1" = "#FF5A5F",
  "Dead-end_2" = "#A50F15"
)

plot1 <- scCustomize::DimPlot_scCustom(
  seurat_obj, 
  group.by = "cell_type",
  colors_use = celltype_cols
)

ggplot2::ggsave(plot1,
                filename = paste0(fig_folder, "Writeup10_lcl_umap_celltype.png"),
                width = 12, height = 14)

table(seurat_obj$cell_type, seurat_obj$dataset)

###########

# largest 5 celltag lineages
tab_vec <- table(seurat_obj@meta.data[seurat_obj$dataset == "celltag","clone"])
tab_vec <- sort(tab_vec, decreasing = TRUE)

lineage_colors <- rep("gray", length = length(tab_vec))
names(lineage_colors) <- unique(seurat_obj$clone)
lineage_colors[names(tab_vec)[1:5]] <- scales::hue_pal()(5)

plot1 <- scCustomize::DimPlot_scCustom(
  seurat_obj, 
  group.by = "clone",
  colors_use = lineage_colors
) + Seurat::NoLegend()

ggplot2::ggsave(plot1,
                filename = paste0(fig_folder, "Writeup10_lcl_umap_clone_top5-test.png"),
                width = 10, height = 10)

#############

# largest 5 celltag-multi lineages
tab_vec <- table(seurat_obj@meta.data[seurat_obj$dataset == "celltag-Multi","clone"])
tab_vec <- sort(tab_vec, decreasing = TRUE)

lineage_colors <- rep("gray", length = length(tab_vec))
names(lineage_colors) <- unique(seurat_obj$clone)
lineage_colors[names(tab_vec)[1:5]] <- scales::hue_pal()(5)

plot1 <- scCustomize::DimPlot_scCustom(
  seurat_obj, 
  group.by = "clone",
  colors_use = lineage_colors
) + Seurat::NoLegend()

ggplot2::ggsave(plot1,
                filename = paste0(fig_folder, "Writeup10_lcl_umap_clone_top5-train.png"),
                width = 10, height = 10)


