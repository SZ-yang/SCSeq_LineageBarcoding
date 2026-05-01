rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/shaffer_clonal-treatment/all_data_final_lineages.RData")

df <- all_data@meta.data

write.csv(df, 
          file = "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup7/Writeup7_dylan_metadata.csv")
