# https://github.com/SZ-yang/SCSeq_LineageBarcoding/blob/kevin/kevin/Writeup5_GEMLI/Writeup5_LARRY_gsea-plotting.R

rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)

data_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-hotspot/"

df <- read.csv(paste0(data_folder, "LCL_celltagMulti_hotspot_autocorrelations.csv"))

teststat_vec <- df[,"Z"]
names(teststat_vec) <- df[,"Gene"]
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 1,
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

gse_df <- as.data.frame(gse)
table(gse_df$p.adjust <= 0.05)
table(gse_df$pvalue <= 0.05)

gse_df[which(gse_df$pvalue <= 0.05), "Description"]

########################

df_scvi <- read.csv(paste0(data_folder, "scVI_celltagMulti_hotspot_autocorrelations.csv"))

teststat_vec <- df_scvi[,"Z"]
names(teststat_vec) <- df_scvi[,"Gene"]
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

set.seed(10)
gse_scvi <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 1,
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

gse_df_scvi <- as.data.frame(gse_scvi)
table(gse_df_scvi$p.adjust <= 0.05)
table(gse_df_scvi$pvalue <= 0.05)

gse_df_scvi[which(gse_df_scvi$p.adjust <= 0.05), "Description"]
gse_df_scvi[which(gse_df_scvi$pvalue <= 0.05), "Description"]



