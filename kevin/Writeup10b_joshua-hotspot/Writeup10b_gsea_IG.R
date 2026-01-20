# https://github.com/SZ-yang/SCSeq_LineageBarcoding/blob/kevin/kevin/Writeup5_GEMLI/Writeup5_LARRY_gsea-plotting.R

rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)

df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-IG/IG_analysis_res.csv")

teststat_vec <- df[,"score_001_15"]
names(teststat_vec) <- df[,"gene"]
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

gse_df[which(gse_df$p.adjust <= 0.05), "Description"]
