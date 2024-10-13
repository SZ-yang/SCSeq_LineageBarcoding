rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)

csv_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup5/"
df <- read.csv(paste0(csv_folder, "Writeup5_LARRY_LCL_hotspot_day2_autocorrelations.csv"))

teststat_vec <- df[,"Z"]
names(teststat_vec) <- df[,"Gene"]
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

head(as.data.frame(gse)) # nothing

gse_df <- as.data.frame(gse)

write.csv(gse_df, 
          file = paste0(csv_folder, "Writeup5_LARRY_LCL_hotspot_day2_GSEA.csv"))

############

set.seed(10)
gse_full <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 1,
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

gse_full_df <- as.data.frame(gse_full)

write.csv(gse_full_df, 
          file = paste0(csv_folder, "Writeup5_LARRY_LCL_hotspot_day2_GSEA.csv"))


