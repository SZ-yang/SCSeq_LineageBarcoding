rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)

csv_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup5/"
df <- read.csv(paste0(csv_folder, "Writeup5_LARRY_LCL_hotspot_autocorrelations.csv"))

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

#########################################

module <- read.csv(paste0(csv_folder, "Writeup5_LARRY_LCL_hotspot_module-clustering.csv"))
max_cluster <- max(module$Module)

for(kk in 1:max_cluster){
  print(paste0("Cluster ", kk))
  idx <- which(module$Module == kk)
  genes <- module$X[idx]
  
  writeLines(genes, paste0(csv_folder, "Writeup5_LARRY_LCL_hotspot_module_", kk, ".csv"))
}


