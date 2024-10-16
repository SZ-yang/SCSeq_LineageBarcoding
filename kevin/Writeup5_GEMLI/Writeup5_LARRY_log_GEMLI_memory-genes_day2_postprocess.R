rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_GEMLI-memory-genes_day2.RData"))

head(GEMLI_items$memory_genes)
quantile(GEMLI_items$memory_genes[,"var"])

teststat_vec <- GEMLI_items$memory_genes[,"var"]
names(teststat_vec) <- rownames(GEMLI_items$memory_genes)
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

gse_df <- as.data.frame(gse)
head(gse_df)

#################

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

csv_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup5/"
write.csv(gse_full_df, 
          file = paste0(csv_folder, "Writeup5_LARRY_log_GEMLI_memory_day2_GSEA.csv"))


