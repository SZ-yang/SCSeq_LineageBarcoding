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
  pvalueCutoff = 1,
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

gse_df_LCL <- as.data.frame(gse)

################

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_GEMLI-memory-genes_day2.RData"))

teststat_vec <- GEMLI_items$memory_genes[,"var"]
names(teststat_vec) <- rownames(GEMLI_items$memory_genes)
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

gse_df_memory <- as.data.frame(gse)

##################

gse_df_LCL[which(gse_df_LCL$p.adjust <= 0.05), "Description"]

# find the GO terms
go_terms <- c("immune system process",
              "response to external stimulus",
              "response to oxygen-containing compound")
go_id <- sapply(go_terms, function(x){
  gse_df_LCL[which(gse_df_LCL$Description == x),"ID"]
})
go_full <- sapply(1:3, function(i){
  paste0(go_id[i], ": ", go_terms[i])
})

lcl_log10pvalue <- -log10(gse_df_LCL[go_id, "pvalue"])
memory_log10pvalue <- -log10(gse_df_memory[go_id, "pvalue"])


# Create a data frame in long format
df <- data.frame(
  pathway = rep(go_full, 2),
  log10_p = c(lcl_log10pvalue, memory_log10pvalue),
  method = rep(c("Hotspot on LCL", "GEMLI memory score"), each = 3)
)

library(ggplot2)

# Plot the horizontal barplot
plot1 <- ggplot(df, aes(x = pathway, y = log10_p, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # "dodge" to compare side-by-side bars
  coord_flip() +  # Flip coordinates to make bars horizontal
  labs(x = "Pathways", y = "-log10(p-value)", 
       title = "-log10(p-value) of Gene Pathways for Two Methods") +
  scale_fill_manual(values = c(rgb(212, 63, 136, maxColorValue = 255), 
                               rgb(117, 164, 58, maxColorValue = 255))) +  # Customize bar colors
  theme_minimal() +
  theme(axis.title.y = element_blank())  # Remove the y-axis title (optional)


plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup5/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup5_LARRY_GSEA-barplots.png"),
                height = 900, width = 3000, units = "px")



