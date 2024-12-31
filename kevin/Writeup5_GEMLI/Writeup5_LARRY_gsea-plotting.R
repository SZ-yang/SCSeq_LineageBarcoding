rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup5/"

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
length(which(gse_df_LCL$p.adjust <= 0.05))

gse_short <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 0.05,
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

plot1 <- enrichplot::dotplot(gse_short, showCategory=nrow(gse_short)) + ggplot2::ggtitle("GSEA on LCL's Hotspot statistic")
plot1 <- plot1 + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

ggplot2::ggsave(filename = paste0(plot_folder, "Writeup5_LARRY_LCL-hotspot_dotplot.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

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

################

df <- read.csv(paste0(csv_folder, "Writeup5_LARRY_scVI_hotspot_day2_autocorrelations.csv"))

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

gse_df_scVI <- as.data.frame(gse)
length(which(gse_df_scVI$p.adjust <= 0.05))

go_scVI <- gse_df_scVI[which(gse_df_scVI$p.adjust <= 0.05),"ID"]
go_LCL <- gse_df_LCL[which(gse_df_LCL$p.adjust <= 0.05),"ID"]
setdiff(go_LCL, go_scVI)
gse_df_LCL[setdiff(go_LCL, go_scVI),1:2]
gse_df_LCL["GO:0002520",]

##################

# now do cospar
cospar_df <- read.csv("~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup4/Writeup4_cospar-LARRY_postprocess_adata-obs.csv")
cospar_df <- cospar_df[cospar_df$time_info == 2,]
quantile(cospar_df$fate_bias_transition_map_Neutrophil.Monocyte)

tmp <- cospar_df[,c("Library", "Cell.barcode")]
tmp$Cell.barcode <- gsub("-", "", tmp$Cell.barcode)
cell_id <- apply(tmp, 1, function(x){paste0(x, collapse = ":")})
fate_vec <- cospar_df$fate_bias_transition_map_Neutrophil.Monocyte
names(fate_vec) <- cell_id

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))
mat <- SeuratObject::LayerData(
  seurat_obj,
  layer = "data",
  assay = "RNA"
)
mat <- mat[,names(fate_vec)]
teststat_vec <- apply(mat, 1, function(x){
  stats::cor(fate_vec, x)
})
teststat_vec <- teststat_vec[!is.na(teststat_vec)]
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

gse_df_cospar <- as.data.frame(gse)
length(which(gse_df_cospar$p.adjust <= 0.05))

##################

gse_df_LCL[which(gse_df_LCL$p.adjust <= 0.05), "Description"]

# find the GO terms
go_terms <- c("response to oxygen-containing compound",
              "defense response to bacterium",
              "response to external stimulus")
go_id <- sapply(go_terms, function(x){
  gse_df_LCL[which(gse_df_LCL$Description == x),"ID"]
})
go_full <- sapply(1:3, function(i){
  paste0(go_id[i], ": ", go_terms[i])
})

LCL_log10pvalue <- -log10(gse_df_LCL[go_id, "pvalue"])
scVI_log10pvalue <- -log10(gse_df_scVI[go_id, "pvalue"])
memory_log10pvalue <- -log10(gse_df_memory[go_id, "pvalue"])
cospar_log10pvalue <- -log10(gse_df_cospar[go_id, "pvalue"])

# Create a data frame in long format
df <- data.frame(
  pathway = rep(go_full, 4),
  log10_p = c(LCL_log10pvalue, scVI_log10pvalue, memory_log10pvalue, cospar_log10pvalue),
  method = rep(c("Hotspot on LCL",
                 "Hotspot on scVI",
                 "GEMLI memory score", 
                 "CoSPAR"), each = 3)
)

# Specify the order of the pathways (from top to bottom in the plot)
df$pathway <- factor(df$pathway, levels = go_full)  

# Specify the order of the methods (the order within each pathway)
df$method <- factor(df$method, levels = c("GEMLI memory score", "CoSPAR", "Hotspot on scVI", "Hotspot on LCL"))

# Plot the horizontal barplot
plot1 <- ggplot(df, aes(x = pathway, y = log10_p, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # "dodge" to compare side-by-side bars
  coord_flip() +  # Flip coordinates to make bars horizontal
  labs(x = "Pathways", y = "-log10(p-value)", 
       title = "-log10(p-value) of Gene Pathways for Two Methods") +
  scale_fill_manual(values = c(rgb(212, 63, 136, maxColorValue = 255), 
                               rgb(117, 164, 58, maxColorValue = 255),
                               rgb(221, 173, 59, maxColorValue = 255),
                               rgb(0, 123, 206, maxColorValue = 255))) +  # Customize bar colors
  theme_minimal() +
  theme(axis.title.y = element_blank())  # Remove the y-axis title (optional)

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup5_LARRY_GSEA-barplots.png"),
                height = 900, width = 3000, units = "px")

####

color_vec <-  c("Hotspot on LCL" = rgb(212, 63, 136, maxColorValue = 255),
                "CoSPAR" = rgb(0, 123, 206, maxColorValue = 255), 
                "GEMLI memory score" = rgb(117, 164, 58, maxColorValue = 255),
                "Hotspot on scVI" = rgb(221, 173, 59, maxColorValue = 255))

# Plot the horizontal barplot
plot1 <- ggplot(df, aes(x = pathway, y = log10_p, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # "dodge" to compare side-by-side bars
  coord_flip() +  # Flip coordinates to make bars horizontal
  labs(x = "", y = "", 
       title = "") +
  scale_fill_manual(values = color_vec) +  # Customize bar colors
  theme_minimal() +
  theme(axis.title.y = element_blank(),  # Remove the y-axis title
        axis.text.y = element_blank(),   # Remove the y-axis text (pathway names)
        axis.ticks.y = element_blank()) +   # Remove the y-axis ticks
  Seurat::NoLegend()

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup5_LARRY_GSEA-barplots_cleaned.png"),
                height = 900, width = 847, units = "px")




