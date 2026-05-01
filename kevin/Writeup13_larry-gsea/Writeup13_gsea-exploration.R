# https://github.com/SZ-yang/Lineage_aware_ContraLearn/blob/master/analysis/CrossEntropy_sup/train_test/LARRY_top200/CE_GSEA_IG_GeneScores.R
# https://github.com/SZ-yang/SCSeq_LineageBarcoding/blob/kevin/kevin/Writeup11_cospar_celltagmulti/Writeup11_cospar_gsea.R
# https://github.com/SZ-yang/SCSeq_LineageBarcoding/blob/kevin/kevin/Writeup10b_joshua-hotspot/Writeup10b_gsea_IG.R

rm(list=ls())

library(clusterProfiler)
library(org.Mm.eg.db)

# 2. Load your data (Change to ce_ig_scores.csv for the baseline)
lcl_df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/csv/kevin/Writeup13/joshua_output/lcl_ig_scores.csv")

# 3. Create a named vector and sort it descending
# This is the exact format clusterProfiler requires
teststat_vec <- lcl_df$Score
names(teststat_vec) <- lcl_df$Gene
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

# 4. Run GSEA
set.seed(10)
lcl_gse <- clusterProfiler::gseGO(
  geneList     = teststat_vec,
  ont          = "BP",           # "BP" = Biological Process pathways
  keyType      = "SYMBOL",       # Assuming LARRY genes are symbols (e.g., "Sca1")
  OrgDb        = "org.Mm.eg.db", # Mouse database
  pvalueCutoff = 1,           # Only keep statistically significant pathways
  minGSSize    = 10,
  maxGSSize    = 500,
  scoreType    = "pos"           # Because your IG scores are absolute/positive
)

# 5. Extract results to a data frame and view the top hits
lcl_gse_df <- as.data.frame(lcl_gse)
lcl_gse_df <- lcl_gse_df[order(lcl_gse_df$p.adjust, decreasing = FALSE), ]
lcl_gse_df[1:50, c("Description", "p.adjust")]

#####################

# 2. Load your data (Change to ce_ig_scores.csv for the baseline)
ce_df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/csv/kevin/Writeup13/joshua_output/ce_ig_scores.csv")

# 3. Create a named vector and sort it descending
# This is the exact format clusterProfiler requires
teststat_vec <- ce_df$Score
names(teststat_vec) <- ce_df$Gene
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

# 4. Run GSEA
set.seed(10)
ce_gse <- clusterProfiler::gseGO(
  geneList     = teststat_vec,
  ont          = "BP",           # "BP" = Biological Process pathways
  keyType      = "SYMBOL",       # Assuming LARRY genes are symbols (e.g., "Sca1")
  OrgDb        = "org.Mm.eg.db", # Mouse database
  pvalueCutoff = 1,           # Only keep statistically significant pathways
  minGSSize    = 10,
  maxGSSize    = 500,
  scoreType    = "pos"           # Because your IG scores are absolute/positive
)

# 5. Extract results to a data frame and view the top hits
ce_gse_df <- as.data.frame(ce_gse)
ce_gse_df <- ce_gse_df[order(ce_gse_df$p.adjust, decreasing = FALSE), ]
ce_gse_df[1:50, c("Description", "p.adjust")]

###############

colnames(lcl_df)[2] <- "lcl_Score"
colnames(ce_df)[2] <- "ce_Score"

merge_gene_df <- merge(lcl_df, ce_df, by = "Gene")
rownames(merge_gene_df) <- merge_gene_df$Gene

plot(merge_gene_df[,2], merge_gene_df[,3])

###############

library(ggplot2)
library(ggrepel)

ggplot(merge_gene_df, aes(x = lcl_Score, y = ce_Score, label = Gene)) +
  # 1. Add a dotted diagonal reference line (y = x)
  # Putting this first keeps it behind the data points
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey40") +
  
  # 2. Add the data points
  geom_point(color = "steelblue", alpha = 0.5) +
  
  # 3. Add the labels with overlap protection
  geom_text_repel(
    # 1. Reduce the padding to 0 or a very small number
    point.padding = 0.1, 
    box.padding = 0.25,   
    
    # 2. Increase the attraction to the data point
    # Higher force_pull (default is 1) makes labels "snap" closer to the points
    force_pull = 5,
    
    # 3. Decrease the repulsion force (default is 1) 
    # to allow them to sit closer to each other
    force = 0.5,
    
    max.overlaps = 20, 
    size = 3,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.alpha = 0.5
  ) +
  
  # 4. Force the 1:1 aspect ratio
  coord_fixed(ratio = 1) +
  
  # 5. Styling
  theme_minimal() +
  labs(
    title = "Gene Score Comparison",
    subtitle = "Dotted line represents y = x",
    x = "lcl_Score",
    y = "ce_Score"
  )


round(merge_gene_df[c("Mpo", "Ngp", "Ltf", "Camp", "Gata2", "Cebpe", "S100a8", "Itgam", "Klf4", "Irf8"), c(2,3)], 3)

