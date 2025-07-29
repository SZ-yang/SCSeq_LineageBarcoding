# Full modified script with proper PRC / F1 / AUPRC handling

rm(list=ls())

# --- Libraries ---
library(Seurat)
library(GEMLI)
library(devtools)   # for session_info()

# --- Paths & Data Load ---
out_folder <- "/Users/apple/Downloads/"
load(file.path(out_folder, "Larry_41093_2000_norm_log_cleaned.RData"))

date_of_run  <- Sys.time()
session_info <- session_info()

code_folder <- "/Users/apple/Desktop/KB/SCSeq_LineageBarcoding2/SCSeq_LineageBarcoding/SCLineage_ConstrativeLearning/out/plot/feat_LCL_2025/GEMLI_custom/"
source(file.path(code_folder, "predict_lineages_custom.R"))
source(file.path(code_folder, "quantify_clusters_iterative_custom.R"))
source(file.path(code_folder, "test_lineages_custom.R"))

# --- Filter to Top 50 Lineages ---
lineage_tab      <- table(seurat_obj$clone_id)
largest_lineages <- names(lineage_tab)[order(lineage_tab, decreasing=TRUE)[1:50]]
keep_vec         <- rep(FALSE, length(Cells(seurat_obj)))
keep_vec[seurat_obj$clone_id %in% largest_lineages] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj      <- subset(seurat_obj, keep==TRUE)

# --- Run GEMLI on LCL embeddings ---
data_matrix                     <- t(seurat_obj[["LCL"]]@cell.embeddings)
GEMLI_items                     <- list(
  gene_expression = data_matrix,
  barcodes        = seurat_obj$clone_id
)

set.seed(10)
GEMLI_items <- predict_lineages_custom(
  GEMLI_items,
  bool_find_markers    = FALSE,
  desired_cluster_size = c(50,200),
  verbose              = 1
)

# Save intermediate
save(date_of_run, session_info, GEMLI_items,
     file = file.path(out_folder, "Writeup5_Larry_LCL_GEMLI.RData"))

# Quick QC of co-clustering scores
zz <- GEMLI_items$prediction
print(quantile(zz[zz != 0]))

# --- Evaluate Predictions ---
GEMLI_items <- test_lineages_custom(
  GEMLI_items,
  valid_fam_sizes = c(50,200)
)
print(GEMLI_items$testing_results)

# --- Prepare Plot Folder ---
plot_folder <- "/Users/apple/Desktop/KB/SCSeq_LineageBarcoding2/SCSeq_LineageBarcoding/SCLineage_ConstrativeLearning/out/plot/feat_LCL_2025/GEMLI_custom/"
if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive=TRUE)

# --- Extract & Clean Metrics ---
res       <- GEMLI_items$testing_results
precision <- res[,"precision"]
recall    <- res[,"sensitivity"]

# 1) Clamp non-finite → 0
precision[!is.finite(precision)] <- 0
recall   [!is.finite(recall)]    <- 0

# 2) Drop thresholds with no predictions
pred_counts <- res[,"TP"] + res[,"FP"]
keep_idx    <- pred_counts > 0
res         <- res[keep_idx, ]

precision <- res[,"precision"]
recall    <- res[,"sensitivity"]

# 3) Compute F1 per threshold
F1 <- ifelse(
  (precision + recall) == 0,
  0,
  2 * precision * recall / (precision + recall)
)
res <- cbind(res, F1 = F1)

# 4) Build PR-dataframe (with optional anchor at (0,1))
dfpr <- data.frame(recall = recall, precision = precision)
dfpr <- rbind(data.frame(recall=0, precision=1), dfpr)
dfpr <- dfpr[order(dfpr$recall), ]

# 5) Compute AUPRC via trapezoid rule
R     <- dfpr$recall
P     <- dfpr$precision
auprc <- sum((R[-1] - R[-length(R)]) * (P[-1] + P[-length(P)]) / 2)
message("AUPRC = ", round(auprc, 4))

# 6) Save augmented results
GEMLI_items$testing_results <- res

# 7) Plot Precision–Recall Curve
prc_file <- file.path(plot_folder, "LARRY_LCL_GEMLI_PRC.png")
png(prc_file, height=1200, width=1200, units="px", res=300)
plot(dfpr$recall, dfpr$precision,
     type="b", pch=16, lwd=2,
     xlab="Recall", ylab="Precision",
     main=sprintf("PR Curve (AUPRC=%.3f)", auprc))
grid()
dev.off()

# --- Plot ROC Curve ---
total_true  <- res["0","TP"]
total_false <- res["0","FP"]
tpr <- res[,"TP"] / total_true
fpr <- res[,"FP"] / total_false

roc_file <- file.path(plot_folder, "LARRY_LCL_GEMLI_ROC.png")
png(roc_file, height=1200, width=1200, units="px", res=300)
plot(fpr, tpr,
     type="b", pch=16, lwd=2,
     xlab="FPR", ylab="TPR",
     main="ROC Curve")
lines(c(0,1), c(0,1), col=2, lwd=2, lty=2)
grid()
dev.off()