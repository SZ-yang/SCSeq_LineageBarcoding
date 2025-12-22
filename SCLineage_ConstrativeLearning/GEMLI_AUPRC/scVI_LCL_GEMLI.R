############################################
## Fancy combined PR plot (ggplot2)
############################################

# ---- packages ----
library(ggplot2)

# ---- helpers ----
extract_pr <- function(res) {
  recall <- res[, "sensitivity"]
  precision <- res[, "precision"]
  
  ok <- is.finite(recall) & is.finite(precision)
  recall <- recall[ok]
  precision <- precision[ok]
  
  ord <- order(recall)
  recall <- recall[ord]
  precision <- precision[ord]
  
  # enforce endpoints for cleaner curves/integration
  if (length(recall) > 0) {
    if (min(recall) > 0) {
      recall <- c(0, recall)
      precision <- c(precision[1], precision)
    }
    if (max(recall) < 1) {
      recall <- c(recall, 1)
      precision <- c(precision[length(precision)], precision)
    }
  }
  
  data.frame(recall = recall, precision = precision)
}

auprc_trapz <- function(df) {
  # assumes df is sorted by recall
  r <- df$recall
  p <- df$precision
  sum(diff(r) * (p[-1] + p[-length(p)]) / 2, na.rm = TRUE)
}

load_gemli_pr <- function(rdata_path, label) {
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  
  res <- e$GEMLI_items$testing_results
  df <- extract_pr(res)
  df <- df[order(df$recall), , drop = FALSE]
  
  auprc <- auprc_trapz(df)
  baseline <- res["0","TP"] / (res["0","TP"] + res["0","FP"])
  
  df$method <- label
  
  list(df = df, auprc = auprc, baseline = baseline)
}

# ---- paths ----
out_folder <- "/Users/apple/Desktop/KB/SCSeq_LineageBarcoding/SCLineage_ConstrativeLearning/GEMLI_AUPRC"
lcl_path  <- file.path(out_folder, "cellTagMulti_LCL_GEMLI.RData")
scvi_path <- file.path(out_folder, "cellTagMulti_scVI_GEMLI.RData")

# ---- load data ----
lcl  <- load_gemli_pr(lcl_path,  "LCL")
scvi <- load_gemli_pr(scvi_path, "scVI")

df_plot <- rbind(lcl$df, scvi$df)

# baseline should be the same if evaluated on same pairs; take the mean just in case
baseline_line <- mean(c(lcl$baseline, scvi$baseline), na.rm = TRUE)

# labels with AUPRC
legend_labels <- c(
  paste0("LCL (AUPRC=", signif(lcl$auprc, 4), ")"),
  paste0("scVI (AUPRC=", signif(scvi$auprc, 4), ")")
)
names(legend_labels) <- c("LCL", "scVI")

# ---- plot ----
p <- ggplot(df_plot, aes(x = recall, y = precision, color = method)) +
  geom_line(linewidth = 1.1, alpha = 0.95) +
  geom_point(size = 1.3, alpha = 0.85) +
  geom_hline(yintercept = baseline_line, linetype = "dashed", linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_color_manual(values = c("LCL" = "black", "scVI" = "#2C7FB8"),
                     breaks = c("LCL", "scVI"),
                     labels = legend_labels,
                     name = NULL) +
  labs(
    title = "GEMLI Precision–Recall (CellTag multi): LCL vs scVI",
    subtitle = paste0("Baseline precision = ", signif(baseline_line, 3)),
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey80"),
    panel.grid.minor = element_blank()
  )

# ---- save ----
ggsave(
  filename = file.path(out_folder, "cellTagMulti_LCL_vs_scVI_GEMLI_PR_ggplot.png"),
  plot = p,
  width = 7.5, height = 6.0, dpi = 300
)

# Optional: PDF (vector)
ggsave(
  filename = file.path(out_folder, "cellTagMulti_LCL_vs_scVI_GEMLI_PR_ggplot.pdf"),
  plot = p,
  width = 7.5, height = 6.0
)

print(p)

# ---- print numbers ----
cat("LCL  AUPRC:", signif(lcl$auprc, 5), " baseline:", signif(lcl$baseline, 5), "\n")
cat("scVI AUPRC:", signif(scvi$auprc, 5), " baseline:", signif(scvi$baseline, 5), "\n")
