rm(list=ls())

library(EnhancedVolcano)
library(ggplot2)

plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup5/"

csv_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/csv/kevin/Writeup5/"
df <- read.csv(paste0(csv_folder, "Writeup5_LARRY_LCL_hotspot_day2_autocorrelations.csv"))

res <- df
rownames(res) <- res$Gene
pval_vec <- df[,"Pval"]
pval_vec[pval_vec == 0] <- min(pval_vec[pval_vec!=0])
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)
#define thresholds
pCutoff <- max(pval_vec[idx])
FCcutoff <- quantile(abs(res[,"C"]), probs = 0.9)
xlim <- quantile(res[,"C"], probs = c(0.01,0.99))

idx <- intersect(which(res[,"C"] > xlim[1]),
                 which(res[,"C"] > xlim[2]))
ylim <- c(0, max(-log10(pval_vec[idx])))

plot1 <- EnhancedVolcano::EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "C",
  y = "Pval",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  xlim = xlim,
  ylim = ylim
)
plot1 <- plot1 + ggplot2::labs(x = "Hotspot autocorrelation statistic",
                               title = "LCL's Hotspot statistic") 

ggplot2::ggsave(filename = paste0(plot_folder, "Writeup5_LARRY_LCL-hotspot_volcano.png"),
                plot1, device = "png", width = 6, height = 6, units = "in")
