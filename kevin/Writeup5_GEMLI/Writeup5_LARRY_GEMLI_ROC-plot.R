rm(list=ls())

plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup5/"

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Writeup5_Larry_LCL_GEMLI.RData"))
GEMLI_items_LCL <- GEMLI_items

load(paste0(out_folder, "Writeup5_Larry_log_GEMLI.RData"))
GEMLI_items_log <- GEMLI_items

load(paste0(out_folder, "Writeup5_Larry_scVI_GEMLI.RData"))
GEMLI_items_scVI <- GEMLI_items

###############

.compute_tpr_fpr <- function(GEMLI_items){
  total_true <- GEMLI_items$testing_results["0","TP"]
  total_false <- GEMLI_items$testing_results["0","FP"]
  
  tpr <- GEMLI_items$testing_results[,"TP"]/total_true
  fpr <- GEMLI_items$testing_results[,"FP"]/total_false
  
  list(fpr = fpr,
       tpr = tpr)
}

tmp <- .compute_tpr_fpr(GEMLI_items_LCL)
fpr_LCL <- tmp$fpr
tpr_LCL <- tmp$tpr

tmp <- .compute_tpr_fpr(GEMLI_items_log)
fpr_log <- tmp$fpr
tpr_log <- tmp$tpr

tmp <- .compute_tpr_fpr(GEMLI_items_scVI)
fpr_scVI <- tmp$fpr
tpr_scVI <- tmp$tpr

############

# Compute AUC using the trapezoidal rule
.compute_auc <- function(fpr, tpr){
  auc <- -(sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
  return(2*(auc-0.5))
}

LCL_auc <- .compute_auc(fpr_LCL, tpr_LCL)
log_auc <- .compute_auc(fpr_log, tpr_log)
scVI_auc <- .compute_auc(fpr_scVI, tpr_scVI)

############

df <- data.frame(
  tpr = c(tpr_LCL, tpr_log, tpr_scVI),
  fpr = c(fpr_LCL, fpr_log, fpr_scVI),
  model = rep(c("LCL", "log-normalized" , "scVI"),
              each = length(tpr_LCL))
)

library(ggplot2)

# Plot both ROC curves
plot1 <- ggplot(df, aes(x = fpr, y = tpr, color = model)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "red") +  # Diagonal line (random classifier)
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", 
       title = paste0("ROC Curve Comparison",
                      "\nAUC: LCL: ", round(LCL_auc, 2), 
                      ", log: ", round(log_auc, 2), 
                      ", scVI: ", round(scVI_auc, 2))) +
  scale_color_manual(values = c(rgb(212, 63, 136, maxColorValue = 255), 
                                rgb(117, 164, 58, maxColorValue = 255),
                                rgb(221, 173, 59, maxColorValue = 255))) +  # Customize colors if you want
  coord_equal() 

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup5_LARRY_GEMLI_ROC.png"),
                height = 1200, width = 1200, units = "px")


###

# Plot both ROC curves
plot1 <- ggplot(df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(size = 1.5) +
  geom_abline(linetype = "dashed", color = "red") +  # Diagonal line (random classifier)
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c(rgb(212, 63, 136, maxColorValue = 255), 
                                rgb(117, 164, 58, maxColorValue = 255),
                                rgb(221, 173, 59, maxColorValue = 255))) +  # Customize colors if you want
  coord_equal() +
  Seurat::NoLegend()

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup5_LARRY_GEMLI_ROC_cleaned.png"),
                height = 900, width = 900, units = "px")


