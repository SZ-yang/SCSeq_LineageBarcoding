rm(list=ls())

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Writeup5_Larry_LCL_GEMLI.RData"))
GEMLI_items_LCL <- GEMLI_items

load(paste0(out_folder, "Writeup5_Larry_log_GEMLI.RData"))
GEMLI_items_log <- GEMLI_items

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

############

df <- data.frame(
  tpr = c(tpr_LCL, tpr_log),
  fpr = c(fpr_LCL, fpr_log),
  model = rep(c("LCL", "log-normalized"),
              each = length(tpr_LCL))
)

library(ggplot2)

# Plot both ROC curves
plot1 <- ggplot(df, aes(x = fpr, y = tpr, color = model)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "red") +  # Diagonal line (random classifier)
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", 
       title = "ROC Curve Comparison") +
  scale_color_manual(values = c(rgb(212, 63, 136, maxColorValue = 255), 
                                rgb(117, 164, 58, maxColorValue = 255))) +  # Customize colors if you want
  coord_equal()

plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup5/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup5_LARRY_GEMLI_ROC.png"),
                height = 1200, width = 1200, units = "px")


