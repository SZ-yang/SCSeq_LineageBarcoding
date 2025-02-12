rm(list=ls())

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/"
plot_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup7/"

load(paste0(out_folder, "Writeup7_GEMLI_pca.RData"))
GEMLI_items_pca <- GEMLI_items

load(paste0(out_folder, "Writeup7_GEMLI_lcl.RData"))
GEMLI_items_lcl <- GEMLI_items

###############

.compute_tpr_fpr <- function(GEMLI_items){
  total_true <- GEMLI_items$testing_results["0","TP"]
  total_false <- GEMLI_items$testing_results["0","FP"]
  
  tpr <- GEMLI_items$testing_results[,"TP"]/total_true
  fpr <- GEMLI_items$testing_results[,"FP"]/total_false
  
  list(fpr = fpr,
       tpr = tpr)
}

tmp <- .compute_tpr_fpr(GEMLI_items_lcl)
fpr_lcl <- tmp$fpr
tpr_lcl <- tmp$tpr
round(cbind(fpr_lcl, tpr_lcl),2)

tmp <- .compute_tpr_fpr(GEMLI_items_pca)
fpr_pca <- tmp$fpr
tpr_pca <- tmp$tpr
round(cbind(fpr_pca, tpr_pca),2)

############

# Compute AUC using the trapezoidal rule
.compute_auc <- function(fpr, tpr){
  auc <- -(sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
  return(2*(auc-0.5))
}

lcl_auc <- .compute_auc(fpr_lcl, tpr_lcl)
pca_auc <- .compute_auc(fpr_pca, tpr_pca)

############

df <- data.frame(
  tpr = c(tpr_lcl, tpr_pca),
  fpr = c(fpr_lcl, fpr_pca),
  model = rep(c("LCL", "log-normalized"),
              each = length(tpr_lcl))
)

library(ggplot2)

# Plot both ROC curves
plot1 <- ggplot(df, aes(x = fpr, y = tpr, color = model)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "red") +  # Diagonal line (random classifier)
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", 
       title = paste0("ROC Curve Comparison",
                      "\nAUC: LCL: ", round(lcl_auc, 2), 
                      ", log: ", round(pca_auc, 2))) +
  scale_color_manual(values = c(rgb(212, 63, 136, maxColorValue = 255), 
                                rgb(117, 164, 58, maxColorValue = 255))) +  # Customize colors if you want
  coord_equal() 

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup7_GEMLI_ROC.png"),
                height = 1200, width = 1200, units = "px")


###

# Plot both ROC curves
plot1 <- ggplot(df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 1.5) +
  geom_abline(linetype = "dashed", color = "red") +  # Diagonal line (random classifier)
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c(rgb(212, 63, 136, maxColorValue = 255), 
                                rgb(117, 164, 58, maxColorValue = 255))) +  # Customize colors if you want
  coord_equal() +
  Seurat::NoLegend()

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup7_GEMLI_ROC_cleaned.png"),
                height = 900, width = 900, units = "px")

########


df <- data.frame(
  tpr = c(tpr_pca),
  fpr = c(fpr_pca),
  model = rep(c("log-normalized"),
              each = length(tpr_lcl))
)


# Plot both ROC curves
plot1 <- ggplot(df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 1.5) +
  geom_abline(linetype = "dashed", color = "red") +  # Diagonal line (random classifier)
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c(rgb(117, 164, 58, maxColorValue = 255))) +  # Customize colors if you want
  coord_equal() +
  Seurat::NoLegend()

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup7_GEMLI_ROC_cleaned_only-pca.png"),
                height = 900, width = 900, units = "px")

