rm(list=ls())

library(Seurat)
library(variancePartition)
library(ggplot2)
library(ggrepel)
library(dplyr)

data_folder <- "~/kzlinlab/data/shaffer_clonal-treatment/"
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/"
fig_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup7/"

load(paste0(out_folder, "Writeup7_shaffer_variancePartition.RData"))

plot1 <- variancePartition::plotVarPart(varPart)

ggsave(plot1, height = 5, width = 5,
       filename = paste0(fig_folder, "Writeup7_shaffer_variancePartition.png"))

###

df <- varPart

variable_vec <- c("OG_condition", "G2M.Score")
for(variable in variable_vec){
  # Compute quantiles
  cutoff_reprog <- sort(df$Lineage, decreasing = TRUE)[50]
  cutoff_celltag <- sort(df[,variable], decreasing = TRUE)[50]
  
  # Identify genes to label
  df$label <- ifelse(df$Lineage >= cutoff_reprog |
                       df[,variable] >= cutoff_celltag, rownames(df), "")
  
  # Make scatterplot
  plot1 <- eval(parse(text = paste0("ggplot(df, aes(x = ", variable, ", y = Lineage))")))
  plot1 <- plot1 + geom_point() +
    geom_text_repel(aes(label = label), color = "red", max.overlaps = 10) +
    labs(x = variable, y = "Lineage", title = paste0("Scatterplot of Lineage vs ", variable))
  ggsave(plot1, height = 8, width = 8,
         filename = paste0(fig_folder, "Writeup7_vp_", variable, "-vs-Lineage.png"))
  
}

