rm(list=ls())

library(Seurat)
library(GEMLI)

out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/"
load(paste0(out_folder, "Writeup5_Larry_log_GEMLI.RData"))

zz <- GEMLI_items[['prediction']] # we need to look at how they did it on LARRY...
quantile(zz[zz!=0])