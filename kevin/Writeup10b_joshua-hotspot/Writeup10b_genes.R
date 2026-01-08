rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)

data_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-hotspot/"
out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup10b_joshua-hotspot/kevin_results/"

df <- read.csv(paste0(data_folder, "LCL_celltagMulti_hotspot_autocorrelations.csv"))
colnames(df)[-1] <- paste0("LCL_", colnames(df)[-1])
df_scvi <- read.csv(paste0(data_folder, "scVI_celltagMulti_hotspot_autocorrelations.csv"))
colnames(df_scvi)[-1] <- paste0("scVI_", colnames(df_scvi)[-1])

df <- merge(x = df, y = df_scvi, by = "Gene")
df <- df[order(df$LCL_Z, decreasing = TRUE),]
rownames(df) <- df$Gene

plot(df$LCL_Pval, df$scVI_Pval, asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
lines(c(0,1), c(0,1), col = 2, lwd = 2, lty = 2)

plot(df$LCL_FDR, df$scVI_FDR, asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
lines(c(0,1), c(0,1), col = 2, lwd = 2, lty = 2)

plot(df$LCL_Z, df$scVI_Z, asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
lines(c(0,1e6), c(0,1e6), col = 2, lwd = 2, lty = 2)

plot(df$LCL_C, df$scVI_C, asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
lines(c(0,1e6), c(0,1e6), col = 2, lwd = 2, lty = 2)


#####################

idx <- intersect(which(df$LCL_FDR <= 0.05),
                 which(df$scVI_FDR >= 0.05))
df[idx,"Gene"]

write.table(df[idx,"Gene"],
            file = paste0(out_folder, "Writeup10b_in-lcl_not-in-scvi.csv"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

####

idx <- intersect(which(df$scVI_FDR <= 0.05),
                 which(df$LCL_FDR >= 0.05))
df[idx,"Gene"]


write.table(df[idx,"Gene"],
            file = paste0(out_folder, "Writeup10b_in-scvi_not-in-lcl.csv"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

#####################

specific_genes <- c("Cdh1", "Apoa1", "Mettl7a1", "Wnt4", "Spint2",
                    "Dlk1", "Peg3", "Foxd2", "Zfp281")
table(specific_genes %in% df$Gene)
df[which(df$Gene %in% specific_genes),]
