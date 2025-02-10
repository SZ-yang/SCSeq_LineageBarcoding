rm(list=ls())
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/out/kevin/Writeup7_joshua-results_dylan/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Joshua/git/SCSeq_LineageBarcoding/fig/kevin/Writeup7/"

load(paste0(out_folder, "adata_with_lcl.RData"))

color_vec <- c("cis" = rgb(202,111,55, maxColorValue = 255),
               "cistocis" = rgb(248,210,152, maxColorValue = 255),
               "cistococl2" = rgb(240,148,71, maxColorValue = 255),
               "cistodabtram" = rgb(160,61,38, maxColorValue = 255),
               "cocl2" = rgb(69,132,69, maxColorValue = 255),
               "cocl2tocis" = rgb(131,202,163, maxColorValue = 255),
               "cocl2tococl2" = rgb(126,191,90, maxColorValue = 255),
               "cocl2todabtram" = rgb(35,63,58, maxColorValue = 255),
               "dabtram" = rgb(68,49,147, maxColorValue = 255),
               "dabtramtocis" = rgb(147,137,193, maxColorValue = 255),
               "dabtramtococl2" = rgb(145,54,147, maxColorValue = 255),
               "dabtramtodabtram" = rgb(68,32,85, maxColorValue = 255))

tab_mat <- table(seurat_obj$Lineage, seurat_obj$OG_condition)

treatment_vec <- c("cis", "cocl2", "dabtram")
for(treatment in treatment_vec){
  col_idx <- grep(paste0("^", treatment), colnames(tab_mat))
  tab_mat2 <- tab_mat[,col_idx]
  tab_mat2 <- tab_mat2[order(tab_mat2[,treatment], decreasing = TRUE),]
  tab_mat2 <- log10(tab_mat2+1)
  tab_mat2 <- tab_mat2[1:20,]
  
  df <- as.data.frame(tab_mat2)
  
  # Create an index for the x-axis, and keep gene names
  df$index  <- rep(seq_len(nrow(tab_mat2)), times = ncol(tab_mat2))
  
  # We also need a min and max per gene (per index) to make the vertical line
  df_segments <- df %>%
    group_by(index) %>%
    summarize(ymin = min(Freq), ymax = max(Freq), .groups = "drop")
  
  plot1 <- ggplot() +
    # 1) Vertical line connecting the three points for each gene
    geom_segment(
      data = df_segments,
      aes(
        x     = index,
        xend  = index,
        y     = ymin,
        yend  = ymax
      ),
      color = "gray50"
    ) +
    # 2) The points themselves (colored by which column they came from)
    geom_point(
      data = df,
      aes(
        x     = index,
        y     = Freq,
        color = Var2,
        size  = ifelse(Var2 %in% c("cis", "cocl2", "dabtram"), 4, 2)  # Larger size for "cis"
      )
    ) +
    # Optional: set your own colors here
    scale_color_manual(values = color_vec) +
    scale_size_identity() +  # Ensures the manual size values are used directly
    # Style the axes, etc.
    labs(
      x = "Lineage (in sorted order)",
      y = "Number of cells (Log10 scale)"
    ) +
    theme_bw(base_size = 14)
  
  plot1 <- plot1 + Seurat::NoLegend()
  
  ggplot2::ggsave(plot1, 
                  filename = paste0(plot_folder, "Writeup7_lineage-size_barplot_", treatment, ".png"),
                  height = 4, width = 4)
}
