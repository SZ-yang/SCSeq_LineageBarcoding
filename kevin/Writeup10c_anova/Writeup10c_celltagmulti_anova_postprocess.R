rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)

load("~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/Writeup10c_celltagmulti_anova.RData")

anova_mat <- anova_mat[,c("predicted.id_cca_co", "sample", "replicate", "assigned_lineage")]

rowsum_vec <- Matrix::rowSums(anova_mat)
if(max(rowsum_vec) > 1) anova_mat <- anova_mat/(max(rowsum_vec)+0.05)

## Stacked barplot of variance explained (top 50 genes)
## Assumes anova_mat is a gene x covariate numeric matrix/data.frame with rownames = genes

# 1) Pick top 50 genes by total explained variation
top_n <- 50

df_top <- as.data.frame(anova_mat) %>%
  rownames_to_column("gene") %>%
  mutate(total_explained = rowSums(across(-gene), na.rm = TRUE)) %>%
  arrange(desc(total_explained)) %>%
  slice_head(n = top_n)

# (Optional sanity check; keep plot y in [0,1] as requested)
# If some totals exceed 1, coord_cartesian below will still show 0..1, but will clip.
# max(df_top$total_explained, na.rm = TRUE)

# 2) Order covariates: Study_Designation first; others by total contribution among top 50
covars <- setdiff(colnames(anova_mat), "assigned_lineage")
other_order <- df_top %>%
  summarise(across(all_of(covars), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "covariate", values_to = "tot") %>%
  arrange(desc(tot)) %>%
  pull(covariate)

covariate_levels <- rev(c("assigned_lineage", other_order))

# 3) Long format for ggplot
plot_df <- df_top %>%
  select(-total_explained) %>%
  pivot_longer(-gene, names_to = "covariate", values_to = "var_explained") %>%
  mutate(
    gene = factor(gene, levels = df_top$gene),                  # left-to-right by total explained
    covariate = factor(covariate, levels = covariate_levels)    # stacking order
  )

# 4) Colors (Study_Designation fixed to purple; others distinct & colorblind-friendly)
fill_cols <- c(
  "assigned_lineage"                  = "#6a0dad",  # purple (focus)
  "sample"                   = "#009E73",
  "replicate" = "#56B4E9",
  "predicted.id_cca_co"            = "#E69F00"
)

# If your column names differ in capitalization/spaces, you can inspect:
# setdiff(levels(plot_df$covariate), names(fill_cols))

y_max <- max(df_top$total_explained)

plot1 <- ggplot(plot_df, aes(x = gene, y = var_explained, fill = covariate)) +
  geom_col(width = 0.9, color = NA) +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(x = NULL, y = "Variance explained") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.ticks.x = element_blank()
  )

fig_folder <- "~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/fig/kevin/Writeup10c/"

ggplot2::ggsave(plot1,
                filename = paste0(fig_folder, "Writeup10c_celltagmulti_variance.png"),
                width = 6, height = 6)



