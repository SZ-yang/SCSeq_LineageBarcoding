rm(list=ls())

library(Seurat)
library(Matrix)

source("~/kzlinlab/projects/scContrastiveLearn/git/SCSeq_LineageBarcoding_kevin/kevin/Writeup10c_celltagmulti-anova/anova_functions.R")
load("~/kzlinlab/data/larry_hematopoiesis/larry-dataset_KZL.RData")
out_folder <- "~/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup10c/"

# for the lineages w/ more than 50 cells, stratify by lineage AND cell-type. 

keep_vec <- !is.na(seurat_obj$assigned_lineage)
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

tab_vec <- table(seurat_obj$assigned_lineage)
lineage_names <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- seurat_obj$assigned_lineage %in% lineage_names
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# for any lineage x cell-type with less than 5 cells, remove those cells
tab_mat <- table(seurat_obj$assigned_lineage, seurat_obj$Cell.type.annotation)

# ============================================================
# 1) Remove cells in any (lineage x cell-type) with < 5 cells
# ============================================================
md <- seurat_obj@meta.data
combo_key <- paste(md$assigned_lineage, md$Cell.type.annotation, sep = "||")
combo_tab <- table(combo_key)

keep_step1 <- combo_tab[combo_key] >= 5
cells_step1 <- rownames(md)[keep_step1]

seurat_obj <- subset(seurat_obj, cells = cells_step1)

# ============================================================
# 2) Remove any lineage whose remaining cells are only 1 cell-type
# ============================================================
md1 <- seurat_obj@meta.data

# number of distinct cell types per lineage (after step 1)
n_types_per_lineage <- tapply(
  md1$Cell.type.annotation,
  md1$assigned_lineage,
  function(x) length(unique(x[!is.na(x)]))
)

keep_lineages <- names(n_types_per_lineage)[n_types_per_lineage > 1]
cells_step2 <- rownames(md1)[md1$assigned_lineage %in% keep_lineages]

seurat_obj <- subset(seurat_obj, cells = cells_step2)

# look at the final meta.data
tab_mat <- table(seurat_obj$assigned_lineage, seurat_obj$Cell.type.annotation)

seurat_obj <- Seurat::DietSeurat(seurat_obj,
                                 features = Seurat::VariableFeatures(seurat_obj))

####################################
####################################
####################################

# Then, compute the mean and variance of each gene's expression for that {lineage x cell-type} group.
# Store this as two "long-form" matrices (columns: {lineage x cell-type}. rows: {gene}), one for the mean, one for the variance

# ----------------------------
# 0) Choose which expression to summarize
#    "data" = normalized/log data (common default)
#    "counts" = raw counts
# ----------------------------

expr <- SeuratObject::LayerData(seurat_obj, assay = "RNA", layer = "data")
md <- seurat_obj@meta.data
lineage_vec <- md$assigned_lineage
celltype_vec <- md$Cell.type.annotation

# ----------------------------
# 1) Build group membership (cells x groups)
# ----------------------------
group_key <- factor(paste(lineage_vec, celltype_vec, sep = "||"))
mm <- Matrix::sparse.model.matrix(~ 0 + group_key)  # cells x groups
colnames(mm) <- sub("^group_key", "", colnames(mm))

n_per_group <- Matrix::colSums(mm)  # group sizes (length = #groups)

# ----------------------------
# 2) Compute means via sums / n
# ----------------------------
sums <- expr %*% mm                                # genes x groups
mean_mat <- sweep(as.matrix(sums), 2, n_per_group, "/")  # genes x groups (dense matrix)

# ----------------------------
# 3) Compute variances (sample variance by default)
#    Var = (sum(x^2) - n*mean^2)/(n-1)
# ----------------------------
expr_sq <- expr
expr_sq@x <- expr_sq@x^2                           # sparse-friendly square
sumsq <- expr_sq %*% mm                            # genes x groups
sumsq_mat <- as.matrix(sumsq)

var_mat <- sweep(
  sumsq_mat - sweep(mean_mat^2, 2, n_per_group, "*"),
  2,
  (n_per_group - 1),
  "/"
)

# tiny negative numerical noise -> clamp to 0
var_mat[var_mat <= 0 & var_mat > -1e-12] <- 1e-6

########################

df <- data.frame(
  combo = colnames(var_mat)
)
df$lineage <- factor(sapply(colnames(var_mat), function(x){strsplit(x, split = "\\|\\|")[[1]][1]}))
df$celltype <- factor(sapply(colnames(var_mat), function(x){strsplit(x, split = "\\|\\|")[[1]][2]}))

# For one gene:
gene_R2 <- matrix(NA, nrow = nrow(expr), ncol = 2)
rownames(gene_R2) <- rownames(expr)
colnames(gene_R2) <- c("Lineage", "Celltype")
for(gene_idx in 1:nrow(expr)){
  if(gene_idx %% floor(nrow(expr)/10) == 0) cat('*')
  
  Was2_dist <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  for(i in 1:(nrow(df)-1)){
    for(j in (i+1):nrow(df)){
      Was2_dist[i,j] <- W2_gauss_1d(mu1 = mean_mat[gene_idx,i], 
                                    var1 = var_mat[gene_idx,i], 
                                    mu2 = mean_mat[gene_idx,j],  
                                    var2 = var_mat[gene_idx,j])
      Was2_dist[j,i] <- Was2_dist[i,j]
    }
  }
  
  res1 <- manova_var_explained(D = Was2_dist, 
                               group = df$lineage)
  res2 <- manova_var_explained(D = Was2_dist, 
                               group = df$celltype)
  
  gene_R2[gene_idx,"Lineage"] <- res1$R2
  gene_R2[gene_idx,"Celltype"] <- res2$R2
}

gene_R2 <- gene_R2[!is.nan(gene_R2[,1]),]
