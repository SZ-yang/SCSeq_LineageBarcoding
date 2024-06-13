# scVI environment

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import scanpy as sc
import scvi
import seaborn as sns
import torch
import pandas as pd
import numpy as np
import gc
from collections import Counter
import scipy
from scipy.sparse import csr_matrix, diags

# Set global settings
scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

print("Last run with scvi-tools version:", scvi.__version__)

####################

normed_counts = "/home/users/kzlin/kzlinlab/data/larry_hematopoiesis/stateFate_inVitro_normed_counts.mtx.gz"  
gene_names = "/home/users/kzlin/kzlinlab/data/larry_hematopoiesis/stateFate_inVitro_gene_names.txt.gz" 
clone_matrix = "/home/users/kzlin/kzlinlab/data/larry_hematopoiesis/stateFate_inVitro_clone_matrix.mtx.gz" 
metadata = "/home/users/kzlin/kzlinlab/data/larry_hematopoiesis/stateFate_inVitro_metadata.txt.gz" 

# load data
normed_counts_mat = scipy.io.mmread(normed_counts).tocsr()
genes = pd.read_csv(gene_names, sep='\t',header=None).to_numpy().flatten()
clone_mat = scipy.io.mmread(clone_matrix).tocsr()
meta_df = pd.read_csv(metadata, sep='\t')

# create full adata
adata = ad.AnnData(normed_counts_mat, obs=meta_df, var=pd.DataFrame(index=genes), dtype=np.float32)
# optimize dtypes
adata.obs['Library'] = adata.obs['Library'].astype('category')
adata.obs['Time point'] = adata.obs['Time point'].astype(int)
adata.obs['Starting population'] = adata.obs['Starting population'].astype('category')
adata.obs['Cell type annotation'] = adata.obs['Cell type annotation'].astype('category')
adata.obs['Well'] = adata.obs['Well'].astype(int)
# assign clone_id
adata.obs['clone_id'] = (clone_mat @ np.arange(1,1+clone_mat.shape[1])) - 1
print("number of lineages: ", len(adata.obs['clone_id'].unique()))

# adding raw counts for referring to it in the future
adata.layers["counts"] = adata.X.copy()

# get 2000 genes from the 130887(all) cells
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,
                            n_top_genes=2000,
                            subset=True)

# remove the cells that belong to the clone_id with few cells
value_counts = adata.obs['clone_id'].value_counts()
frequency_dict = {}
for value, count in value_counts.items():
    if count in frequency_dict:
        frequency_dict[count].append(value)
    else:
        frequency_dict[count] = [value]

clone_for_remove = frequency_dict[81585]+frequency_dict[2]+frequency_dict[3]+frequency_dict[4] 
adata = adata[~adata.obs['clone_id'].isin(clone_for_remove)]
print("adata.obs.shape:", adata.obs.shape)

# trim the number of genes to the required number using the highly variable gene calculated from the original adata (all cells)
hvgene = (adata.var.highly_variable[adata.var.highly_variable==True]).index
print("number of the highly variable genes:", len(hvgene))
adata = adata[:,hvgene]
print("adata.X.shape:", adata.X.shape)


# modify the counts to be actual counts
# Compute the row sums
row_sums = np.array(adata.layers["counts"].sum(axis=1)).flatten()
# Handle division by zero for rows with sum 0
row_sums[row_sums == 0] = 1
# Create a diagonal matrix with the inverse of the row sums
row_inv_sums = 1.0 / row_sums
D_inv = diags(row_inv_sums)
# Normalize each row by its row sum
normalized_counts = D_inv.dot(adata.layers["counts"])
# Optionally, you can update the counts in adata
adata.layers["counts"] = normalized_counts

# non_zero_values = adata.layers["counts"].data
# # Compute quantiles
# quantiles = np.quantile(non_zero_values, [0, 0.25, 0.5, 0.7, 1])

scaled_counts = normalized_counts * 1e4
rounded_counts = scaled_counts.copy()
rounded_counts.data = np.round(scaled_counts.data)

# Optionally, you can update the counts in adata
adata.layers["counts"] = rounded_counts

# non_zero_values = adata.layers["counts"].data
# # Compute quantiles
# quantiles = np.quantile(non_zero_values, [0, 0.25, 0.5, 0.75, 1])

# Compute the column sums
column_sums = np.array(adata.layers["counts"].sum(axis=0)).flatten()
# Identify columns with sums greater than 1
columns_to_keep = column_sums > 1
# true_false_counts = np.bincount(columns_to_keep.astype(int))
# Subset the matrix to keep only columns with sums > 1
adata = adata[:, columns_to_keep]

adata = adata.copy()
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts"
)

model = scvi.model.SCVI(adata)
print("Starting to train")
model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

model.save(dir_path="/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup2/Writeup2_LARRY-scvi", 
           overwrite=True)

reserved_names = {'_index'}
if adata.raw is not None:
    raw_var_reserved = reserved_names.intersection(adata.raw.var.columns)
    print("Reserved names in raw.var:", raw_var_reserved)
    if '_index' in adata.raw.var.columns:
        adata.raw.var.rename(columns={'_index': 'index_raw_var'}, inplace=True)

adata.write("/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup2/Writeup2_LARRY-scvi_anndata.h5ad")

print("Done! :)")