# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html

import os
import tempfile
import scanpy as sc
import scvi
import seaborn as sns
import torch
import anndata

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)
sc.settings.figdir = "/home/users/kzlin/kzlinlab/projects/lineageBarcodingVAE/git/lineageEmbedding_kevin/fig/kevin/Writeup1/"

adata = anndata.read_h5ad("/home/users/kzlin/kzlinlab/projects/lineageBarcodingVAE/out/kevin/Writeup1/scvi_adata.h5ad")
adata

model = scvi.model.SCVI.load("/home/users/kzlin/kzlinlab/projects/lineageBarcodingVAE/out/kevin/Writeup1/scvi_model", 
                            adata=adata)
model

SCVI_LATENT_KEY = "X_scVI"
latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

SCVI_NORMALIZED_KEY = "scvi_normalized"
adata.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)

# run PCA then generate UMAP plots
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)

filename = "scvi-test_cell-type.png"
plot_save_path_umap = plot_save_path + filename
sc.pl.umap(
    adata,
    color=["cell_type"],
    frameon=False,
    save=filename
)

filename = "scvi-test_donor_cell-source.png"
sc.pl.umap(
    adata,
    color=["donor", "cell_source"],
    ncols=2,
    frameon=False,
    save=filename
)