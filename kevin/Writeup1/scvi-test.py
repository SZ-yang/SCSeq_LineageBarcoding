import os
import tempfile
import scanpy as sc
import scvi
import seaborn as sns
import torch

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

adata = scvi.data.heart_cell_atlas_subsampled(save_path=save_dir.name)
adata

sc.pp.filter_genes(adata, min_counts=3)

adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="cell_source",
)

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["cell_source", "donor"],
    continuous_covariate_keys=["percent_mito", "percent_ribo"],
)

model = scvi.model.SCVI(adata)

model

model.train()

model_dir = os.path.join("/home/users/kzlin/kzlinlab/projects/lineageBarcodingVAE/out/kevin/Writeup1/scvi_model")
model.save(model_dir, overwrite=True)

print("Done! :)")