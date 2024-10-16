# from https://cospar.readthedocs.io/en/latest/20210121_all_hematopoietic_data_v3.html
import cospar as cs
import numpy as np
import anndata as ad

cs.logging.print_version()
cs.settings.verbosity = 2
cs.settings.data_path = "/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup4/LARRY_data"  # A relative path to save data. If not existed before, create a new one.
cs.settings.figure_path = "/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup4/LARRY_figure"  # A relative path to save figures. If not existed before, create a new one.

adata_orig = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/Larry_41093_2000_norm_log_cleaned.h5ad")
print(adata_orig)

# Convert the 'time_info' (or 'Time point') column to string type
adata_orig.obs['time_info'] = adata_orig.obs['time_info'].astype(int).astype(str)

cs.hf.check_available_choices(adata_orig)

adata = cs.tmap.infer_Tmap_from_multitime_clones(
    adata_orig,
    clonal_time_points=["2", "4", "6"],
    later_time_point="6",
    smooth_array=[20, 15, 10],
    sparsity_threshold=0.2,
    max_iter_N=3,
)

adata.write("/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup4/Writeup4_cospar-LARRY.h5ad")