import anndata as ad
import numpy as np
import scipy
import pandas as pd
import scanpy as sc
import DataLoader_combination as dl


def General_DataLoader(file_path, num_lineage_1batch , cells_per_lineage, size_factor, batch_seed):

#----------------------------------------------------Generate Batches--------------------------------------------------------
    # input data
    adata = sc.read_h5ad(file_path)
    count_matrix = adata.X
    cell_lineage = adata.obs['clone_id'].values.reshape(-1, 1)
    
    DLoader = dl.SClineage_DataLoader(count_matrix, cell_lineage, num_lineage_1batch, cells_per_lineage, size_factor, seed=batch_seed )
    batch_all, num_batch, lineage_info = DLoader.batch_generator()
    
    print("number of batches: ", num_batch)
    print("number of lineages in one batch:", num_lineage_1batch)

    return batch_all, lineage_info, num_batch

if __name__ == "__main__":
    file_path = "/Users/apple/Desktop/KB/data/LarryData/Larry_200.h5ad"
    # file_path = "/home/users/syang71/Dataset/Larry_41201_2000.h5ad"
    
    batch_size = 150
    size_factor = 1
    batch_all, lineage_info, num_batch = General_DataLoader(file_path,batch_size, size_factor, 42)
    print(lineage_info.shape)
    


