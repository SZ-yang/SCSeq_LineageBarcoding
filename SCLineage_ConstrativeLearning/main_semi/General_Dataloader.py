import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import DataLoader_combination_final as dl
import SCDataset as ds


def General_DataLoader(file_path, batch_size, size_factor, unlabeled_per_batch, batch_seed):
    """
    Wrapper to generate contrastive learning batches with labeled and unlabeled subsets.

    Args:
        file_path (str): Path to .h5ad file containing count_matrix and 'clone_id' obs.
        batch_size (int): Number of positive lineage-pairs per batch.
        size_factor (float): Fraction of available lineage-pairs to sample per lineage.
        batch_seed (int): Random seed for reproducibility.

    Returns:
        batch_all_label (dict): batch_idx -> list of (tensor_i, tensor_j) positive pairs.
        batch_all_unlabel (dict): batch_idx -> list of tensor_k unlabeled samples.
        lineage_array_label (np.ndarray): (total_labeled, 1) true labels for labeled cells.
        lineage_array_unlabel (dict): batch_idx -> (n_unlabeled,1) true labels for unlabeled cells.
        num_batches (int): total number of batches.
    """
    # Load AnnData and extract raw counts and lineage labels
    adata = sc.read_h5ad(file_path)
    count_matrix = adata.X
    cell_lineage = adata.obs['clone_id'].values.reshape(-1, 1)

    # Initialize custom DataLoader_combination
    DLoader = dl.SClineage_DataLoader(
        count_matrix,
        cell_lineage,
        batch_size=batch_size,
        size_factor=size_factor,
        unlabeled_per_batch=unlabeled_per_batch, 
        seed=batch_seed
        
    )

    # Unpack labeled vs. unlabeled batches and labels
    batch_all_label, batch_all_unlabel,num_batches,lineage_array_label,lineage_array_unlabel = DLoader.batch_generator()

    print(f"Number of batches: {num_batches}")
    print(f"Total labeled pairs: {num_batches} * {batch_size} = {num_batches * batch_size}")
    # show how many unlabeled per batch (constant across batches)
    sample_unlabel = len(batch_all_unlabel[next(iter(batch_all_unlabel.keys()))])
    print(f"Unlabeled per batch: {sample_unlabel}")

    return batch_all_label, batch_all_unlabel, num_batches, lineage_array_label, lineage_array_unlabel



if __name__ == "__main__":
    # Example usage
    file_path = "/Users/apple/Desktop/KB/data/LarryData/Larry_200.h5ad"
    batch_size = 30
    size_factor = 0.5
    unlabeled_per_batch = 10
    seed = 42

    batches_label, batches_unlabel, n_batches, labels_label, labels_unlabel = \
        General_DataLoader(file_path, batch_size, size_factor, unlabeled_per_batch, seed)

    print(f"Labels_label shape: {labels_label.shape}")
    # inspect first batch
    b0 = 0
    print(f"Batch {b0} has {len(batches_label[b0])} labeled pairs and {len(batches_unlabel[b0])} unlabeled samples.")
    print("Labeled example:", batches_label[b0][0])
    print("Unlabeled example shape:", batches_unlabel[b0][0].shape)
