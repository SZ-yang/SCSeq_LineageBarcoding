import numpy as np
from itertools import combinations
import random
import copy
import torch
import scanpy as sc
import scipy.sparse as sp

class SClineage_DataLoader:
    def __init__(
        self,
        count_matrix,
        lineages,
        batch_size=20,
        size_factor=0.5,
        unlabeled_per_batch=5,
        seed=None
    ):
        """
        Args:
            count_matrix (scipy.sparse.csr_matrix | np.ndarray): Expression matrix (n_cells × n_genes).
            lineages (np.ndarray): Array of shape (n_cells,1) with lineage labels.
            batch_size (int): Number of positive pairs per batch.
            size_factor (float): Fraction of all possible pairs per lineage.
            unlabeled_per_batch (int): Number of unlabeled cells to sample per batch.
            seed (int): Random seed for reproducibility.
        """
        # Load count_matrix
        try:
            self.count_matrix = count_matrix.toarray()
        except AttributeError:
            self.count_matrix = count_matrix
        self.lineages = lineages.flatten()
        self.batch_size = batch_size
        self.size_factor = size_factor
        self.unlabeled_per_batch = unlabeled_per_batch
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        # Prepare lineage-pair information
        self.lineage_info = self.generate_lineage_info()
        self.avail_lineage_pairs = self.generate_avail_lineage_pairs()
        self.batch_all_index = self.generate_batch_all_index()

        # Generate labeled batches and labels
        self.batch_all_label = self.generate_batch_all()
        self.lineage_array_label = self.generate_lineage_array()

        # Generate unlabeled batches and labels
        self.batch_all_unlabel, self.unlabeled_indices = self.generate_batch_all_unlabel()
        self.lineage_array_unlabel = self.generate_lineage_array_unlabel()

    def generate_lineage_info(self):
        """
        Map each unique lineage to the list of cell indices.
        """
        info = {}
        for lab in np.unique(self.lineages):
            info[lab] = np.where(self.lineages == lab)[0].tolist()
        return info

    def generate_avail_lineage_pairs(self):
        """
        For each lineage, sample ⌈size_factor * (j choose 2)⌉ positive pairs,
        ensuring each cell appears in at least one pair.
        """
        avail = {}
        for lab, inds in self.lineage_info.items():
            all_pairs = list(combinations(inds, 2))
            k = int(self.size_factor * len(all_pairs))
            # ensure coverage
            used, essential = set(), []
            for i, j in all_pairs:
                if i not in used or j not in used:
                    essential.append((i, j))
                    used.update([i, j])
            # fill to k
            remaining = [p for p in all_pairs if p not in essential]
            random.shuffle(remaining)
            extra = remaining[:max(0, k - len(essential))]
            avail[lab] = essential + extra
        return avail

    def generate_batch_all_index(self):
        """
        Partition all available lineage-pairs into batches of size batch_size.
        Returns dict: batch_idx -> list of (i,j) tuples.
        """
        pool = copy.deepcopy(self.avail_lineage_pairs)
        batches, idx = {}, 0
        while pool:
            # select up to batch_size distinct lineages
            labs = list(pool.keys())
            selected = random.sample(labs, min(self.batch_size, len(labs)))
            batches[idx] = []
            for lab in selected:
                pair = random.choice(pool[lab])
                batches[idx].append(pair)
                pool[lab].remove(pair)
                if not pool[lab]:
                    del pool[lab]
            # if fewer, complement from all
            if len(selected) < self.batch_size:
                complement = list(self.avail_lineage_pairs.keys())
                extra = random.sample(complement, self.batch_size - len(selected))
                for lab in extra:
                    batches[idx].append(random.choice(self.avail_lineage_pairs[lab]))
            idx += 1
        return batches

    def generate_batch_all(self):
        """
        Build labeled batches: dict batch_idx -> list of (tensor_i, tensor_j).
        """
        batch_all = {}
        for bidx, pairs in self.batch_all_index.items():
            batch_all[bidx] = [
                (
                    torch.tensor(self.count_matrix[i], dtype=torch.float),
                    torch.tensor(self.count_matrix[j], dtype=torch.float)
                ) for (i, j) in pairs
            ]
        return batch_all

    def generate_lineage_array(self):
        """
        Build lineage labels for labeled data: array of shape (2*batch_size_total,1).
        """
        labels = []
        for pairs in self.batch_all_index.values():
            for i, j in pairs:
                labels.append(self.lineages[i])
                labels.append(self.lineages[j])
        return np.array(labels, dtype=int).reshape(-1, 1)

    def generate_batch_all_unlabel(self):
        """
        Build unlabeled batches: dict batch_idx -> list of tensor_unlabeled,
        and record their indices for label retrieval.
        """
        batch_unlab, indices_unlab = {}, {}
        n = len(self.lineages)
        for bidx, pairs in self.batch_all_index.items():
            # get used lineages this batch
            used = {self.lineages[i] for i, _ in pairs}
            pool = [idx for idx in range(n) if self.lineages[idx] not in used]
            sample_n = min(self.unlabeled_per_batch, len(pool))
            chosen = random.sample(pool, sample_n)
            batch_unlab[bidx] = [torch.tensor(self.count_matrix[idx], dtype=torch.float) for idx in chosen]
            indices_unlab[bidx] = chosen
        return batch_unlab, indices_unlab

    def generate_lineage_array_unlabel(self):
        """
        Build lineage labels for unlabeled data: dict batch_idx -> ndarray (sample_n,1).
        """
        labels_unlab = {}
        for bidx, inds in self.unlabeled_indices.items():
            labels_unlab[bidx] = np.array([self.lineages[idx] for idx in inds], dtype=int).reshape(-1, 1)
        return labels_unlab

    def batch_generator(self):
        """
        Returns:
          batch_all_label, batch_all_unlabel,
          num_batches,
          lineage_array_label, lineage_array_unlabel
        """
        num_batches = len(self.batch_all_label)
        return (
            self.batch_all_label,
            self.batch_all_unlabel,
            num_batches,
            self.lineage_array_label,
            self.lineage_array_unlabel
        )

if __name__ == "__main__":
    # Example usage
    adata = sc.read_h5ad("/Users/apple/Desktop/KB/data/LarryData/Larry_200.h5ad")
    matrix = adata.X
    lineage = adata.obs['clone_id'].values.reshape(-1, 1)
    loader = SClineage_DataLoader(
        matrix,
        lineage,
        batch_size=30,
        size_factor=0.5,
        unlabeled_per_batch=10,
        seed=42
    )
    batches_label, batches_unlabel, n_batches, labels_label, labels_unlabel = loader.batch_generator()
    print(f"Batches: {n_batches}")
    for b in range(min(1, n_batches)):
        print(f"Batch {b}: {len(batches_label[b])} pairs, {len(batches_unlabel[b])} unlabeled")
        print("Label shapes:", labels_label[b].shape, labels_unlabel[b].shape)
        print("batches_label[b]: ", batches_label[b])
