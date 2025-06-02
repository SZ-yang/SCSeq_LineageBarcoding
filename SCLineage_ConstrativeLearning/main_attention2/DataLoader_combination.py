import numpy as np
import random
import torch
import scipy.sparse as sp
import scanpy as sc


class SClineage_DataLoader:
    def __init__(self, count_matrix, lineages, 
                 num_lineage_1batch=30, 
                 cells_per_lineage=4, 
                 size_factor=0.5, 
                 seed=None):
        """
        :param count_matrix: Expression matrix (dense or sparse). Rows = cells, columns = features.
        :param lineages: Array of lineage identifiers (same length as number of cells).
        :param num_lineage_1batch: Number of lineages to sample in each batch.
        :param cells_per_lineage: Number of cells to sample per lineage.
        :param size_factor: Controls how many cell-groups we sample per lineage 
                            (relative to the number of cells in that lineage).
        :param seed: Random seed for reproducibility.
        """
        # Convert to dense numpy array if it's sparse
        try:
            self.count_matrix = count_matrix.toarray()
        except AttributeError:
            self.count_matrix = count_matrix
        
        self.lineages = lineages.flatten()
        self.num_lineage_1batch = num_lineage_1batch
        self.cells_per_lineage = cells_per_lineage
        self.size_factor = size_factor

        # Set random seed if specified
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        # Precompute lineage information and sample groups
        self.lineage_info = self.generate_lineage_info()
        self.avail_lineage_groups = self.generate_avail_lineage_groups()

        # Generate all batches (each batch is a list of (lineage_id, cell_group) pairs)
        self.batch_all_index = self.generate_batch_all_index()

        # Convert each batch into a torch tensor
        self.batch_all = self.generate_batch_all()

        # Flatten lineage labels in the same order as batch_all
        self.lineage_array = self.generate_lineage_array()

    def generate_lineage_info(self):
        """
        Create a dictionary mapping each unique lineage ID to 
        a list of row indices (cells) belonging to that lineage.
        """
        lineage_info = {}
        unique_elements = np.unique(self.lineages)
        for element in unique_elements:
            lineage_info[element] = np.where(self.lineages == element)[0].tolist()
        return lineage_info

    def sample_random_groups(self, indices, cells_per_lineage, group_count):
        """
        Randomly sample 'group_count' groups, each of size 'cells_per_lineage', 
        from the given 'indices'. This avoids enumerating all possible combinations.
        """
        groups = []
        # If there aren't enough cells to form a group, skip
        if len(indices) < cells_per_lineage:
            return groups

        for _ in range(group_count):
            # Randomly choose 'cells_per_lineage' cells without replacement
            group = tuple(np.random.choice(indices, size=cells_per_lineage, replace=False))
            groups.append(group)
        return groups

    def generate_avail_lineage_groups(self):
        """
        For each lineage, randomly sample some number of cell-groups (each group
        has 'cells_per_lineage' cells). The number of groups per lineage is 
        determined by 'size_factor * len(indices)' (with a minimum of 1).
        """
        avail_lineage_groups = {}
        for lineage_id, indices in self.lineage_info.items():
            # Decide how many groups to sample based on 'size_factor'
            group_count = max(int(self.size_factor * len(indices)), 1)
            # Randomly sample these groups
            groups = self.sample_random_groups(indices, self.cells_per_lineage, group_count)
            avail_lineage_groups[lineage_id] = groups
        return avail_lineage_groups

    def generate_batch_all_index(self):
        """
        Generate a dictionary of all batches:
            batch_all_index[batch_idx] = [(lineage_id, group), (lineage_id, group), ...]
        Each 'group' is a tuple of cell indices.
        """
        # Make a working copy so we can remove groups as they're sampled
        lineage_groups_cp = {k: v.copy() for k, v in self.avail_lineage_groups.items()}
        batches = {}
        batch_idx = 0

        # Continue until all sampled groups are exhausted
        while lineage_groups_cp:
            available_lineages = list(lineage_groups_cp.keys())
            if len(available_lineages) >= self.num_lineage_1batch:
                selected_lineages = random.sample(available_lineages, self.num_lineage_1batch)
            else:
                # If not enough lineages remain, sample the remainder from lineages that already exist
                remaining = self.num_lineage_1batch - len(available_lineages)
                complement_lineages = list(
                    set(self.avail_lineage_groups.keys()) - set(available_lineages)
                )
                # If complement_lineages is empty, we might randomly pick from the original set
                # to fill up the batch. Here we do a simple approach:
                if not complement_lineages:
                    complement_lineages = available_lineages
                selected_lineages = available_lineages + random.choices(complement_lineages, k=remaining)

            batch = []
            for lineage in selected_lineages:
                # If this lineage still has unsampled groups, pick one at random
                if lineage in lineage_groups_cp and lineage_groups_cp[lineage]:
                    group = random.choice(lineage_groups_cp[lineage])
                    lineage_groups_cp[lineage].remove(group)
                    if not lineage_groups_cp[lineage]:
                        del lineage_groups_cp[lineage]
                else:
                    # Otherwise, just pick from the original (this means repeating some groups)
                    group = random.choice(self.avail_lineage_groups[lineage])
                batch.append((lineage, group))

            batches[batch_idx] = batch
            batch_idx += 1

        return batches

    def generate_batch_all(self):
        """
        Convert each batch of cell indices into a single torch.Tensor of shape:
            (num_lineages_in_batch, cells_per_lineage, num_features)
        """
        batch_tensors = {}
        for idx, batch in self.batch_all_index.items():
            groups_tensor = []
            for lineage, group in batch:
                # Vectorized approach: slice rows from count_matrix, convert to torch
                cells_tensor = torch.from_numpy(self.count_matrix[list(group), :])
                groups_tensor.append(cells_tensor)
            # Stack along the first dimension => shape: (num_lineages_in_batch, cells_per_lineage, num_features)
            batch_tensors[idx] = torch.stack(groups_tensor)
        return batch_tensors

    def generate_lineage_array(self):
        """
        Create a numpy array of shape (total_lineages_sampled, 1), 
        in the same order as the batches in 'batch_all_index'.
        """
        lineage_labels = []
        for batch in self.batch_all_index.values():
            for lineage, _ in batch:
                lineage_labels.append(lineage)
        return np.array(lineage_labels).reshape(-1, 1)

    def batch_generator(self):
        """
        Print some stats about how many cells belong to each lineage
        and return (batch_dictionary, total_batches, lineage_label_array).
        """
        lengths = [len(indices) for indices in self.lineage_info.values()]
        min_length, max_length, mean_length = min(lengths), max(lengths), round(np.mean(lengths), 2)
        print(f"Cells per lineage: range ({min_length}, {max_length}), average: {mean_length}")

        return self.batch_all, len(self.batch_all), self.lineage_array


if __name__ == "__main__":
    do = 1
    if do == 0:
        # Simple test example
        n, p = 11373, 10
        data = sp.csr_matrix(np.random.randn(n, p))
        lineages = np.random.choice(np.arange(1, 200), size=n)

        num_lineage_1batch = 30
        cells_per_lineage = 4
        size_factor = 1
        seed = 42

        generator = SClineage_DataLoader(
            data, 
            lineages, 
            num_lineage_1batch=num_lineage_1batch,
            cells_per_lineage=cells_per_lineage, 
            size_factor=size_factor, 
            seed=seed
        )

        batch_all, num_batches, lineage_array = generator.batch_generator()

        print("Number of batches:", num_batches)
        for batch_idx, batch_tensor in batch_all.items():
            print(f"Batch {batch_idx} shape:", batch_tensor.shape)

        print("Lineage array shape:", lineage_array.shape)
    
    else:
        file_path = "/Users/apple/Desktop/KB/data/LarryData/Larry_200.h5ad"
        # "/Users/apple/Desktop/KB/data/LarryData/Larry_41201_2000.h5ad"
        adata = sc.read_h5ad(file_path)
        count_matrix = adata.X
        print(count_matrix.shape)
        cell_lineage = adata.obs['clone_id'].values.reshape(-1, 1)
        num_lineage_1batch = 150
        cells_per_lineage = 4
        size_factor = 1
        seed = 42

        generator = SClineage_DataLoader(
            count_matrix, 
            cell_lineage, 
            num_lineage_1batch=num_lineage_1batch,
            cells_per_lineage=cells_per_lineage, 
            size_factor=size_factor, 
            seed=seed
        )

        batch_all, num_batches, lineage_array = generator.batch_generator()

        print("num_lineage_1batch:", num_lineage_1batch)
        print("Number of batches:", num_batches)
        for batch_idx, batch_tensor in batch_all.items():
            print(f"Batch {batch_idx} shape:", batch_tensor.shape)

        print("Lineage array shape:", lineage_array.shape)