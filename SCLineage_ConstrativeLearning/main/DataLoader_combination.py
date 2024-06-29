import numpy as np
from itertools import combinations
import random
import copy
import torch
import scanpy as sc
import scipy.sparse as sp

class SClineage_DataLoader:
    def __init__(self, count_matrix, lineages, batch_size=20, size_factor=0.5, seed=None):
        """
        Args:
            count_matrix (scipy.sparse.csr_matrix): Data array of shape (n, p) with n samples and p features.
            lineage (numpy.ndarray): Array of shape (n, 1) with group labels for each sample.
            batch_size(integer): usually takes value 10 or 20 
            size_factor(float): range from 0 to 1.
            seed (int): Random seed for reproducibility.
        """
        try:
        # Attempt to call toarray() method
            self.count_matrix = count_matrix.toarray()
        except AttributeError:
        # Fallback if count_matrix has no toarray() method
            self.count_matrix = count_matrix
        self.lineages = lineages
        self.batch_size = batch_size
        self.size_factor = size_factor
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        
        self.lineage_info = self.generate_lineage_info()
        self.avail_lineage_pairs = self.generate_avail_lineage_pairs()
        self.batch_all_index = self.generate_batch_all_index()
        self.batch_all = self.generate_batch_all()
        self.lineage_array = self.generate_lineage_array()
        
    def generate_lineage_info(self):
        """
        Return a dictionary with key of the unique element of the lineage and the elements a list contains the index.
        """
        lineage_info = {}
        unique_elements = np.unique(self.lineages)
        
        for element in unique_elements:
            lineage_info[element] = np.where(self.lineages == element)[0].tolist()
        
        return lineage_info

    def generate_avail_lineage_pairs(self):
        """
        Return the combinations of the j-choose-2, with the length of size_factor*j*(j-1)/2,
        ensuring that all elements are used at least once.
        """
        avail_lineage_pairs = {}
        
        for key, indices in self.lineage_info.items():
            j = len(indices)
            all_pairs = list(combinations(indices, 2))
            pair_count = int(self.size_factor * (j * (j - 1) / 2))
            
            # Ensure all elements are used at least once
            used_indices = set()
            essential_pairs = []
            
            for pair in all_pairs:
                if pair[0] not in used_indices or pair[1] not in used_indices:
                    essential_pairs.append(pair)
                    used_indices.update(pair)
            
            # Shuffle the remaining pairs and select the required number to fill up to pair_count
            remaining_pairs = [pair for pair in all_pairs if pair not in essential_pairs]
            random.shuffle(remaining_pairs)
            additional_pairs = remaining_pairs[:max(0, pair_count - len(essential_pairs))]
            
            avail_lineage_pairs[key] = essential_pairs + additional_pairs
        
        return avail_lineage_pairs


    def get_min_max_length(self):
        lengths = [len(indices) for indices in self.lineage_info.values()]
        min_length = min(lengths)
        max_length = max(lengths)
        mean_length = round(sum(lengths) / len(lengths), 2)
        
        return min_length, max_length, mean_length

    def generate_batch_all_index(self):
        avail_lineage_pairs_cp = copy.deepcopy(self.avail_lineage_pairs)
        batch_all_index = {}
        i = 0
        
        while len(avail_lineage_pairs_cp.keys()) != 0:
            batch_all_index[i] = []
            if len(avail_lineage_pairs_cp.keys()) >= self.batch_size:
                selected_keys = random.sample(list(avail_lineage_pairs_cp.keys()), self.batch_size)
            else:
                selected_keys = list(avail_lineage_pairs_cp.keys())
            
            for key in selected_keys:
                selected_tuple = random.choice(avail_lineage_pairs_cp[key])
                batch_all_index[i].append(selected_tuple)
                avail_lineage_pairs_cp[key].remove(selected_tuple)
                if not avail_lineage_pairs_cp[key]:
                    del avail_lineage_pairs_cp[key]
            
            if len(selected_keys) < self.batch_size:
                complement_keys = list(set(self.avail_lineage_pairs.keys()) - set(avail_lineage_pairs_cp.keys()))
                remaining_keys = random.sample(complement_keys, self.batch_size - len(selected_keys))
                for key in remaining_keys:
                    selected_tuple = random.choice(self.avail_lineage_pairs[key])
                    batch_all_index[i].append(selected_tuple)
            
            i += 1
        
        return batch_all_index

    def generate_batch_all(self):
        batch_all = {}
        for key, pairs in self.batch_all_index.items():
            batch_all[key] = [(torch.tensor(self.count_matrix[pair[0]]), torch.tensor(self.count_matrix[pair[1]])) for pair in pairs]
        return batch_all

    def generate_lineage_array(self):
        batch_size = len(next(iter(self.batch_all_index.values())))
        m = len(self.batch_all_index) * batch_size
        lineage_array = np.zeros((m, 1), dtype=int)  # Specify dtype as int
        index = 0
        for pairs in self.batch_all_index.values():
            for pair in pairs:
                lineage_array[index] = self.lineages[pair[0]]
                index += 1
        return lineage_array

    def batch_generator(self):
        min_length, max_length, mean_length = self.get_min_max_length()
        print(f"The range of number of cells in a lineage: {min_length, max_length}, average of number of cells in a lineage {mean_length}")
        num_batch = len(self.batch_all.keys())

        return self.batch_all, num_batch, self.lineage_array
        
if __name__ == "__main__":
    do_ = 1
    if do_ == 0:
        n, p= 22, 5, 15 
        data = np.random.randn(n, p)  # Random data
        data = sp.csr_matrix(data)
        lineages = np.array([1, 2, 1, 2, 3, 1, 3, 2,1, 2, 1, 2,1, 2, 1, 2,3, 1, 3,3, 1, 3])
        size_factor = 0.5
        batch_size = 2
        
        generator = SClineage_DataLoader(data, lineages, batch_size, size_factor, seed=42)
        batch_all, num_batch, lineage_info = generator.batch_generator()
        # batch_all = generator.generate_batch_all()
        # generator.get_min_max_length()
        print("num_batch: ", num_batch)
        print("---------------------")
        print(batch_all[1])

    else:
        file_path = "/Users/apple/Desktop/KB/data/LarryData/Larry_41201_2000.h5ad"
        adata = sc.read_h5ad(file_path)
        count_matrix = adata.X
        cell_lineage = adata.obs['clone_id'].values.reshape(-1, 1)
        size_factor = 0.5
        batch_size = 30
        generator = SClineage_DataLoader(count_matrix, cell_lineage, batch_size, size_factor, seed=42)
        batch_all, num_batch, lineage_info = generator.batch_generator()
        print("num_batch: ", num_batch)