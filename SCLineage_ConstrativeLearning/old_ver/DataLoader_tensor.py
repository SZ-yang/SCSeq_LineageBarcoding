# general package
import tempfile
import os
import scipy
import numpy as np
import pandas as pd
import copy
import random

# single cell package
import scanpy as sc
import anndata as ad

# deep learning package
import torch
import torchvision.models as models
import torchvision.transforms as T
import torch.nn as nn
import torch.nn.functional as F
from torchvision.datasets import STL10
from torch.utils.data import DataLoader
from torch.multiprocessing import cpu_count
import torchvision.transforms as T


'''
left for improvements:
1. Randomness Control for Reproductivity
2. Efficiency
3. Outputdata Type -> tensor
4. Documentation (actually batch size = self.batch_size*2)
5. it won't handle if a lineage only has one cell

'''

class SClineage_DataLoader:
    def __init__(self, count_matrix, lineages, batch_size = 10):
        """
        Args:
            count_matrix (numpy.ndarray): Data array of shape (n, p) with n samples and p features.
            lineage (numpy.ndarray): Array of shape (n, 1) with group labels for each sample.
        
        Vars:
            avail_lineage (dictionary): there's at least a cell in a lineage that has not yet assigned.
        """
        self.count_matrix = count_matrix
        self.lineages = lineages
        self.batch_size = batch_size
        self.batch_all = {}
        self.num_batch = 0 # denotes both as the key for self.batch_all and the total number of batches

        # lineage indices for sampling
        self.lineage_indices = {}
        for idx, group in enumerate(lineages):
            if group[0] not in self.lineage_indices:
                self.lineage_indices[group[0]] = []
            self.lineage_indices[group[0]].append(idx)

        self.avail_lineages = copy.deepcopy(self.lineage_indices)
    def __len__(self):
        return len(self.count_matrix)
    
    def batch_generator(self):
        while len(self.avail_lineages) >= 1:
            if len(self.avail_lineages) >= self.batch_size:

                # uniformly sample m out of available lineages
                sampled_lineages_idx = random.sample(list(self.avail_lineages.keys()), self.batch_size) 
                # sampled_lineages = {key: self.avail_lineages[key] for key in sampled_lineages_idx} 
                single_batch = []
                for key in sampled_lineages_idx:
                    if len(self.avail_lineages[key])>1: 
                        sampled_2cell = random.sample(self.avail_lineages[key],2)
                        single_batch += [tuple(sampled_2cell)]
                        
                        # remove the chosen cells in this lineage
                        self.avail_lineages[key].remove(sampled_2cell[0])
                        self.avail_lineages[key].remove(sampled_2cell[1]) #self.avail_lineages[key] = [item for item in self.avail_lineages[key] if item not in sampled_2cell]
                        # check if this lineage is still available 
                        if len(self.avail_lineages[key]) == 0:
                            del self.avail_lineages[key]
                              
                    else: # only 1 cell in this lineage
                        lefted_cell = self.avail_lineages[key][0]

                        # Randomly sample one cell from the overall lineage dict
                        filtered_cells = [cell for cell in self.lineage_indices[key] if cell != lefted_cell]
                        selected_cell = random.choice(filtered_cells) if filtered_cells else None


                        sampled_2cell = [lefted_cell,selected_cell]
                        single_batch += [tuple(sampled_2cell)]

                        # remove the chosen cells in this lineage and this lineage
                        # self.avail_lineages[key].remove(sampled_2cell[0])
                        del self.avail_lineages[key]
                        
                
                self.num_batch += 1
                
                # Convert the entire batch to a tensor before storing
                # Assuming each element of single_batch is a tuple (tensor1, tensor2)
                tensor_batch = [(torch.from_numpy(self.count_matrix[i]).float(), torch.from_numpy(self.count_matrix[j]).float()) for i, j in single_batch]
                self.batch_all[self.num_batch] = tensor_batch
                
                
            else:
                # generate a set of linages with a full batch size
                num_lgs_to_add = self.batch_size - len(self.avail_lineages)
                #filtered_lgs_ls = [lg for lg in list(self.lineage_indices.keys()) if lg not in list(self.avail_lineages.keys())]
                filtered_lgs_ls = [lg for lg in self.lineage_indices if lg not in self.avail_lineages]
                selected_lgs = random.sample(filtered_lgs_ls, num_lgs_to_add)

                single_batch = []
                
                # for the available lineages (same as above):
                for key in list(self.avail_lineages.keys()):
                    if len(self.avail_lineages[key])>1: 
                        sampled_2cell = random.sample(self.avail_lineages[key],2)
                        single_batch += [tuple(sampled_2cell)]
                        
                        # remove the chosen cells in this lineage
                        self.avail_lineages[key].remove(sampled_2cell[0])
                        self.avail_lineages[key].remove(sampled_2cell[1]) #self.avail_lineages[key] = [item for item in self.avail_lineages[key] if item not in sampled_2cell]
                        # check if this lineage is still available 
                        if len(self.avail_lineages[key]) == 0:
                            del self.avail_lineages[key]
                              
                    else: # only 1 cell in this lineage
                        lefted_cell = self.avail_lineages[key][0]

                        # Randomly sample one cell from the overall lineage dict
                        filtered_cells = [cell for cell in self.lineage_indices[key] if cell != lefted_cell]
                        selected_cell = random.choice(filtered_cells) if filtered_cells else None


                        sampled_2cell = [lefted_cell,selected_cell]
                        single_batch += [tuple(sampled_2cell)]

                        # remove the chosen cells in this lineage and this lineage
                        # self.avail_lineages[key].remove(sampled_2cell[0])
                        del self.avail_lineages[key] 

                # for the unavailable lineages:
                for key in selected_lgs:
                    if len(self.lineage_indices[key])>1: 
                        sampled_2cell = random.sample(self.lineage_indices[key],2)
                        single_batch += [tuple(sampled_2cell)]


                self.num_batch += 1
                # Convert the entire batch to a tensor before storing
                # Assuming each element of single_batch is a tuple (tensor1, tensor2)
                tensor_batch = [(torch.from_numpy(self.count_matrix[i]).float(), torch.from_numpy(self.count_matrix[j]).float()) for i, j in single_batch]
                self.batch_all[self.num_batch] = tensor_batch
                

            
        return self.batch_all, self.num_batch 




if __name__ == "__main__":
    # Dummy data
    n, p, M = 500, 5, 15  # 100 samples, 5 features, 4 groups
    data = np.random.randn(n, p)  # Random data
    lineages = np.random.randint(0, M, size=(n, 1))  # Random group labels

    for i in range(15):
        loader = SClineage_DataLoader(data, lineages, batch_size=10)
        batches, num_batches = loader.batch_generator()

        print(num_batches, loader.avail_lineages)
