import torch
from torch.utils.data import Dataset

# class SCDataset(Dataset):
#     def __init__(self, batches):
#         """
#         Initialize with preprocessed batches.
        
#         Args:
#             batches (dict): A dictionary with batch indices as keys and tensors as values.
#                             Each tensor shape: (num_lineage_1batch, cells_per_lineage, num_features)
#         """
#         self.batch_keys = sorted(batches.keys())
#         self.batches = batches

#     def __len__(self):
#         """Returns the total number of batches."""
#         return len(self.batch_keys)

#     def __getitem__(self, idx):
#         """Returns a batch tensor given an index."""
#         batch_idx = self.batch_keys[idx]
#         batch_tensor = self.batches[batch_idx]  # tensor shape (num_lineage_1batch, cells_per_lineage, num_features)
#         return batch_tensor
    





# old version
# deep learning package
# import torch

# from torch.utils.data import Dataset
# from torch.multiprocessing import cpu_count
# import torchvision.transforms as T


# class SCDataset(Dataset):
#     """
#     A PyTorch Dataset class that wraps around the batches generated by SClineage_DataLoader.
#     This class makes the batches generated by SClineage_DataLoader accessible to PyTorch's DataLoader class.
#     """
#     def __init__(self, batches):
#         """
#         Initializes the dataset with the batches.

#         Args:
#             batches (dict): first output of the 
#             A dictionary where keys are batch indices and values are the batches (list of tuples of tensors).
#         """
#         self.batches = batches
#         self.all_batches = self._flatten_batches(batches)

#     def _flatten_batches(self, batches):
#         """
#         Flattens the batch dictionary into a list of samples.
#         Each sample is a tuple consisting of two tensors.

#         Args:
#             batches (dict): The batches as generated by SClineage_DataLoader.

#         Returns:
#             List of tuples: A flattened list where each item is a tuple of tensors.
#         """
#         flattened = []
#         for batch in batches.values():
#             for sample in batch:
#                 flattened.append(sample)
#         return flattened

#     def __len__(self):
#         """
#         Returns the total number of samples in the dataset.
#         """
#         return len(self.all_batches)

#     def __getitem__(self, idx):
#         """
#         Retrieves the nth sample from the dataset.

#         Args:
#             idx (int): The index of the sample to retrieve.

#         Returns:
#             tuple: A tuple containing two tensors.
#         """
#         # Assuming each item in all_batches is a tuple of tensors
#         return self.all_batches[idx]



import torch
from torch.utils.data import Dataset

class SCDataset(Dataset):
    def __init__(self, batches):
        """
        Custom dataset to wrap batches of tensors.

        Args:
            batches (dict): Dictionary of shape {batch_index: torch.Tensor}
                            Each value should be of shape (cells_per_lineage, input_dim)
        """
        self.batches = batches
        self.keys = list(batches.keys())  # Track order of batches

    def __len__(self):
        return len(self.batches)

    def __getitem__(self, idx):
        """
        Returns a single batch.

        Output:
            Tensor of shape: (cells_per_lineage, input_dim)
            PyTorch's DataLoader will automatically add a batch dimension, so:
            Final shape when passed to the model: (batch_size, cells_per_lineage, input_dim)
        """
        key = self.keys[idx]
        batch_tensor = self.batches[key]  # shape: (cells_per_lineage, input_dim)

        # Optional debug print (comment out if not needed)
        # print(f"[DEBUG] __getitem__ key {key} -> batch shape: {batch_tensor.shape}")


        return batch_tensor