from sklearn.metrics import classification_report, confusion_matrix
import numpy as np
import anndata as ad
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score

# Load data
# input_dir = "/home/users/syang71/kzlinlab/projects/scContrastiveLearn/out/joshua/data/feat_414/feat_414__bs25_temp0.5"
# X = np.load(input_dir+"/scBaseEncoderFeat_Z_bs25_tau0.5.npy")
input_dir = "/home/users/syang71/Dataset"
X = np.load(input_dir+"/Larry_41201_2000_pca.npy")
adata_subset = ad.read_h5ad('/home/users/syang71/Dataset/Larry_41201_2000.h5ad')
y = adata_subset.obs["clone_id"].to_numpy()

# Initialize the KNN classifier with 31 neighbors
# We use 31 because it includes the point itself
knn = KNeighborsClassifier(n_neighbors=31)

# Fit the classifier on the entire dataset
knn.fit(X, y)

# Finding the 31 nearest neighbors for each point in X
neighbors = knn.kneighbors(X, return_distance=False)

# Predict labels by excluding the point itself and taking the mode of the labels of the 30 nearest neighbors
y_pred = np.array([np.argmax(np.bincount(y[neighbors[i][neighbors[i] != i]])) for i in range(len(X))])

# Calculate the accuracy
accuracy = accuracy_score(y, y_pred)
print("Accuracy:", accuracy)
