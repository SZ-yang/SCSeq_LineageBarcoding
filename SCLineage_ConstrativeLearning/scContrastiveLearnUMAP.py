import numpy as np
import umap
import matplotlib.pyplot as plt

# Load your dataset
data = np.load('/home/users/syang71/kzlinlab/projects/lineageBarcodingCL/git/scContrastiveLearn_Joshua/scBaseEncoderFeat.npy')  
# Optionally, if you have labels # Comment this line if you don't have labels
labels = np.load("/home/users/syang71/kzlinlab/projects/lineageBarcodingCL/git/scContrastiveLearn_Joshua/scBaseEncoderFeat.npy")

# Initialize UMAP
reducer = umap.UMAP()

# Fit the model to your data
embedding = reducer.fit_transform(data)

# Plotting
plt.figure(figsize=(10, 7))
# If you have labels for coloring
plt.scatter(embedding[:, 0], embedding[:, 1], c=labels, cmap='Spectral', s=5)
# If you don't have labels, just plot the embedding
# plt.scatter(embedding[:, 0], embedding[:, 1], s=5)
plt.colorbar()  # Comment this line if you don't have labels
plt.title('UMAP projection')

plt.savefig('/home/users/syang71/kzlinlab/projects/lineageBarcodingCL/git/scContrastiveLearn_Joshua/umap_plot.png', dpi=300) 
plt.savefig('/home/users/syang71/kzlinlab/projects/lineageBarcodingCL/git/scContrastiveLearn_Joshua/umap_plot.svg', format='svg')

plt.show()

