import numpy as np
import umap
import matplotlib.pyplot as plt
import anndata as ad
from sklearn.metrics import calinski_harabasz_score, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split

class Eval:
    def __init__(self, embedding, anndata, num_top_clone_ids=5, custom_palette=None):
        """
        Initialize UMAPAnalyzer with embedding and anndata.
        
        Parameters:
        - embedding: np.ndarray
            The feature matrix for UMAP embedding.
        - anndata: anndata.AnnData
            The AnnData object containing metadata, including clone_id labels.
        - num_top_clone_ids: int, optional (default=5)
            Number of most frequent clone IDs to highlight.
        - custom_palette: list, optional
            List of color hex codes for the highlighted clone IDs.
        """
        self.embedding = embedding
        self.anndata = anndata
        self.num_top_clone_ids = num_top_clone_ids
        self.custom_palette = custom_palette or ['#D1392C', '#4A7CB3', '#67AD57', '#8E529F', '#EE8632']
        self.labels = anndata.obs["clone_id"].to_numpy()
        self.clone_id_counts = anndata.obs['clone_id'].value_counts()
        self.top_clone_ids = self.clone_id_counts.index[:num_top_clone_ids]
        self.knn_train_model = None
        self.knn_test_model = None
        
    def plot_umap_top_lin(self, title="UMAP Plot"):
        """
        Plots a UMAP embedding with the top clone IDs highlighted.
        
        Parameters:
        - title: str, optional (default="UMAP Plot")
            Title for the plot.
        """
        # Initialize and fit UMAP
        reducer = umap.UMAP()
        embedding = reducer.fit_transform(self.embedding)
        
        # Create color map
        color_map = {clone_id: 'lightgray' for clone_id in self.labels}
        for i, clone_id in enumerate(self.top_clone_ids):
            if i < len(self.custom_palette):
                color_map[clone_id] = self.custom_palette[i]
        
        # Separate indices for gray and colored cells
        gray_indices = [i for i, clone_id in enumerate(self.labels) if color_map[clone_id] == 'lightgray']
        color_indices = [i for i, clone_id in enumerate(self.labels) if color_map[clone_id] != 'lightgray']
        
        # Plot
        plt.figure(figsize=(10, 5.85), dpi=300)
        plt.grid(True, alpha=0.9, linestyle=':')
        plt.gca().tick_params(axis='both', which='both', length=0)
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
        
        # Plot gray cells first
        plt.scatter(embedding[gray_indices, 0], embedding[gray_indices, 1], c='gray', s=1, label='_nolegend_', alpha=0.7)
        
        # Overlay colored cells
        for clone_id in self.top_clone_ids:
            color_mask = np.array([clone_id == cid for cid in self.labels])
            plt.scatter(embedding[color_mask, 0], embedding[color_mask, 1], 
                        c=[color_map[clone_id]], s=40, label=clone_id)
        
        # Create custom legend
        legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[clone_id], 
                                    markersize=10, label=f'{clone_id}') for clone_id in self.top_clone_ids]
        plt.legend(handles=legend_handles, title="Top Clone IDs")
        
        # Set title
        plt.title(title)
        
        plt.show()
    
    def calculate_calinski_harabasz_score(self):
        """
        Calculates the Calinski-Harabasz score for the given embedding and labels.
        
        Returns:
        - score: float
            The Calinski-Harabasz score.
        """
        score = calinski_harabasz_score(self.embedding, self.labels)
        print("Calinski-Harabasz Score:", score)
        return score
    
    def KNN_train(self, n_neighbors=10, test_size=0.2, random_state=42):
        """
        Train a KNN classifier on a subset of the training embeddings.
        """
        X_train, X_test, y_train, y_test = train_test_split(self.embedding, self.labels, test_size=test_size, random_state=random_state)
        self.knn_train_model = KNeighborsClassifier(n_neighbors=n_neighbors)
        self.knn_train_model.fit(X_train, y_train)
        
        y_pred = self.knn_train_model.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        print(f"KNN classifier training accuracy: {accuracy * 100:.2f}%")
        return accuracy
    
    def KNN_test(self, test_embedding, anndata_test, n_neighbors=10):
        """
        Train a KNN classifier on the full training set and evaluate on the test set.
        """
        test_labels = anndata_test.obs["clone_id"].to_numpy()
        self.knn_test_model = KNeighborsClassifier(n_neighbors=n_neighbors)
        self.knn_test_model.fit(self.embedding, self.labels)
        
        y_pred = self.knn_test_model.predict(test_embedding)
        accuracy = accuracy_score(test_labels, y_pred)
        print(f"KNN classifier testing accuracy: {accuracy * 100:.2f}%")
        return accuracy
