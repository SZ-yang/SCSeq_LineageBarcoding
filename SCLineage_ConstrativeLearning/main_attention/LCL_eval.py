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
        Initialize Eval with embedding and anndata.
        """
        self.embedding = embedding
        self.anndata = anndata
        self.num_top_clone_ids = num_top_clone_ids
        self.custom_palette = custom_palette or ['#D1392C', '#4A7CB3', '#67AD57', '#8E529F', '#EE8632']
        self.labels = anndata.obs["clone_id"].to_numpy()
        self.clone_id_counts = anndata.obs['clone_id'].value_counts()
        self.top_clone_ids = self.clone_id_counts.index[:num_top_clone_ids]
        # placeholders for KNN scores
        self.train_score = None
        self.test_score = None
        self.knn_train_model = None
        self.knn_test_model = None
        
    def plot_umap_top_lin(self, title="UMAP Plot"):
        # ... your plotting code unchanged ...
        reducer = umap.UMAP()
        embedding = reducer.fit_transform(self.embedding)
        color_map = {cid: 'lightgray' for cid in self.labels}
        for i, cid in enumerate(self.top_clone_ids):
            if i < len(self.custom_palette):
                color_map[cid] = self.custom_palette[i]
        gray_idx = [i for i, cid in enumerate(self.labels) if color_map[cid]=='lightgray']
        plt.figure(figsize=(10,5.85), dpi=300)
        plt.grid(True, alpha=0.9, linestyle=':')
        plt.gca().tick_params(axis='both', which='both', length=0)
        plt.gca().set_xticklabels([]); plt.gca().set_yticklabels([])
        plt.scatter(embedding[gray_idx,0], embedding[gray_idx,1], c='gray', s=1, alpha=0.7)
        for cid in self.top_clone_ids:
            mask = self.labels == cid
            plt.scatter(embedding[mask,0], embedding[mask,1],
                        c=[color_map[cid]], s=40, label=cid)
        handles = [plt.Line2D([0],[0], marker='o', color='w',
                              markerfacecolor=color_map[cid],
                              markersize=10, label=str(cid))
                   for cid in self.top_clone_ids]
        plt.legend(handles=handles, title="Top Clone IDs")
        plt.title(title)
        plt.show()
    
    def calculate_calinski_harabasz_score(self):
        score = calinski_harabasz_score(self.embedding, self.labels)
        print("Calinski-Harabasz Score:", score)
        return score
    
    def KNN_train(self, n_neighbors=10, test_size=0.2, random_state=42):
        """
        Train a KNN on a train/test split of self.embedding/labels,
        stores train_score and returns it.
        """
        X_tr, X_va, y_tr, y_va = train_test_split(
            self.embedding, self.labels,
            test_size=test_size, random_state=random_state
        )
        self.knn_train_model = KNeighborsClassifier(n_neighbors=n_neighbors)
        self.knn_train_model.fit(X_tr, y_tr)
        y_pred = self.knn_train_model.predict(X_va)
        acc = accuracy_score(y_va, y_pred)
        print(f"KNN training accuracy: {acc*100:.2f}%")
        self.train_score = acc
        return acc
    
    def KNN_test(self, test_embedding, anndata_test, n_neighbors=10):
        """
        Fit KNN on the entire training set, evaluate on test_embedding/test_labels,
        stores test_score and returns it.
        """
        test_labels = anndata_test.obs["clone_id"].to_numpy()
        self.knn_test_model = KNeighborsClassifier(n_neighbors=n_neighbors)
        self.knn_test_model.fit(self.embedding, self.labels)
        y_pred = self.knn_test_model.predict(test_embedding)
        acc = accuracy_score(test_labels, y_pred)
        print(f"KNN testing accuracy: {acc*100:.2f}%")
        self.test_score = acc
        return acc
    
    def scores(self):
        """
        Return (train_score, test_score). Either may be None if not yet computed.
        """
        return (self.train_score, self.test_score)