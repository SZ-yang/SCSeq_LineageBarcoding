import numpy as np
import umap
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from typing import Tuple, List, Dict

class LCL_Eval:
    """
    Given a concatenated AnnData (train+test with embedding in .obsm["LCL_embedding"]), this class provides:
      1) compute_adjusted_knn_stats(...) → returns train/test KNN accuracy + rank‐based stats
      2) plot_top_clones_umap(...)     → UMAP highlighting top‐N clones (train vs test)
      3) plot_test_accuracy_umap(...)  → UMAP where all train are gray, and test cells = green/red by correctness
      4) run_all(...)                  → convenience wrapper to run all three at once
    """

    def __init__(
        self,
        adata,
        clone_key: str = "clone_id",
        dataset_key: str = "dataset",
        num_top: int = 5,
        palette: list = None,
        umap_kwargs: dict = None
    ):
        """
        Parameters:
        -----------
        adata : anndata.AnnData
            Must contain:
              - adata.obs[clone_key]      (per‐cell lineage/clone ID)
              - adata.obs[dataset_key]    (either "train" or "test")
              - adata.obsm["LCL_embedding"] (high‐D LCL embedding: shape [n_cells, D], D≥3)
        clone_key : str
            Column in adata.obs that holds each cell’s lineage ID.
        dataset_key : str
            Column in adata.obs that indicates "train" vs. "test".
        num_top : int
            Number of top lineages (by count) to highlight in `plot_top_clones_umap`.
        palette : list[str], optional
            List of ≥ num_top color names/hex codes; if None, uses seaborn “tab10.”
        umap_kwargs : dict, optional
            Keyword‐args for `umap.UMAP(**umap_kwargs)` when computing 2D.
        """
        # Copy AnnData so as not to modify the original in‐place:
        self.adata = adata.copy()
        self.embedding_key = "LCL_embedding"
        self.clone_key = clone_key
        self.dataset_key = dataset_key
        self.num_top = num_top
        self.umap_kwargs = umap_kwargs or {"random_state": 42}

        # Enforce seaborn style
        sns.set_style("whitegrid")

        # 1) Verify high‐D embedding is present:
        if self.embedding_key not in self.adata.obsm:
            raise KeyError(f"Could not find embedding under adata.obsm['{self.embedding_key}'].")

        hd = self.adata.obsm[self.embedding_key]
        if hd.ndim != 2 or hd.shape[1] < 3:
            raise ValueError(f"adata.obsm['{self.embedding_key}'] must be shape [n_cells, D] with D≥3.")

        # 2) Compute a 2D UMAP projection from that high‐D, store into obsm["UMAP_embedding"]:
        reducer = umap.UMAP(**self.umap_kwargs)
        umap2d = reducer.fit_transform(hd)  # (n_cells, 2)
        self.adata.obsm["UMAP_embedding"] = umap2d

        # 3) Identify the top‐N clones by frequency across the entire AnnData
        counts = self.adata.obs[self.clone_key].value_counts()
        self.top_clones = counts.index[: self.num_top].tolist()

        # 4) Build a new column "clone_group": top‐N clones remain as is; all others → "Other"
        def _group_fn(x):
            return x if x in self.top_clones else "Other"

        self.adata.obs["clone_group"] = (
            self.adata.obs[self.clone_key]
            .apply(_group_fn)
            .astype("category")
        )

        # 5) Reorder categories so that the top‐N clones appear first, then “Other”
        cats = self.top_clones + ["Other"]
        self.adata.obs["clone_group"] = (
            self.adata.obs["clone_group"]
                .cat
                .reorder_categories(cats, ordered=True)
        )

        # 6) Build a color map for those top‐N clones. If no palette given, use seaborn "tab10"
        if palette is None:
            base_pal = sns.color_palette("tab10", n_colors=self.num_top)
            self.color_map = {self.top_clones[i]: base_pal[i] for i in range(self.num_top)}
        else:
            if len(palette) < self.num_top:
                raise ValueError(f"Provided palette has length {len(palette)} but num_top={self.num_top}.")
            self.color_map = {self.top_clones[i]: palette[i] for i in range(self.num_top)}

        # Always draw "Other" as light gray:
        self.color_map["Other"] = "lightgray"


    # ──────────────────── KNN & rank‐based evaluation methods ──────────────────────────

    @staticmethod
    def compute_global_freq(labels: np.ndarray) -> Dict[int, float]:
        """
        Utility: compute global frequency of each label in the array `labels`.
        """
        unique, counts = np.unique(labels, return_counts=True)
        total = len(labels)
        return {lab: cnt / total for lab, cnt in zip(unique, counts)}


    @staticmethod
    def adjusted_knn_predict_with_rank(
        knn_model: KNeighborsClassifier,
        train_labels: np.ndarray,
        embeddings: np.ndarray,
        labels_true: np.ndarray,
        global_freq: Dict[int, float],
        k: int
    ) -> Tuple[np.ndarray, List[int], float, float, List[int]]:
        """
        For each point in `embeddings` (which may be training‐set or test‐set embeddings):
          1) find its k nearest TRAINING neighbors (via knn_model.kneighbors)
          2) compute local_freq[L] = (# neighbors with label L) / k
          3) adjusted_score[L] = local_freq[L] - global_freq[L]
          4) predict = argmax_L (adjusted_score[L])
          5) compute rank_of_true_label = position of true_label in descending local_freq order
             (if true_label not among neighbors, rank = num_unique_labels + 1)

        Returns:
          preds          : np.ndarray, length = n_points
          ranks          : List[int], each the rank of the true_label among that point’s neighbors
          avg_rank       : float = mean(ranks)
          avg_unique_lbl : float = average # distinct labels among each point’s k neighbors
          correct_flags  : List[int], 1 if predicted == true_label else 0 each point
        """
        neigh_indices = knn_model.kneighbors(embeddings, return_distance=False)  # shape (n_points, k)
        n_points = embeddings.shape[0]

        preds = np.zeros(n_points, dtype=train_labels.dtype)
        ranks = []
        unique_label_counts = []
        correct_flags = []

        for i, nbrs in enumerate(neigh_indices):
            nbr_labels = train_labels[nbrs]
            unique, counts = np.unique(nbr_labels, return_counts=True)
            num_unique = len(unique)
            unique_label_counts.append(num_unique)

            # local frequency among neighbors
            local_freq = {lab: cnt / k for lab, cnt in zip(unique, counts)}

            # build adjusted scores (for every label that appears globally—zero if not in local_freq)
            scores = {
                lab: local_freq.get(lab, 0.0) - global_freq.get(lab, 0.0)
                for lab in global_freq.keys()
            }

            # predicted label = argmax adjusted score
            best_label = max(scores.items(), key=lambda x: x[1])[0]
            preds[i] = best_label

            # find rank of the true_label among local_freq
            sorted_labels = sorted(local_freq.items(), key=lambda x: x[1], reverse=True)
            ranked_labs = [lab for lab, _ in sorted_labels]
            true_lab = labels_true[i]
            if true_lab in ranked_labs:
                rank_pos = ranked_labs.index(true_lab) + 1  # 1-based
            else:
                rank_pos = num_unique + 1
            ranks.append(rank_pos)

            # correct? 1 if pred == true
            correct_flags.append(1 if best_label == true_lab else 0)

        avg_rank = float(np.mean(ranks))
        avg_unique_lbl = float(np.mean(unique_label_counts))

        return preds, ranks, avg_rank, avg_unique_lbl, correct_flags


    def evaluate_adjusted_knn(
        self,
        train_embeddings: np.ndarray,
        train_labels: np.ndarray,
        test_embeddings: np.ndarray,
        test_labels: np.ndarray,
        k: int = 30
    ) -> Dict[str, Dict]:
        """
        Fit an adjusted‐KNN on (train_embeddings, train_labels), then evaluate on both train & test.

        Returns a dict of the form:
        {
          "train": {
              "overall_accuracy"     : float,  # Layer‐2 average of per‐lineage accuracies
              "overall_avg_rank"     : float,  # Layer‐2 average of per‐lineage average ranks
              "overall_avg_unique"   : float,  # Layer‐2 average of per‐lineage avg unique labels
              "per_lineage_accuracy" : Dict[label, float],
              "per_lineage_avg_rank" : Dict[label, float],
              "rank_quantiles"       : Dict[qname, float]  # {"q25":…, "q50":…, "q75":…}
          },
          "test": { … same keys … }
        }
        """
        # 1) Global frequencies from training labels only
        global_freq = self.compute_global_freq(train_labels)

        # 2) Fit a standard KNN on the training data
        knn = KNeighborsClassifier(n_neighbors=k)
        knn.fit(train_embeddings, train_labels)

        # 3a) Run adjusted KNN prediction on the TRAIN set itself
        (
          train_preds,
          train_ranks,
          train_avg_rank,
          train_avg_unique,
          train_correct_flags
        ) = self.adjusted_knn_predict_with_rank(
                knn_model=knn,
                train_labels=train_labels,
                embeddings=train_embeddings,
                labels_true=train_labels,
                global_freq=global_freq,
                k=k
            )

        # 3b) Run adjusted KNN prediction on the TEST set
        (
          test_preds,
          test_ranks,
          test_avg_rank,
          test_avg_unique,
          test_correct_flags
        ) = self.adjusted_knn_predict_with_rank(
                knn_model=knn,
                train_labels=train_labels,
                embeddings=test_embeddings,
                labels_true=test_labels,
                global_freq=global_freq,
                k=k
            )

        # 4) Helper to perform two‐layer aggregation (cells → lineages → overall)
        def layer2_aggregate(
            cell_labels: np.ndarray,
            correct_flags: List[int],
            ranks: List[int]
        ) -> Tuple[
               Dict[int, float],    # per‐lineage accuracy
               Dict[int, float],    # per‐lineage avg_rank
               float,               # overall accuracy (layer‐2)
               float,               # overall avg_rank (layer‐2)
               List[float]          # vector of per‐lineage avg_rank, sorted by lineage label
           ]:
            per_lin_acc = {}
            per_lin_rank = {}

            unique_lineages = np.unique(cell_labels)
            for lab in unique_lineages:
                idxs = np.where(cell_labels == lab)[0]
                lab_acc = np.mean([correct_flags[i] for i in idxs])
                lab_avg_rank = np.mean([ranks[i] for i in idxs])
                per_lin_acc[lab] = float(lab_acc)
                per_lin_rank[lab] = float(lab_avg_rank)

            # overall (layer‐2) averages (each lineage counts equally)
            overall_acc = float(np.mean(list(per_lin_acc.values())))
            overall_rank = float(np.mean(list(per_lin_rank.values())))

            # produce a list of per‐lineage average ranks in sorted‐label order
            lineage_rank_list = [per_lin_rank[lab] for lab in sorted(unique_lineages)]

            return per_lin_acc, per_lin_rank, overall_acc, overall_rank, lineage_rank_list

        # 5) Aggregate train‐side results
        (
          train_per_lin_acc,
          train_per_lin_rank,
          train_overall_acc,
          train_overall_rank,
          train_lineage_rank_list
        ) = layer2_aggregate(
                cell_labels=train_labels,
                correct_flags=train_correct_flags,
                ranks=train_ranks
            )

        # 6) Aggregate test‐side results
        (
          test_per_lin_acc,
          test_per_lin_rank,
          test_overall_acc,
          test_overall_rank,
          test_lineage_rank_list
        ) = layer2_aggregate(
                cell_labels=test_labels,
                correct_flags=test_correct_flags,
                ranks=test_ranks
            )

        # 7) Compute 25%, 50%, 75% quantiles over per‐lineage average ranks (rounded to 3 decimals)
        def compute_quantiles(rank_list: List[float]) -> Dict[str, float]:
            arr = np.array(rank_list)
            q25 = round(float(np.quantile(arr, 0.25)), 3)
            q50 = round(float(np.quantile(arr, 0.50)), 3)
            q75 = round(float(np.quantile(arr, 0.75)), 3)
            return {"q25": q25, "q50": q50, "q75": q75}

        train_quantiles = compute_quantiles(train_lineage_rank_list)
        test_quantiles  = compute_quantiles(test_lineage_rank_list)

        # 8) Package into one dict
        results = {
          "train": {
            "overall_accuracy":        round(train_overall_acc, 3),
            "overall_avg_rank":        round(train_overall_rank, 3),
            "overall_avg_unique":      round(train_avg_unique, 3),
            "rank_quantiles":          train_quantiles
          },
          "test": {
            "overall_accuracy":        round(test_overall_acc, 3),
            "overall_avg_rank":        round(test_overall_rank, 3),
            "overall_avg_unique":      round(test_avg_unique, 3),
            "rank_quantiles":          test_quantiles
          }
        }

        return results


    # ──────────────────── UMAP PLOTTING METHODS ──────────────────────────

    def plot_top_clones_umap(self, figsize=(8, 6), title=None, savepath: str = None):
        """
        Draw a UMAP scatterplot (train+test) highlighting top‐N lineages:
          - "Other" cells in light gray (low alpha)
          - For each top lineage:
                • train‐cells: marker=".", size=30, alpha=0.8
                • test‐cells:  marker="x", size=40, alpha=0.9
        """
        df = self.adata.obs
        coords = self.adata.obsm["UMAP_embedding"]

        is_train = (df[self.dataset_key] == "train")
        is_test  = (df[self.dataset_key] == "test")

        fig, ax = plt.subplots(figsize=figsize, dpi=200)

        # 1) Draw “Other” second‐tier behind everything (train + test)
        mask_train_other = (df["clone_group"] == "Other") & is_train
        mask_test_other  = (df["clone_group"] == "Other") & is_test

        ax.scatter(
            coords[mask_train_other, 0],
            coords[mask_train_other, 1],
            color=self.color_map["Other"],
            s=8,
            marker=".",
            alpha=0.2,
            label="Train Other"
        )
        ax.scatter(
            coords[mask_test_other, 0],
            coords[mask_test_other, 1],
            color=self.color_map["Other"],
            s=12,
            marker="x",
            alpha=0.2,
            label="Test Other"
        )

        # 2) Overlay each of the top‐N clones
        for clone in self.top_clones:
            mask_train_clone = (df["clone_group"] == clone) & is_train
            mask_test_clone  = (df["clone_group"] == clone) & is_test

            # train = dot
            ax.scatter(
                coords[mask_train_clone, 0],
                coords[mask_train_clone, 1],
                color=self.color_map[clone],
                s=30,
                marker=".",
                alpha=0.8,
                label=f"Train {clone}"
            )
            # test = cross
            ax.scatter(
                coords[mask_test_clone, 0],
                coords[mask_test_clone, 1],
                color=self.color_map[clone],
                s=40,
                marker="x",
                alpha=0.9,
                label=f"Test {clone}"
            )

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title or f"UMAP – Top {self.num_top} Clones")

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(
            by_label.values(),
            by_label.keys(),
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            fontsize="small"
        )

        plt.tight_layout()
        if savepath:
            fig.savefig(savepath, bbox_inches="tight")
        return fig, ax


    def plot_test_accuracy_umap(
        self,
        train_embeddings: np.ndarray,
        train_labels: np.ndarray,
        test_embeddings: np.ndarray,
        test_labels: np.ndarray,
        k: int = 30,
        figsize=(8, 6),
        title=None,
        savepath: str = None
    ):
        """
        Draw a UMAP scatterplot (train+test), coloring:
          - All TRAIN cells in light gray
          - TEST cells in green if adjusted‐KNN predicted correctly, or red if incorrectly.

        Requires calling `evaluate_adjusted_knn(...)` internally to get per‐cell correctness flags.

        Returns fig, ax.
        """
        # 1) Run adjusted KNN to get correctness flags (only need per‐cell correct_flags for test set)
        stats = self.evaluate_adjusted_knn(
            train_embeddings=train_embeddings,
            train_labels=train_labels,
            test_embeddings=test_embeddings,
            test_labels=test_labels,
            k=k
        )

        # We need `test_correct_flags` per‐cell. Unfortunately, evaluate_adjusted_knn
        # returns only aggregated lineage‐level data.
        # So we’ll re‐run adjusted_knn_predict_with_rank(...) directly to get per‐cell flags:
        global_freq = self.compute_global_freq(train_labels)
        knn_model = KNeighborsClassifier(n_neighbors=k)
        knn_model.fit(train_embeddings, train_labels)

        _, _, _, _, test_correct_flags = self.adjusted_knn_predict_with_rank(
            knn_model=knn_model,
            train_labels=train_labels,
            embeddings=test_embeddings,
            labels_true=test_labels,
            global_freq=global_freq,
            k=k
        )
        test_correct_flags = np.array(test_correct_flags)  # length = n_test

        # 2) Mark each cell in the concatenated adata whether it’s train/test and if test‐cell correct
        df = self.adata.obs.copy()
        coords = self.adata.obsm["UMAP_embedding"]

        # create a boolean array of length = total cells in adata:
        is_train = (df[self.dataset_key] == "train")
        is_test  = (df[self.dataset_key] == "test")

        # Build an array of the same length: 0 = train, 1 = test‐correct, 2 = test‐incorrect
        status = np.zeros(len(df), dtype=int)
        status[is_train] = 0

        # Now fill in test positions in the order they appear:
        test_indices = np.where(is_test)[0]
        for idx_in_ad, correct_flag in zip(test_indices, test_correct_flags):
            status[idx_in_ad] = 1 if correct_flag == 1 else 2

        # 3) Prepare colors & plotting
        #   0  → train = light gray
        #   1  → test‐correct = green
        #   2  → test‐incorrect = red
        color_map = {0: "lightgray", 1: "green", 2: "red"}
        size_map  = {0: 12,         1: 30,      2: 30}
        marker_map= {0: ".",       1: "x",     2: "x"}
        alpha_map = {0: 0.2,       1: 0.9,     2: 0.9}

        fig, ax = plt.subplots(figsize=figsize, dpi=200)
        for status_val in [0, 1, 2]:
            mask = (status == status_val)
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=color_map[status_val],
                s=size_map[status_val],
                marker=marker_map[status_val],
                alpha=alpha_map[status_val],
                label=(
                    "Train" if status_val == 0 else
                    "Test Correct" if status_val == 1 else
                    "Test Incorrect"
                )
            )

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title or "UMAP – Test‐cell Accuracy (green=correct, red=incorrect)")

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(
            by_label.values(),
            by_label.keys(),
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            fontsize="small"
        )

        plt.tight_layout()
        if savepath:
            fig.savefig(savepath, bbox_inches="tight")
        return fig, ax


    # ──────────────────── Convenience wrapper ──────────────────────────

    def run_all(
        self,
        train_embeddings: np.ndarray,
        train_labels: np.ndarray,
        test_embeddings: np.ndarray,
        test_labels: np.ndarray,
        k: int = 30,
        umap1_kwargs: dict = None,
        umap2_kwargs: dict = None
    ) -> Tuple[Dict[str, Dict], Tuple[plt.Figure, plt.Axes], Tuple[plt.Figure, plt.Axes]]:
        """
        Convenience wrapper that runs, in order:
          1) KNN‐based evaluation → returns stats dict
          2) top‐clones UMAP → returns (fig1, ax1)
          3) test‐accuracy UMAP → returns (fig2, ax2)

        Parameters:
        -----------
        train_embeddings, train_labels, test_embeddings, test_labels : np.ndarray, np.ndarray
            The embeddings + labels needed by evaluate_adjusted_knn and plot_test_accuracy_umap.
        k : int
            # neighbors for the adjusted KNN.
        umap1_kwargs, umap2_kwargs: dict or None
            (Unused here—both plotting methods share the same UMAP that was pre‐computed
             in __init__, so you can ignore them. They exist in case you want to pass further
             styling instructions in the future.)

        Returns:
        --------
        (knn_stats, (fig1, ax1), (fig2, ax2))
          - knn_stats:       the dict returned by evaluate_adjusted_knn(...)
          - (fig1, ax1):     figure/axis from plot_top_clones_umap(...)
          - (fig2, ax2):     figure/axis from plot_test_accuracy_umap(...)
        """
        knn_stats = self.evaluate_adjusted_knn(
            train_embeddings=train_embeddings,
            train_labels=train_labels,
            test_embeddings=test_embeddings,
            test_labels=test_labels,
            k=k
        )

        fig1, ax1 = self.plot_top_clones_umap()
        fig2, ax2 = self.plot_test_accuracy_umap(
            train_embeddings=train_embeddings,
            train_labels=train_labels,
            test_embeddings=test_embeddings,
            test_labels=test_labels,
            k=k
        )

        return knn_stats, (fig1, ax1), (fig2, ax2)


# EOF of lcl_full_eval.py