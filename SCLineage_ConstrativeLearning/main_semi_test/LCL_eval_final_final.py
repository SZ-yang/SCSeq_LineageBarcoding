# lcl_full_eval.py

import os
import numpy as np
import umap
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.neighbors import KNeighborsClassifier
from typing import Tuple, List, Dict

class LCL_Eval:
    """
    Given a concatenated AnnData (train+test with embedding in .obsm["LCL_embedding"]), this class provides:
      1) evaluate_adjusted_knn(...) → returns train/test KNN accuracy, containment rate,
                                        conditional‐rank stats + unique‐label stats + accuracy‐quantiles
      2) plot_top_clones_umap(...)   → UMAP highlighting top‐N clones (train vs test)
      3) plot_test_accuracy_umap(...)→ UMAP: train=gray, test correct=green (on top), test incorrect=red
      4) run_all(...)                → convenience wrapper to run all three at once
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
        # copy so we don’t modify original
        self.adata = adata.copy()
        self.embedding_key = "LCL_embedding"
        self.clone_key    = clone_key
        self.dataset_key  = dataset_key
        self.num_top      = num_top
        self.umap_kwargs  = umap_kwargs or {"random_state": 42}

        sns.set_style("whitegrid")

        # 1) verify high‐D embedding
        if self.embedding_key not in self.adata.obsm:
            raise KeyError(f"Could not find embedding under adata.obsm['{self.embedding_key}'].")

        hd = self.adata.obsm[self.embedding_key]
        if hd.ndim != 2 or hd.shape[1] < 3:
            raise ValueError(f"adata.obsm['{self.embedding_key}'] must be shape [n_cells, D] with D≥3.")

        # 2) compute 2D UMAP
        reducer = umap.UMAP(**self.umap_kwargs)
        self.adata.obsm["UMAP_embedding"] = reducer.fit_transform(hd)

        # 3) find top‐N clones
        counts = self.adata.obs[self.clone_key].value_counts()
        self.top_clones = counts.index[: self.num_top].tolist()

        # 4) build clone_group column
        def _group_fn(x):
            return x if x in self.top_clones else "Other"

        self.adata.obs["clone_group"] = (
            self.adata.obs[self.clone_key]
            .apply(_group_fn)
            .astype("category")
        )
        # reorder so top clones first
        cats = self.top_clones + ["Other"]
        self.adata.obs["clone_group"] = (
            self.adata.obs["clone_group"]
                .cat
                .reorder_categories(cats, ordered=True)
        )

        # 5) build color map
        if palette is None:
            base_pal = sns.color_palette("tab10", n_colors=self.num_top)
            self.color_map = {self.top_clones[i]: base_pal[i] for i in range(self.num_top)}
        else:
            if len(palette) < self.num_top:
                raise ValueError(f"Palette length {len(palette)} < num_top={self.num_top}")
            self.color_map = {self.top_clones[i]: palette[i] for i in range(self.num_top)}
        self.color_map["Other"] = "lightgray"


    # ─── KNN & rank/containment evaluation ───────────────────────────

    @staticmethod
    def compute_global_freq(labels: np.ndarray) -> Dict[int, float]:
        unique, counts = np.unique(labels, return_counts=True)
        total = len(labels)
        return {lab: cnt/total for lab, cnt in zip(unique, counts)}


    @staticmethod
    def adjusted_knn_predict_with_rank(
        knn_model: KNeighborsClassifier,
        train_labels: np.ndarray,
        embeddings: np.ndarray,
        labels_true: np.ndarray,
        global_freq: Dict[int, float],
        k: int
    ) -> Tuple[
         np.ndarray,        # preds
         List[float],       # ranks (np.nan if not contained)
         List[int],         # containment flags (1 if true_label in k-NN)
         List[int],         # correct flags
         List[int]          # unique_counts per-cell
    ]:
        """
        Returns:
          preds              : np.ndarray, length = n_points
          ranks              : List[float], rank of true label if contained, else np.nan
          containment_flags  : List[int], 1 if true_label in k‐NN else 0
          correct_flags      : List[int], 1 if predicted==true_label else 0
          unique_counts      : List[int], # distinct train‐labels in each k‐NN
        """
        neigh = knn_model.kneighbors(embeddings, return_distance=False)
        n = embeddings.shape[0]

        preds            = np.zeros(n, dtype=train_labels.dtype)
        ranks            = []
        containment      = []
        correct          = []
        unique_counts    = []

        for i, nbrs in enumerate(neigh):
            nbr_lab = train_labels[nbrs]
            uq, ct  = np.unique(nbr_lab, return_counts=True)
            unique_counts.append(len(uq))

            # local freq
            local_freq = {lab: c/k for lab, c in zip(uq, ct)}

            # adjusted scores
            scores = {
                lab: local_freq.get(lab,0.0) - global_freq.get(lab,0.0)
                for lab in global_freq.keys()
            }
            best = max(scores.items(), key=lambda x: x[1])[0]
            preds[i] = best

            # containment & rank
            sorted_lbl = sorted(local_freq.items(), key=lambda x: x[1], reverse=True)
            ranked     = [lab for lab,_ in sorted_lbl]
            t          = labels_true[i]
            if t in ranked:
                containment.append(1)
                ranks.append(ranked.index(t) + 1)
            else:
                containment.append(0)
                ranks.append(np.nan)

            # correctness
            correct.append(1 if best == t else 0)

        return preds, ranks, containment, correct, unique_counts


    def evaluate_adjusted_knn(
        self,
        train_embeddings: np.ndarray,
        train_labels: np.ndarray,
        test_embeddings: np.ndarray,
        test_labels: np.ndarray,
        k: int = 30
    ) -> Dict[str, Dict]:
        """
        Fit adjusted‐KNN on train, then evaluate:

        Returns dict:
        {
          "train": {
            "accuracy",    # standard adjusted‐knn accuracy
            "avg_unique"   # avg #unique labels in k‐NN (layer‐1)
          },
          "test": {
            "accuracy",          # adjusted‐knn accuracy
            "containment_rate",  # P(true lineage in k‐NN)
            "overall_avg_rank",  # avg rank *conditional on containment* (layer‐2 over lineages)
            "avg_unique",        # avg #unique labels in k‐NN (layer‐1)
            "unique_quantiles",  # q0,q25,q50,q75,q100 over per‐lineage avg‐unique
            "accuracy_quantiles" # q0,q25,q50,q75,q100 over per‐lineage accuracy
          }
        }
        """
        gf  = self.compute_global_freq(train_labels)
        knn = KNeighborsClassifier(n_neighbors=k).fit(train_embeddings, train_labels)

        # ── train side ───────────────────────────────────────────────────────────
        _, train_ranks, train_cont, train_corr, train_uniq = (
            self.adjusted_knn_predict_with_rank(
                knn, train_labels, train_embeddings, train_labels, gf, k
            )
        )
        train_accuracy   = float(np.mean(train_corr))
        train_avg_unique = float(np.mean(train_uniq))

        # ── test side ────────────────────────────────────────────────────────────
        _, test_ranks, test_cont, test_corr, test_uniq = (
            self.adjusted_knn_predict_with_rank(
                knn, train_labels, test_embeddings, test_labels, gf, k
            )
        )
        test_accuracy    = float(np.mean(test_corr))
        containment_rate = float(np.mean(test_cont))
        avg_unique_test  = float(np.mean(test_uniq))

        # ── conditional‐rank (layer‐2 over lineages) ───────────────────────────────
        mask = np.array(test_cont, dtype=bool)
        # build per‐lineage avg‐rank conditional on containment
        per_rank = {}
        for lab in np.unique(test_labels[mask]):
            idxs = np.where((test_labels == lab) & mask)[0]
            per_rank[lab] = float(np.nanmean([test_ranks[i] for i in idxs]))
        overall_avg_rank = float(np.mean(list(per_rank.values())))

        # quantiles helper
        def qdict(arr: np.ndarray) -> Dict[str, float]:
            return {
              "q0":   round(float(np.quantile(arr, 0.0)), 3),
              "q25":  round(float(np.quantile(arr, 0.25)),3),
              "q50":  round(float(np.quantile(arr, 0.50)),3),
              "q75":  round(float(np.quantile(arr, 0.75)),3),
              "q100": round(float(np.quantile(arr, 1.0)),3)
            }

        # unique‐quantiles over per‐lineage avg‐unique (layer‐1 → layer‐2)
        # (we already have test_uniq per‐cell, so first aggregate per‐lineage)
        per_uniq = {}
        for lab in np.unique(test_labels):
            idxs = np.where(test_labels == lab)[0]
            per_uniq[lab] = float(np.mean([test_uniq[i] for i in idxs]))
        unique_quantiles = qdict(np.array(list(per_uniq.values())))

        # accuracy‐quantiles over per‐lineage accuracy
        per_acc = {}
        for lab in np.unique(test_labels):
            idxs = np.where(test_labels == lab)[0]
            per_acc[lab] = float(np.mean([test_corr[i] for i in idxs]))
        accuracy_quantiles = qdict(np.array(list(per_acc.values())))

        return {
          "train": {
            "accuracy":   round(train_accuracy,3),
            "avg_unique": round(train_avg_unique,3)
          },
          "test": {
            "accuracy":           round(test_accuracy,3),
            "containment_rate":   round(containment_rate,3),
            "overall_avg_rank":   round(overall_avg_rank,3),
            "avg_unique":         round(avg_unique_test,3),
            "unique_quantiles":   unique_quantiles,
            "accuracy_quantiles": accuracy_quantiles
          }
        }


    # ─── UMAP PLOTS ────────────────────────────────────────────────────

    def plot_top_clones_umap(
        self,
        figsize=(8,6),
        title=None,
        savepath: str = None
    ):
        df     = self.adata.obs
        coords = self.adata.obsm["UMAP_embedding"]
        is_tr  = df[self.dataset_key] == "train"
        is_te  = df[self.dataset_key] == "test"

        fig, ax = plt.subplots(figsize=figsize, dpi=200)

        # Other first
        m1 = (df["clone_group"]=="Other") & is_tr
        m2 = (df["clone_group"]=="Other") & is_te
        ax.scatter(coords[m1,0], coords[m1,1],
                   c=self.color_map["Other"], s=8,  marker=".", alpha=0.2, label="Train Other")
        ax.scatter(coords[m2,0], coords[m2,1],
                   c=self.color_map["Other"], s=12, marker="x", alpha=0.2, label="Test Other")

        # top clones
        for c in self.top_clones:
            mt = (df["clone_group"]==c) & is_tr
            mx = (df["clone_group"]==c) & is_te
            ax.scatter(coords[mt,0], coords[mt,1],
                       c=self.color_map[c], s=30, marker=".", alpha=0.8, label=f"Train {c}")
            ax.scatter(coords[mx,0], coords[mx,1],
                       c=self.color_map[c], s=40, marker="x", alpha=0.9, label=f"Test {c}")

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title or f"UMAP – Top {self.num_top} Clones")

        handles, labels = ax.get_legend_handles_labels()
        by = dict(zip(labels, handles))
        ax.legend(by.values(), by.keys(),
                  bbox_to_anchor=(1.02, 1), loc="upper left",
                  frameon=False, fontsize="small")

        plt.tight_layout()
        if savepath:
            savepath = os.path.splitext(savepath)[0] + ".png"
            fig.savefig(savepath, format="png", dpi=300, bbox_inches="tight")
        return fig, ax


    def plot_test_accuracy_umap(
        self,
        train_embeddings: np.ndarray,
        train_labels: np.ndarray,
        test_embeddings: np.ndarray,
        test_labels: np.ndarray,
        k: int = 30,
        figsize=(8,6),
        title=None,
        savepath: str = None
    ):
        # get per-cell containment & correctness
        gf  = self.compute_global_freq(train_labels)
        knn = KNeighborsClassifier(n_neighbors=k).fit(train_embeddings, train_labels)
        _, _, cont, corr, _ = self.adjusted_knn_predict_with_rank(
            knn, train_labels, test_embeddings, test_labels, gf, k
        )
        cont = np.array(cont)
        corr = np.array(corr)

        df     = self.adata.obs.copy()
        coords = self.adata.obsm["UMAP_embedding"]
        is_tr  = df[self.dataset_key] == "train"
        is_te  = df[self.dataset_key] == "test"

        # status: 0=train, 1=test‐correct+contained, 2=test‐incorrect/ not contained
        status = np.zeros(len(df), dtype=int)
        status[is_tr] = 0
        ti = np.where(is_te)[0]
        for idx, (cflag, flag) in zip(ti, zip(cont, corr)):
            status[idx] = 1 if (cflag==1 and flag==1) else 2

        cmap  = {0:"lightgray", 1:"green", 2:"red"}
        smap  = {0:12,         1:30,       2:30}
        mmap  = {0:".",        1:"x",      2:"x"}
        amap  = {0:0.2,        1:0.9,      2:0.9}

        fig, ax = plt.subplots(figsize=figsize, dpi=200)
        # draw train + incorrect behind
        for sv in [0,2]:
            mask = status == sv
            ax.scatter(coords[mask,0], coords[mask,1],
                       c=cmap[sv], marker=mmap[sv],
                       s=smap[sv], alpha=amap[sv],
                       label=("Train" if sv==0 else "Test Incorrect"))
        # correct on top
        mask = status == 1
        ax.scatter(coords[mask,0], coords[mask,1],
                   c=cmap[1], marker=mmap[1],
                   s=smap[1], alpha=amap[1],
                   label="Test Correct")

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title or "UMAP – Test‐cell Accuracy")

        handles, labels = ax.get_legend_handles_labels()
        by = dict(zip(labels, handles))
        ax.legend(by.values(), by.keys(),
                  bbox_to_anchor=(1.02,1), loc="upper left",
                  frameon=False, fontsize="small")

        plt.tight_layout()
        if savepath:
            savepath = os.path.splitext(savepath)[0] + ".png"
            fig.savefig(savepath, format="png", dpi=300, bbox_inches="tight")
        return fig, ax


    # ─── convenience wrapper ───────────────────────────────────────

    def run_all(
        self,
        train_embeddings,
        train_labels,
        test_embeddings,
        test_labels,
        k: int = 30
    ):
        stats = self.evaluate_adjusted_knn(
            train_embeddings, train_labels,
            test_embeddings,  test_labels,
            k
        )
        f1, a1 = self.plot_top_clones_umap()
        f2, a2 = self.plot_test_accuracy_umap(
            train_embeddings, train_labels,
            test_embeddings,  test_labels,
            k
        )
        return stats, (f1, a1), (f2, a2)
