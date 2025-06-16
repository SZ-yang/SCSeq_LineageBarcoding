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
      1) compute_adjusted_knn_stats(...) → returns train/test KNN accuracy + rank‐based + unique‐label stats
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
        # copy so we don’t modify original
        self.adata = adata.copy()
        self.embedding_key = "LCL_embedding"
        self.clone_key = clone_key
        self.dataset_key = dataset_key
        self.num_top = num_top
        self.umap_kwargs = umap_kwargs or {"random_state": 42}

        sns.set_style("whitegrid")

        # verify high-D embedding
        if self.embedding_key not in self.adata.obsm:
            raise KeyError(f"Could not find embedding under adata.obsm['{self.embedding_key}'].")

        hd = self.adata.obsm[self.embedding_key]
        if hd.ndim != 2 or hd.shape[1] < 3:
            raise ValueError(f"adata.obsm['{self.embedding_key}'] must be shape [n_cells, D] with D≥3.")

        # compute 2D UMAP
        reducer = umap.UMAP(**self.umap_kwargs)
        umap2d = reducer.fit_transform(hd)
        self.adata.obsm["UMAP_embedding"] = umap2d

        # top-N clones
        counts = self.adata.obs[self.clone_key].value_counts()
        self.top_clones = counts.index[:self.num_top].tolist()

        # clone_group column
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

        # color map
        if palette is None:
            base_pal = sns.color_palette("tab10", n_colors=self.num_top)
            self.color_map = {self.top_clones[i]: base_pal[i] for i in range(self.num_top)}
        else:
            if len(palette) < self.num_top:
                raise ValueError(f"Palette length {len(palette)} < num_top={self.num_top}")
            self.color_map = {self.top_clones[i]: palette[i] for i in range(self.num_top)}

        self.color_map["Other"] = "lightgray"


    # ─── KNN & rank/unique evaluation ───────────────────────────

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
         np.ndarray, List[int], float, float, List[int], List[int]
    ]:
        """
        Returns preds, ranks, avg_rank, avg_unique_lbl, correct_flags, unique_label_counts.
        ⚠️ If true_label not among neighbors, rank = k+1
        """
        neigh = knn_model.kneighbors(embeddings, return_distance=False)
        n = embeddings.shape[0]

        preds = np.zeros(n, dtype=train_labels.dtype)
        ranks = []
        unique_label_counts = []
        correct_flags = []

        for i, nbrs in enumerate(neigh):
            nbr_lab = train_labels[nbrs]
            uq, ct = np.unique(nbr_lab, return_counts=True)
            num_unique = len(uq)
            unique_label_counts.append(num_unique)

            # local freq
            local_freq = {lab: c/k for lab, c in zip(uq, ct)}

            # adjusted scores
            scores = {lab: local_freq.get(lab,0)-global_freq.get(lab,0)
                      for lab in global_freq.keys()}

            best = max(scores.items(), key=lambda x: x[1])[0]
            preds[i] = best

            # rank of true
            sorted_lbl = sorted(local_freq.items(), key=lambda x: x[1], reverse=True)
            ranked = [lab for lab,_ in sorted_lbl]
            t = labels_true[i]
            if t in ranked:
                rpos = ranked.index(t)+1
            else:
                rpos = k+1    
            ranks.append(rpos)

            correct_flags.append(1 if best==t else 0)

        avg_rank = float(np.mean(ranks))
        avg_unique = float(np.mean(unique_label_counts))

        return preds, ranks, avg_rank, avg_unique, correct_flags, unique_label_counts


    def evaluate_adjusted_knn(
        self,
        train_embeddings: np.ndarray,
        train_labels: np.ndarray,
        test_embeddings: np.ndarray,
        test_labels: np.ndarray,
        k: int = 30
    ) -> Dict[str, Dict]:
        global_freq = self.compute_global_freq(train_labels)

        knn = KNeighborsClassifier(n_neighbors=k)
        knn.fit(train_embeddings, train_labels)

        # train side
        ( _,
          train_ranks,
          _train_avg_rank_layer1,
          _train_avg_unique_layer1,
          train_flags,
          train_uniq_counts
        ) = self.adjusted_knn_predict_with_rank(
                knn, train_labels, train_embeddings, train_labels, global_freq, k
            )

        # test side
        ( _,
          test_ranks,
          _test_avg_rank_layer1,
          _test_avg_unique_layer1,
          test_flags,
          test_uniq_counts
        ) = self.adjusted_knn_predict_with_rank(
                knn, train_labels, test_embeddings, test_labels, global_freq, k
            )

        # two-layer aggregate function
        def layer2_aggregate(labels, flags, ranks, uniqs):
            per_acc, per_rank, per_uniq = {}, {}, {}
            ulabels = np.unique(labels)
            for lab in ulabels:
                idx = np.where(labels==lab)[0]
                per_acc[lab]  = float(np.mean([flags[i] for i in idx]))
                per_rank[lab] = float(np.mean([ranks[i] for i in idx]))
                per_uniq[lab] = float(np.mean([uniqs[i] for i in idx]))
            overall_acc  = float(np.mean(list(per_acc.values())))
            overall_rank = float(np.mean(list(per_rank.values())))
            overall_uniq = float(np.mean(list(per_uniq.values())))
            lineage_ranks = [per_rank[lab] for lab in sorted(ulabels)]
            return per_acc, per_rank, per_uniq, overall_acc, overall_rank, overall_uniq, lineage_ranks

        # aggregate
        (train_pa, train_pr, train_pu,
         train_oa, train_or, train_ou, train_lr) = layer2_aggregate(
            train_labels, train_flags, train_ranks, train_uniq_counts
        )
        (test_pa, test_pr, test_pu,
         test_oa, test_or, test_ou, test_lr) = layer2_aggregate(
            test_labels, test_flags, test_ranks, test_uniq_counts
        )

        # quantiles
        def qtls(arr):
            arr = np.array(arr)
            return {
              "q25": round(float(np.quantile(arr,0.25)),3),
              "q50": round(float(np.quantile(arr,0.50)),3),
              "q75": round(float(np.quantile(arr,0.75)),3)
            }

        return {
          "train": {
            "overall_accuracy":  round(train_oa,3),
            "overall_avg_rank":  round(train_or,3),
            "overall_avg_unique":round(train_ou,3),
            "rank_quantiles":    qtls(train_lr)
          },
          "test": {
            "overall_accuracy":  round(test_oa,3),
            "overall_avg_rank":  round(test_or,3),
            "overall_avg_unique":round(test_ou,3),
            "rank_quantiles":    qtls(test_lr)
          }
        }


    # ─── UMAP PLOTS ────────────────────────────────────────────────────

    def plot_top_clones_umap(
        self,
        figsize=(8,6),
        title=None,
        savepath: str = None
    ):
        df = self.adata.obs
        coords = self.adata.obsm["UMAP_embedding"]
        is_tr = (df[self.dataset_key]=="train")
        is_te = (df[self.dataset_key]=="test")

        fig,ax = plt.subplots(figsize=figsize,dpi=200)

        # Other first
        m1 = (df["clone_group"]=="Other") & is_tr
        m2 = (df["clone_group"]=="Other") & is_te

        ax.scatter(coords[m1,0],coords[m1,1],c=self.color_map["Other"],
                   s=8,marker=".",alpha=0.2,label="Train Other")
        ax.scatter(coords[m2,0],coords[m2,1],c=self.color_map["Other"],
                   s=12,marker="x",alpha=0.2,label="Test Other")

        # top clones
        for c in self.top_clones:
            mt = (df["clone_group"]==c)&is_tr
            mx = (df["clone_group"]==c)&is_te
            ax.scatter(coords[mt,0],coords[mt,1],c=self.color_map[c],
                       s=30,marker=".",alpha=0.8,label=f"Train {c}")
            ax.scatter(coords[mx,0],coords[mx,1],c=self.color_map[c],
                       s=40,marker="x",alpha=0.9,label=f"Test {c}")

        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_title(title or f"UMAP – Top {self.num_top} Clones")

        h,l = ax.get_legend_handles_labels()
        by = dict(zip(l,h))
        ax.legend(by.values(),by.keys(),
                  bbox_to_anchor=(1.02,1),loc="upper left",
                  frameon=False,fontsize="small")

        plt.tight_layout()
        if savepath:
            if not savepath.lower().endswith(".png"):
                savepath = os.path.splitext(savepath)[0]+".png"
            fig.savefig(savepath,format="png",dpi=300,bbox_inches="tight")
        return fig,ax


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
        # get per-cell correctness
        gf = self.compute_global_freq(train_labels)
        knn = KNeighborsClassifier(n_neighbors=k).fit(train_embeddings, train_labels)
        _,_,_,_, flags, _ = self.adjusted_knn_predict_with_rank(
            knn, train_labels, test_embeddings, test_labels, gf, k
        )
        flags = np.array(flags)

        df = self.adata.obs.copy()
        coords = self.adata.obsm["UMAP_embedding"]
        is_tr = (df[self.dataset_key]=="train")
        is_te = (df[self.dataset_key]=="test")

        status = np.zeros(len(df),dtype=int)
        status[is_tr] = 0
        ti = np.where(is_te)[0]
        for idx,fl in zip(ti, flags):
            status[idx] = 1 if fl==1 else 2

        cmap = {0:"lightgray",1:"green",2:"red"}
        smap = {0:12,1:30,2:30}
        mmap = {0:".",1:"x",2:"x"}
        amap = {0:0.2,1:0.9,2:0.9}

        fig,ax = plt.subplots(figsize=figsize,dpi=200)
        # draw train+incorrect first
        for sv in [0,2]:
            m = (status==sv)
            ax.scatter(coords[m,0],coords[m,1],
                       c=cmap[sv],s=smap[sv],
                       marker=mmap[sv],alpha=amap[sv],
                       label=("Train" if sv==0 else "Test Incorrect"))
        # correct on top
        m = (status==1)
        ax.scatter(coords[m,0],coords[m,1],
                   c=cmap[1],s=smap[1],
                   marker=mmap[1],alpha=amap[1],
                   label="Test Correct")

        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_title(title or "UMAP – Test‐cell Accuracy")

        h,l = ax.get_legend_handles_labels()
        by = dict(zip(l,h))
        ax.legend(by.values(),by.keys(),
                  bbox_to_anchor=(1.02,1),loc="upper left",
                  frameon=False,fontsize="small")

        plt.tight_layout()
        if savepath:
            if not savepath.lower().endswith(".png"):
                savepath = os.path.splitext(savepath)[0]+".png"
            fig.savefig(savepath,format="png",dpi=300,bbox_inches="tight")
        return fig,ax


    # ─── convenience wrapper ───────────────────────────────────────

    def run_all(
        self,
        train_embeddings, train_labels,
        test_embeddings,  test_labels,
        k: int = 30
    ):
        stats = self.evaluate_adjusted_knn(
            train_embeddings, train_labels,
            test_embeddings,  test_labels,
            k
        )
        f1,a1 = self.plot_top_clones_umap()
        f2,a2 = self.plot_test_accuracy_umap(
            train_embeddings, train_labels,
            test_embeddings,  test_labels,
            k
        )
        return stats, (f1,a1), (f2,a2)