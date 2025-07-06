import numpy as np
import umap
import seaborn as sns
import matplotlib.pyplot as plt

class LCL_Eval:
    """
    Given a concatenated AnnData (train+test), this class will:
      1) Read the high-dimensional LCL embedding from adata.obsm[embedding_key].
      2) Compute a 2D UMAP (and store it in adata.obsm["UMAP_embedding"]).
      3) Identify the top-N lineages by count.
      4) Plot:
         - "Other" cells (not in top-N) in light gray, low opacity
         - Top-N lineages in distinct colors, with train cells='.' and test cells='x'
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
              - adata.obs[clone_key] (lineage/clone ID for each cell)
              - adata.obs[dataset_key] (either "train" or "test")
              - adata.obsm[embedding_key] (high-dimensional LCL embedding, shape: n_cells × D)
        embedding_key : str
            Key under which the high-dimensional LCL embedding is stored in adata.obsm.
        clone_key : str
            Column name in adata.obs that holds each cell’s clone/lineage ID.
        dataset_key : str
            Column name in adata.obs that indicates "train" vs. "test."
        num_top : int
            How many top lineages (by frequency) to highlight. All others become “Other.”
        palette : list[str], optional
            A list of at least `num_top` color hex codes (or seaborn-compatible names). If None,
            uses seaborn’s default “tab10” palette.
        umap_kwargs : dict, optional
            Passed to `umap.UMAP(**umap_kwargs)` when computing 2D. Defaults to {"random_state": 42}.
        """
        # Copy AnnData so we don't modify the original in place
        self.adata = adata.copy()
        self.embedding_key = "LCL_embedding"
        self.clone_key = clone_key
        self.dataset_key = dataset_key
        self.num_top = num_top
        self.umap_kwargs = umap_kwargs or {"random_state": 42}

        # Use seaborn styling
        sns.set_style("whitegrid")

        # 1) Verify that the high-D LCL embedding exists:
        if self.embedding_key not in self.adata.obsm:
            raise KeyError(f"Could not find embedding under adata.obsm['{self.embedding_key}'].")

        hd = self.adata.obsm[self.embedding_key]
        if hd.ndim != 2 or hd.shape[1] < 3:
            raise ValueError(f"adata.obsm['{self.embedding_key}'] must be n_cells×D with D≥3.")

        # 2) Compute a 2D UMAP and store it in "UMAP_embedding"
        reducer = umap.UMAP(**self.umap_kwargs)
        umap2d = reducer.fit_transform(hd)  # (n_cells, 2)
        self.adata.obsm["UMAP_embedding"] = umap2d

        # 3) Find the top-N clones by frequency
        counts = self.adata.obs[self.clone_key].value_counts()
        self.top_clones = counts.index[: self.num_top].tolist()

        # 4) Create a new categorical column "clone_group": top-N clones stay as is, all else→"Other"
        def _group_fn(x):
            return x if x in self.top_clones else "Other"

        self.adata.obs["clone_group"] = (
            self.adata.obs[self.clone_key]
            .apply(_group_fn)
            .astype("category")
        )

        # 5) Reorder categories so that top_clones come first, then "Other"
        cats = self.top_clones + ["Other"]
        self.adata.obs["clone_group"] = (
            self.adata.obs["clone_group"]
                .cat
                .reorder_categories(cats, ordered=True)
        )

        # 6) Build a color map for top-N clones. If no palette provided, use seaborn tab10
        if palette is None:
            base_pal = sns.color_palette("tab10", n_colors=self.num_top)
            self.color_map = {self.top_clones[i]: base_pal[i] for i in range(self.num_top)}
        else:
            if len(palette) < self.num_top:
                raise ValueError(f"Palette has length {len(palette)} but num_top={self.num_top}.")
            self.color_map = {self.top_clones[i]: palette[i] for i in range(self.num_top)}

        # Always render "Other" in light gray
        self.color_map["Other"] = "lightgray"


    def plot_umap(self, figsize=(8,6), title=None, savepath: str = None):
        """
        Draws the UMAP scatterplot:
          - Plot "Other" cells first (train + test) in light gray, alpha=0.2
          - Then overlay each of the top-N clones:
                * Train cells: marker=".", size=30, alpha=0.8
                * Test  cells: marker="x", size=40, alpha=0.9
          - Deduplicate legend entries and place them outside.

        Returns:
        --------
        fig, ax : matplotlib Figure and Axes, in case you want further tweaks or re-saving.
        """
        df = self.adata.obs
        coords = self.adata.obsm["UMAP_embedding"]

        is_train = (df[self.dataset_key] == "train")
        is_test  = (df[self.dataset_key] == "test")

        fig, ax = plt.subplots(figsize=figsize, dpi=200)

        # 1) Draw "Other" (train & test) in light gray, low opacity
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

        # 2) Overlay each of the top-N clones
        for clone in self.top_clones:
            mask_train_clone = (df["clone_group"] == clone) & is_train
            mask_test_clone  = (df["clone_group"] == clone) & is_test

            # Train cells = dot
            ax.scatter(
                coords[mask_train_clone, 0],
                coords[mask_train_clone, 1],
                color=self.color_map[clone],
                s=30,
                marker=".",
                alpha=0.8,
                label=f"Train {clone}"
            )
            # Test cells = cross
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

        # Deduplicate legend entries
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