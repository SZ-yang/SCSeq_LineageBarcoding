# LCL_eval_modifiedKNN.py

import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from typing import Tuple, List, Dict

# Utility function to calculate global frequencies of labels in an array
def compute_global_freq(labels: np.ndarray) -> Dict[int, float]:
    unique, counts = np.unique(labels, return_counts=True)
    total = len(labels)
    return {lab: cnt / total for lab, cnt in zip(unique, counts)}

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
      1) find its k nearest TRAINING neighbors
      2) compute local_freq[L] = (# neighbors with label L) / k
      3) adjusted_score[L] = local_freq[L] - global_freq[L]
      4) predict = argmax_{L}(adjusted_score[L])
      5) compute rank_of_true_label = rank position of labels_true[i] in descending local_freq order
         (if true_label not among neighbors, rank = num_unique_labels + 1)

    Returns:
      preds          : np.ndarray, shape (n_points,)
      ranks          : List[int],   length n_points, each the rank of the true label among that point’s neighbors
      avg_rank       : float = mean(ranks)
      avg_unique_lbl : float = average number of unique labels seen among each point’s k neighbors
      correct_flags  : List[int], length n_points, 1 if predicted == true_label else 0
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

        # build adjusted scores
        scores = {lab: local_freq.get(lab, 0.0) - global_freq.get(lab, 0.0)
                  for lab in global_freq.keys()}

        # predict label = argmax adjusted score
        best_label = max(scores.items(), key=lambda x: x[1])[0]
        preds[i] = best_label

        # compute the rank of the true label among local_freq
        sorted_labels = sorted(local_freq.items(), key=lambda x: x[1], reverse=True)
        ranked_labs = [lab for lab, _ in sorted_labels]
        true_lab = labels_true[i]
        if true_lab in ranked_labs:
            rank_pos = ranked_labs.index(true_lab) + 1  # 1-based
        else:
            rank_pos = num_unique + 1
        ranks.append(rank_pos)

        # correct? 1 if pred==true, else 0
        correct_flags.append(1 if best_label == true_lab else 0)

    avg_rank = float(np.mean(ranks))
    avg_unique_lbl = float(np.mean(unique_label_counts))

    return preds, ranks, avg_rank, avg_unique_lbl, correct_flags

def evaluate_adjusted_knn(
    train_embeddings: np.ndarray,
    train_labels: np.ndarray,
    test_embeddings: np.ndarray,
    test_labels: np.ndarray,
    k: int = 30
) -> Dict[str, Dict]:
    """
    Wrapper that fits an adjusted‐KNN on (train_embeddings, train_labels), then evaluates both on train & test.

    Internally:
      - Computes global_freq from train_labels
      - Fits a KNN (n_neighbors=k) on train set
      - Calls adjusted_knn_predict_with_rank(...) on both the train_set itself and the test_set.
      - Then reorganizes the raw per‐cell results into two layers of averaging:
          Layer 1: average metrics over cells *within the same lineage*
          Layer 2: average the Layer‐1 quantities across all lineages (so each lineage counts equally).
      - Also collects each lineage’s “average rank” so you can compute quantiles over that vector of per‐lineage averages.

    Returns a dict of the form:
    {
      "train": {
          "overall_accuracy"     : float,  # Layer‐2 average of per‐lineage accuracies
          "overall_avg_rank"     : float,  # Layer‐2 average of per‐lineage average ranks
          "overall_avg_unique"   : float,  # Layer‐2 average of per‐lineage avg_unique_labels
          "per_lineage_accuracy" : Dict[label, float],  # per‐lineage accuracies
          "per_lineage_avg_rank" : Dict[label, float],  # per‐lineage average ranks
          "rank_quantiles"       : Dict[qname, float]   # e.g. {"q25":…, "q50":…, "q75":…}
      },
      "test": { … same keys … }
    }
    """
    # 1) Global frequencies come from training labels only
    global_freq = compute_global_freq(train_labels)

    # 2) Fit standard KNN on train_embeddings
    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(train_embeddings, train_labels)

    # 3) Run adjusted predictions & collect per‐cell metrics for train and test
    # 3a) Train side
    (
      train_preds,
      train_ranks,
      train_avg_rank,
      train_avg_unique,
      train_correct_flags
    ) = adjusted_knn_predict_with_rank(
            knn_model=knn,
            train_labels=train_labels,
            embeddings=train_embeddings,
            labels_true=train_labels,
            global_freq=global_freq,
            k=k
        )

    # 3b) Test side
    (
      test_preds,
      test_ranks,
      test_avg_rank,
      test_avg_unique,
      test_correct_flags
    ) = adjusted_knn_predict_with_rank(
            knn_model=knn,
            train_labels=train_labels,
            embeddings=test_embeddings,
            labels_true=test_labels,
            global_freq=global_freq,
            k=k
        )

    # 4) Now we have for each cell in train/test:
    #    - correct_flags[i] (0/1)
    #    - ranks[i]
    #    - unique_label_count[i]  (not returned here but the function computed it)
    # We also need each cell’s lineage (label) to aggregate within each lineage.

    # Build a helper to do the two‐layer averaging:
    def layer2_aggregate(cell_labels: np.ndarray,
                         correct_flags: List[int],
                         ranks: List[int]
                       ) -> Tuple[
                              Dict[int,float],    # per-lineage accuracy
                              Dict[int,float],    # per-lineage avg_rank
                              float,              # two‐layer overall accuracy
                              float,              # two‐layer overall avg_rank
                              List[float]         # vector of per-lineage avg_rank (in label‐order)
                         ]:
        """
        Given arrays of per‐cell flags and per‐cell ranks, plus each cell’s lineage label,
        produce:
          1) per_lineage_accuracy: { lineage_label -> mean(correct_flags of cells in that lineage) }
          2) per_lineage_avg_rank: { lineage_label -> mean(ranks of cells in that lineage) }
          3) overall_accuracy:   mean over all lineages of (per_lineage_accuracy[l])
          4) overall_avg_rank:  mean over all lineages of (per_lineage_avg_rank[l])
          5) lineage_rank_list: a list of per-lineage average ranks in the same label‐order
        """
        per_lineage_accuracy = {}
        per_lineage_avg_rank = {}

        unique_lineages = np.unique(cell_labels)
        for lab in unique_lineages:
            idxs = np.where(cell_labels == lab)[0]
            lab_acc = np.mean([correct_flags[i] for i in idxs])
            lab_avg_rank = np.mean([ranks[i] for i in idxs])
            per_lineage_accuracy[lab] = float(lab_acc)
            per_lineage_avg_rank[lab] = float(lab_avg_rank)

        # Two‐layer overall averages (treat each lineage equally)
        overall_acc = float(np.mean(list(per_lineage_accuracy.values())))
        overall_rank = float(np.mean(list(per_lineage_avg_rank.values())))

        # Build a list of per‐lineage average ranks in sorted‐label order
        lineage_rank_list = [per_lineage_avg_rank[lab] for lab in sorted(unique_lineages)]

        return per_lineage_accuracy, per_lineage_avg_rank, overall_acc, overall_rank, lineage_rank_list

    # 5) Aggregate train‐side
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

    # 6) Aggregate test‐side
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

    # 7) Compute quantiles over the vector of per‐lineage average ranks (for train and test)
    #   e.g. 25%, 50%, 75%
    def compute_quantiles(rank_list: List[float]) -> Dict[str, float]:
        arr = np.array(rank_list)
        q25 = round(float(np.quantile(arr, 0.25)), 3)
        q50 = round(float(np.quantile(arr, 0.50)), 3)
        q75 = round(float(np.quantile(arr, 0.75)), 3)
        return {"q25": q25, "q50": q50, "q75": q75}

    train_quantiles = compute_quantiles(train_lineage_rank_list)
    test_quantiles  = compute_quantiles(test_lineage_rank_list)

    # 8) Package everything into a single return dict
    results = {
      "train": {
        "overall_accuracy":        train_overall_acc,
        "overall_avg_rank":        train_overall_rank,
        "per_lineage_accuracy":    train_per_lin_acc,
        "per_lineage_avg_rank":    train_per_lin_rank,
        "rank_quantiles":          train_quantiles
      },
      "test": {
        "overall_accuracy":        test_overall_acc,
        "overall_avg_rank":        test_overall_rank,
        "per_lineage_accuracy":    test_per_lin_acc,
        "per_lineage_avg_rank":    test_per_lin_rank,
        "rank_quantiles":          test_quantiles
      }
    }

    return results