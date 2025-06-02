import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import pairwise_distances, accuracy_score
from typing import Tuple, List, Dict

# Utility function to calculate global frequencies
def compute_global_freq(labels: np.ndarray) -> Dict[int, float]:
    unique, counts = np.unique(labels, return_counts=True)
    total = len(labels)
    return {lab: cnt / total for lab, cnt in zip(unique, counts)}

# Modified KNN classifier (sklearn-based)
def adjusted_knn_predict_with_rank(
    knn_model: KNeighborsClassifier,
    train_labels: np.ndarray,
    embeddings: np.ndarray,
    labels_true: np.ndarray,
    global_freq: Dict[int, float],
    k: int
) -> Tuple[np.ndarray, List[int], float, float, float]:
    neigh_indices = knn_model.kneighbors(embeddings, return_distance=False)

    preds = []
    ranks = []
    unique_label_counts = []

    for idx, nbrs in enumerate(neigh_indices):
        nbr_labels = train_labels[nbrs]
        unique, counts = np.unique(nbr_labels, return_counts=True)
        num_unique_labels = len(unique)
        unique_label_counts.append(num_unique_labels)

        local_freq = {lab: cnt / k for lab, cnt in zip(unique, counts)}

        scores = {lab: local_freq.get(lab, 0.0) - global_freq.get(lab, 0.0) for lab in global_freq}

        pred = max(scores.items(), key=lambda x: x[1])[0]
        preds.append(pred)

        sorted_labels = sorted(local_freq.items(), key=lambda x: x[1], reverse=True)
        label_ranks = [lab for lab, _ in sorted_labels]

        true_label = labels_true[idx]
        if true_label in label_ranks:
            rank = label_ranks.index(true_label) + 1
        else:
            rank = num_unique_labels + 1

        ranks.append(rank)

    avg_rank = np.mean(ranks)
    avg_unique_labels = np.mean(unique_label_counts)
    accuracy = accuracy_score(labels_true, preds)

    return np.array(preds), ranks, avg_rank, avg_unique_labels, accuracy

# Modified KNN classifier (distance-based)
def adjusted_knn_predict_with_rank_v2(
    train_embeddings: np.ndarray,
    train_labels: np.ndarray,
    embeddings: np.ndarray,
    labels_true: np.ndarray,
    global_freq: Dict[int, float],
    k: int = 30
) -> Tuple[np.ndarray, List[int], float, float, float]:
    distances = pairwise_distances(embeddings, train_embeddings)

    predictions = []
    ranks = []
    unique_label_counts = []

    for idx, dist in enumerate(distances):
        nearest_indices = np.argsort(dist)[:k]
        neighbor_labels = train_labels[nearest_indices]

        unique_labels, counts = np.unique(neighbor_labels, return_counts=True)
        num_unique_labels = len(unique_labels)
        unique_label_counts.append(num_unique_labels)

        local_freq = {lab: cnt / k for lab, cnt in zip(unique_labels, counts)}

        adjusted_scores = {lab: local_freq.get(lab, 0) - global_freq.get(lab, 0)
                           for lab in global_freq}

        predicted_label = max(adjusted_scores, key=adjusted_scores.get)
        predictions.append(predicted_label)

        sorted_labels = sorted(local_freq.items(), key=lambda x: x[1], reverse=True)
        label_ranks = [lab for lab, _ in sorted_labels]

        true_label = labels_true[idx]
        if true_label in label_ranks:
            rank = label_ranks.index(true_label) + 1
        else:
            rank = num_unique_labels + 1

        ranks.append(rank)

    avg_rank = np.mean(ranks)
    avg_unique_labels = np.mean(unique_label_counts)
    accuracy = accuracy_score(labels_true, predictions)

    return np.array(predictions), ranks, avg_rank, avg_unique_labels, accuracy

# Wrapper to easily calculate both train and test accuracy for adjusted KNN models
def evaluate_adjusted_knn(
    train_embeddings: np.ndarray,
    train_labels: np.ndarray,
    test_embeddings: np.ndarray,
    test_labels: np.ndarray,
    k: int = 30
):
    global_freq = compute_global_freq(train_labels)

    knn_model = KNeighborsClassifier(n_neighbors=k).fit(train_embeddings, train_labels)

    # Train results
    _, train_ranks, train_avg_rank, train_avg_unique, train_accuracy = adjusted_knn_predict_with_rank(
        knn_model, train_labels, train_embeddings, train_labels, global_freq, k)

    # Test results
    _, test_ranks, test_avg_rank, test_avg_unique, test_accuracy = adjusted_knn_predict_with_rank(
        knn_model, train_labels, test_embeddings, test_labels, global_freq, k)

    # Distance-based version results (train)
    _, train_ranks_v2, train_avg_rank_v2, train_avg_unique_v2, train_accuracy_v2 = adjusted_knn_predict_with_rank_v2(
        train_embeddings, train_labels, train_embeddings, train_labels, global_freq, k)

    # Distance-based version results (test)
    _, test_ranks_v2, test_avg_rank_v2, test_avg_unique_v2, test_accuracy_v2 = adjusted_knn_predict_with_rank_v2(
        train_embeddings, train_labels, test_embeddings, test_labels, global_freq, k)

    return {
        "train": {"accuracy": train_accuracy, "avg_rank": train_avg_rank, "avg_unique_labels": train_avg_unique},
        "test": {"accuracy": test_accuracy, "avg_rank": test_avg_rank, "avg_unique_labels": test_avg_unique},
        "train_v2": {"accuracy": train_accuracy_v2, "avg_rank": train_avg_rank_v2, "avg_unique_labels": train_avg_unique_v2},
        "test_v2": {"accuracy": test_accuracy_v2, "avg_rank": test_avg_rank_v2, "avg_unique_labels": test_avg_unique_v2}
    }