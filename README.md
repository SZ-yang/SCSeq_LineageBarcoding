# Deep-learning embedding for static lineage barcoded single-cell RNA-seq data
Under the supervision of Prof. Kevin Lin at the University of Washington

## Some Statistics of the Dataset

| Dataset                | Total cells | Num of lineages (clones) | Num of genes | Cell Types           | Time points |
|------------------------|-------------|--------------------------|--------------|----------------------|-------------|
| Larry Dataset          | 130,887     | 5,864                    | 25,289       | 11 (including undiff)| 2, 4, 6     |
| Larry Dataset for CL   | 41,201      | 2,817                    | 2,000        | 11 (including undiff)| 2, 4, 6     |

## Loss Function

The loss function is defined as:

$$
loss_{m,n} = - \log \left( \frac{\exp(\text{sim}(z_m, z_n)/\tau)}{\sum ^{N}_{k=1} \mathbb{1}_{[k \neq i]} \exp(\text{sim}(z_m, z_k)/\tau)} \right)
$$

Where:
- \( \text{sim}(z_m, z_n) \) is the similarity between the embeddings of cells \( m \) and \( n \),
- \( \tau \) is the temperature parameter that scales the similarity,
- \( N \) is the total number of cells in the dataset,
- \( \mathbb{1}_{[k \neq i]} \) is the indicator function that is 1 if \( k \neq i \) and 0 otherwise.

## Pairwise Comparisons Between Lineages

|                  | Lineage 1 Cell 1 | Lineage 2 Cell 1 | Lineage 1 Cell 2 | Lineage 2 Cell 2 |
|------------------|------------------|------------------|------------------|------------------|
| Lineage 1 Cell 1 |                  | Negative Pair    | Positive Pair    | Negative Pair    |
| Lineage 2 Cell 1 | Negative Pair    |                  | Negative Pair    | Positive Pair    |
| Lineage 1 Cell 2 | Positive Pair    | Negative Pair    |                  | Negative Pair    |
| Lineage 2 Cell 2 | Negative Pair    | Positive Pair    | Negative Pair    |                  |
