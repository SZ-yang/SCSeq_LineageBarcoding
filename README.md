Deep-learning embedding for static lineage barcoded single-cell RNA-seq data under the supervision of Prof. Kevin Lin at UW


|                      | Total cells | num of lineages (clones) | num of genes | Cell Types           | time pints |
| -------------------- | ----------- | ------------------------ | ------------ | -------------------- | ---------- |
| Larry Dataset        | 130887      | 5864                     | 25289        | 11(including undiff) | 2,4,6      |
| Larry Dataset for CL | 41201       | 2817                     | 2000         | 11(including undiff) | 2,4,6      |

Loss Function:

$$
loss_{m,n} = - \log \frac{\exp(\text{sim}(z_m, z_n)/\tau)}
{\sum ^{N}_{k=1} \mathbb{1}_{[k \neq i]}
\exp(\text{sim}(z_m, z_k)/\tau)}
$$


|                  | Lineage 1 Cell 1 | Lineage 2 Cell 1 | Lineage 1 Cell 2 | Lineage 2 Cell 2 |
| ---------------- | ---------------- | ---------------- | ---------------- | ---------------- |
| Lineage 1 Cell 1 |                  | Negative Pair    | Positive Pair    | Negative Pair    |
| Lineage 2 Cell 1 | Negative Pair    |                  | Negative Pair    | Positive Pair    |
| Lineage 1 Cell 2 | Positive Pair    | Negative Pair    |                  | Negative Pair    |
| Lineage 2 Cell 2 | Negative Pair    | Positive Pair    | Negative Pair    |                  |