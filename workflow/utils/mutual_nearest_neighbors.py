from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.optimize import linear_sum_assignment
import numpy as np
import pandas as pd

def compute_weight(rna_rank, atac_rank, k):
    rank_sum = rna_rank + atac_rank
    min_val, max_val = 2, 2 * (k + 1)
    return np.clip((max_val - rank_sum) / (max_val - min_val), 0.0, 1.0)


def mutual_nearest_neighbors(rna_embeddings, atac_embeddings, rna_names, atac_names, k=50):
    """
    Perform mutual nearest neighbors (MNN) matching between RNA and ATAC embeddings.

    This function aligns cells from two modalities (RNA and ATAC) by:
    - Reducing the embeddings with PCA.
    - Computing the nearest neighbors from RNA to ATAC and vice versa.
    - Constructing rank matrices based on neighbor relationships.
    - Computing a weight matrix from these ranks.
    - Applying the Hungarian algorithm to find an optimal matching between cells.

    Parameters
    ----------
    rna_embeddings : np.ndarray
        A 2D array of shape (n_rna_cells, n_features) representing RNA embeddings.
    atac_embeddings : np.ndarray
        A 2D array of shape (n_atac_cells, n_features) representing ATAC embeddings.
    rna_names : array-like
        A list or 1D array of RNA cell names with length matching the number of rows in `rna_embeddings`.
    atac_names : array-like
        A list or 1D array of ATAC cell names with length matching the number of rows in `atac_embeddings`.
    k : int, optional (default: 20)
        The number of nearest neighbors to consider.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns:
        - 'RNA': the RNA cell names
        - 'ATAC': the matched ATAC cell names
        - 'Weight': the matching weight between the cells
    """

    n_pcs = min(rna_embeddings.shape[0], rna_embeddings.shape[1], 50)
    pca = PCA(n_components=n_pcs)
    all_embeddings = np.vstack([rna_embeddings, atac_embeddings])
    pca_result = pca.fit_transform(all_embeddings)
    rna_pca = pca_result[:rna_embeddings.shape[0]]
    atac_pca = pca_result[rna_embeddings.shape[0]:]

    rna_nn = NearestNeighbors(n_neighbors=k).fit(atac_pca)
    atac_nn = NearestNeighbors(n_neighbors=k).fit(rna_pca)
    _, rna_to_atac_idx = rna_nn.kneighbors(rna_pca)
    _, atac_to_rna_idx = atac_nn.kneighbors(atac_pca)

    n_rna, n_atac = len(rna_names), len(atac_names)

    rna_rank_matrix = np.full((n_rna, n_atac), k + 1)
    atac_rank_matrix = np.full((n_rna, n_atac), k + 1)

    for i in range(n_rna):
        for rank, atac_idx in enumerate(rna_to_atac_idx[i]):
            rna_rank_matrix[i, atac_idx] = rank + 1
    for j in range(n_atac):
        for rank, rna_idx in enumerate(atac_to_rna_idx[j]):
            atac_rank_matrix[rna_idx, j] = rank + 1

    weight_matrix = compute_weight(rna_rank_matrix, atac_rank_matrix, k=k)
    cost_matrix = -weight_matrix
    rna_idx, atac_idx = linear_sum_assignment(cost_matrix)

    return pd.DataFrame({
        "RNA": np.array(rna_names)[rna_idx],
        "ATAC": np.array(atac_names)[atac_idx],
        "Weight": weight_matrix[rna_idx, atac_idx],
        "RNA_embedding": list(rna_embeddings[rna_idx]),
        "ATAC_embedding": list(atac_embeddings[atac_idx])
    })
