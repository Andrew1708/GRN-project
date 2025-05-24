import os
import sys
import gc
import torch
import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.optimize import linear_sum_assignment

# === Local imports ===

# === Configuration ===
DEVICE = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
SEEDS = [0]
LATENT_DIM = 4
LEARNING_RATE = 3e-4
N_EPOCHS = 500
USE_ANCHOR = False
REG_D = REG_G = REG_MMD = 1
TS = [30, 50, 70]
USE_POTENTIAL = True

def compute_weight(rna_rank, atac_rank, k):
    rank_sum = rna_rank + atac_rank
    min_val, max_val = 2, 2 * (k + 1)
    return np.clip((max_val - rank_sum) / (max_val - min_val), 0.0, 1.0)


def mutual_nearest_neighbors(rna_embeddings, atac_embeddings, rna_counts, atac_counts, k=20):
    n_pcs = min(rna_embeddings.shape[0], rna_embeddings.shape[1])
    pca = PCA(n_components=n_pcs)
    all_embeddings = np.vstack([rna_embeddings, atac_embeddings])
    pca_result = pca.fit_transform(all_embeddings)
    rna_pca = pca_result[:rna_embeddings.shape[0]]
    atac_pca = pca_result[rna_embeddings.shape[0]:]

    rna_nn = NearestNeighbors(n_neighbors=k).fit(atac_pca)
    atac_nn = NearestNeighbors(n_neighbors=k).fit(rna_pca)
    _, rna_to_atac_idx = rna_nn.kneighbors(rna_pca)
    _, atac_to_rna_idx = atac_nn.kneighbors(atac_pca)

    rna_names = rna_counts.index.to_numpy()
    atac_names = atac_counts.index.to_numpy()
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
        "RNA": rna_names[rna_idx],
        "ATAC": atac_names[atac_idx],
        "Weight": weight_matrix[rna_idx, atac_idx]
    })


def process_project(project_name, rna_path, atac_path, reg_path, out_dir):
    import scDART.scDART.utils as utils
    import scDART.scDART.TI as ti
    from scDART.scDART import scDART

    print(f"\n=== Processing project: {project_name} ===")
    # === Load data ===
    counts_rna = pd.read_csv(rna_path, index_col=0)
    counts_atac = pd.read_csv(atac_path, index_col=0)
    coarse_reg = pd.read_csv(reg_path, index_col=0)

    # === Run scDART ===
    model = scDART(
        n_epochs=N_EPOCHS,
        latent_dim=LATENT_DIM,
        ts=TS,
        use_anchor=USE_ANCHOR,
        use_potential=USE_POTENTIAL,
        k=10,
        reg_d=REG_D,
        reg_g=REG_G,
        reg_mmd=REG_MMD,
        l_dist_type='kl',
        seed=SEEDS[0],
        device=DEVICE
    )

    model = model.fit(
        rna_count=counts_rna.values,
        atac_count=counts_atac.values,
        reg=coarse_reg.values,
        rna_anchor=None,
        atac_anchor=None
    )

    z_rna, z_atac = model.transform(
        rna_count=counts_rna.values,
        atac_count=counts_atac.values,
        rna_anchor=None,
        atac_anchor=None
    )

    # === Match cells via MNN ===
    matches = mutual_nearest_neighbors(z_rna, z_atac, counts_rna, counts_atac)
    matches.to_csv(os.path.join(out_dir, f"{project_name}_matches.csv"), index=False)
    print(f"✅ Saved matches for {project_name}")

    # === Cleanup ===
    del model, z_rna, z_atac, counts_rna, counts_atac, coarse_reg
    gc.collect()
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    torch.cuda.reset_accumulated_memory_stats()
    torch.cuda.synchronize()


def parse_args():
    parser = argparse.ArgumentParser(description="Run scDART + MNN matching for a single project.")
    parser.add_argument("--project_name", required=True, help="Project name (used for naming outputs)")
    parser.add_argument("--rna", required=True, help="Path to RNA counts CSV")
    parser.add_argument("--atac", required=True, help="Path to ATAC counts CSV")
    parser.add_argument("--reg", required=True, help="Path to region2gene CSV")
    parser.add_argument("--out_dir", required=True, help="Output directory for embeddings and matches")
    parser.add_argument("--scdart_dir", required=True, help="Path to scDART directory (for imports)")

    return parser.parse_args()


def main():
    args = parse_args()

    sys.path.insert(0, args.scdart_dir)

    process_project(
        project_name=args.project_name,
        rna_path=args.rna,
        atac_path=args.atac,
        reg_path=args.reg,
        out_dir=args.out_dir
    )


if __name__ == "__main__":
    main()
