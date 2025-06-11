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
from pathlib import Path

# === Local imports ===
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.mutual_nearest_neighbors import *

# === Configuration ===
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
SEEDS = [0]
LATENT_DIM = 4
LEARNING_RATE = 3e-4
N_EPOCHS = 500
USE_ANCHOR = False
REG_D = REG_G = REG_MMD = 1
TS = [30, 50, 70]
USE_POTENTIAL = True

def process_project(project_name, rna_path, atac_path, reg_path, temp_dir, out_dir):
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

    # === Save embeddings ===
    np.save(os.path.join(temp_dir,"z_rna.npy"), z_rna)
    np.save(os.path.join(temp_dir,"z_atac.npy"), z_atac)

    # === Match cells via MNN ===
    rna_names = counts_rna.index.to_numpy()
    atac_names = counts_atac.index.to_numpy()

    matches = mutual_nearest_neighbors(z_rna, z_atac, rna_names, atac_names)
    matches.to_csv(os.path.join(out_dir, "matches.csv"), index=False)
    print(f"✅ Saved matches for {project_name}")


def parse_args():
    parser = argparse.ArgumentParser(description="Run scDART + MNN matching for a single project.")
    parser.add_argument("--project_name", required=True, help="Project name (used for naming outputs)")
    parser.add_argument("--rna", required=True, help="Path to RNA counts CSV")
    parser.add_argument("--atac", required=True, help="Path to ATAC counts CSV")
    parser.add_argument("--reg", required=True, help="Path to region2gene CSV")
    parser.add_argument("--temp_dir", help="Temporary directory for intermediate files")
    parser.add_argument("--out_dir", required=True, help="Output directory for embeddings and matches")
    parser.add_argument("--scdart_dir", required=True, help="Path to scDART directory (for imports)")
    parser.add_argument("--cuda", type=int, default=0, help="CUDA device index (default: 0)")

    return parser.parse_args()


def main():
    args = parse_args()

    sys.path.insert(0, args.scdart_dir)

    cuda = f"cuda:{args.cuda}" if args.cuda != "cpu" else "cpu"
    if cuda != "cpu":
        DEVICE = torch.device(f"cuda:{args.cuda}" if torch.cuda.is_available() else "cpu")
    
    process_project(
        project_name=args.project_name,
        rna_path=args.rna,
        atac_path=args.atac,
        reg_path=args.reg,
        temp_dir=args.temp_dir,
        out_dir=args.out_dir
    )


if __name__ == "__main__":
    main()
