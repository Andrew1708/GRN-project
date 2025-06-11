import os
import subprocess
from pathlib import Path
from collections import defaultdict
import argparse
import sys
import scanpy as sc

# Add the parent directory (assuming utils is in ../utils)
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.mutual_nearest_neighbors import *

def parse_args():
    parser = argparse.ArgumentParser(description="Generate dense matrices for diagonal integration.")
    parser.add_argument("--rna_path", required=True, help="Path to RNA .h5ad file")
    parser.add_argument("--atac_path", required=True, help="Path to ATAC .pkl file")
    parser.add_argument("--out_dir", required=True, help="Output directory for dense matrices")
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Load RNA data
    rna_adata = sc.read_h5ad(args.rna_path)
    rna_embeddings = rna_adata.obsm['Embedding']
    rna_names = rna_adata.obs_names

    # Load ATAC data
    atac_adata = sc.read_h5ad(args.atac_path)
    atac_embeddings = atac_adata.obsm['Embedding']
    atac_names = atac_adata.obs_names

    # Perform mutual nearest neighbors matching
    mnn_results = mutual_nearest_neighbors(rna_embeddings, atac_embeddings, rna_names, atac_names)

    # Save results to output directory
    mnn_results.to_csv(os.path.join(args.out_dir, 'matches.csv'), index=False)