#!/usr/bin/env python3

import argparse
from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd

# Optional: add project utils path if needed
# from utils.scrna_preprocessing import * 

def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocess scRNA-seq data and compute Diffusion Map and PAGA graph using Scanpy."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to input .h5ad AnnData file."
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output .h5ad AnnData file."
    )
    parser.add_argument(
        "--n_top_genes",
        type=int,
        default=2000,
        help="Number of highly variable genes to retain (default: 2000)."
    )
    return parser.parse_args()


def preprocess_adata(adata, n_top_genes):
    print("Filtering genes with min count ≥ 1...")
    sc.pp.filter_genes(adata, min_counts=1)

    print("Normalizing per cell (library size normalization)...")
    sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

    print(f"Selecting top {n_top_genes} highly variable genes...")
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X,
        flavor='cell_ranger',
        n_top_genes=n_top_genes,
        log=False
    )
    adata = adata[:, filter_result.gene_subset]

    print("Renormalizing per cell after gene filtering...")
    sc.pp.normalize_per_cell(adata)

    print("Storing raw counts...")
    adata.raw = adata
    adata.layers["raw_count"] = adata.raw.X.copy()

    print("Applying log1p transformation and scaling...")
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    return adata


def compute_diffmap_paga(adata):
    print("Running PCA...")
    sc.tl.pca(adata, svd_solver='arpack')

    print("Computing initial neighbors for diffusion map...")
    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

    print("Computing diffusion map...")
    sc.tl.diffmap(adata)

    print("Recomputing neighbors based on diffusion map...")
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

    print("Running Louvain clustering...")
    sc.tl.louvain(adata, resolution=0.8)

    print("Constructing PAGA graph...")
    sc.tl.paga(adata, groups='louvain')

    print("Computing PAGA layout positions...")
    sc.pl.paga(adata, show=False)  # This populates adata.uns['paga']['pos']

    print("Drawing layout with PAGA initialization...")
    sc.tl.draw_graph(adata, init_pos='paga', random_state=123)

    print("Running ForceAtlas2 layout...")
    sc.tl.draw_graph(adata, layout="fa")  # ForceAtlas2 layout
    
    return adata


def main():
    args = parse_args()

    print(f"Loading input AnnData: {args.input}")
    adata = sc.read_h5ad(args.input)

    adata = preprocess_adata(adata, n_top_genes=args.n_top_genes)
    adata = compute_diffmap_paga(adata)

    print(f"Saving output AnnData to: {args.output}")
    adata.write_h5ad(args.output)

    print("Processing completed.")


if __name__ == "__main__":
    main()
