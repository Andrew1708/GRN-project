import os
import sys
import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.io import mmread
from mudata import MuData

# Add path to custom utilities
sys.path.append(str(Path(__file__).resolve().parent.parent))


def prepare_scconfluence_mdata(rna_adata, atac_fragments, atac_peaks, atac_barcodes, project_name):
    """
    Prepares a MuData object for scConfluence.
    - Leaves RNA structure untouched except for cm_features
    - Uses HVGs for RNA PCA
    - Uses all peaks for ATAC PCA
    """

    # --- 1. Construct ATAC AnnData
    atac_adata = ad.AnnData(X=atac_fragments)
    atac_adata.var_names = atac_peaks
    atac_adata.obs_names = [f"{bc}___{project_name}" for bc in atac_barcodes]
    atac_adata.var_names.name = None

    # --- 2. RNA: keep only HVGs for PCA
    sc.pp.normalize_total(rna_adata)
    sc.pp.log1p(rna_adata)
    sc.pp.highly_variable_genes(rna_adata, flavor="seurat_v3", n_top_genes=2000, subset=True)
    sc.pp.pca(rna_adata, n_comps=50)
    rna_adata.obsm["cm_features"] = rna_adata.obsm["X_pca"]

    # --- 3. ATAC: PCA on log1p normalized counts
    sc.pp.normalize_total(atac_adata)
    sc.pp.log1p(atac_adata)
    sc.pp.pca(atac_adata, n_comps=50)
    atac_adata.obsm["cm_features"] = atac_adata.obsm["X_pca"]

    # --- 4. Combine into MuData
    mdata = MuData({"rna": rna_adata, "atac": atac_adata})
    return mdata


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare MuData object for scConfluence.")
    parser.add_argument("--rna_path", required=True, help="Path to RNA .h5ad file (AnnData with .X as counts)")
    parser.add_argument("--atac_fragments", required=True, help="Path to ATAC matrix (.mtx file)")
    parser.add_argument("--atac_peaks", required=True, help="Path to ATAC peak BED file")
    parser.add_argument("--atac_barcodes", required=True, help="Path to ATAC cell barcode .tsv file")
    parser.add_argument("--project_name", required=True, help="Project name for saving the MuData object")
    parser.add_argument("--out_path", required=True, help="Output path to store .h5mu file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(os.path.dirname(args.out_path), exist_ok=True)

    # --- Load RNA AnnData
    rna_adata = sc.read_h5ad(args.rna_path)

    # --- Load ATAC matrix (.mtx format)
    atac_fragments = mmread(args.atac_fragments).transpose().tocsr()

    # --- Load peaks from BED file (skip header line if present)
    with open(args.atac_peaks) as f:
        first_line = next(f)
        if first_line.lower().startswith("chrom"):
            atac_peaks = [f"{line.split()[0]}:{line.split()[1]}-{line.split()[2]}" for line in f]
        else:
            atac_peaks = [f"{first_line.split()[0]}:{first_line.split()[1]}-{first_line.split()[2]}"] + \
                         [f"{line.split()[0]}:{line.split()[1]}-{line.split()[2]}" for line in f]

    # --- Load barcodes from TSV
    with open(args.atac_barcodes) as f:
        atac_barcodes = [line.strip() for line in f]

    # --- Sanity check dimensions
    n_cells_frag, n_peaks_frag = atac_fragments.shape
    assert n_cells_frag == len(atac_barcodes), (
        f"Fragment matrix has {n_cells_frag} rows but {len(atac_barcodes)} barcodes."
    )
    assert n_peaks_frag == len(atac_peaks), (
        f"Fragment matrix has {n_peaks_frag} columns but {len(atac_peaks)} peaks."
    )

    # --- Prepare MuData object for scConfluence
    mdata = prepare_scconfluence_mdata(rna_adata, atac_fragments, atac_peaks, atac_barcodes, args.project_name)

    # --- Save MuData
    mdata.write(args.out_path)
    print(f"Saved MuData to: {args.out_path}")
