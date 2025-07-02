import os
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from scipy.io import mmread

from scvi.model import MULTIVI
from mudata import MuData

# Add path to custom utilities
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.mutual_nearest_neighbors import mutual_nearest_neighbors


def run_multivi_integration(rna_adata, atac_fragments, atac_peaks, atac_barcodes):
    """
    Performs integration of scRNA and scATAC using MultiVI and returns embeddings + cell names.
    """

    # --- 1. Prepare ATAC AnnData
    atac_adata = ad.AnnData(X=atac_fragments)
    atac_adata.var_names = atac_peaks
    atac_adata.obs_names = atac_barcodes

    # --- 2. Fix RNA .var structure
    rna_adata.var["feature_type"] = "gene"

    if rna_adata.var_names.name == "gene_ids":
        rna_adata.var["gene_ids"] = rna_adata.var.index
        rna_adata.var_names = rna_adata.var["gene_ids"]

    if "gene_ids" not in rna_adata.var.columns:
        rna_adata.var["gene_ids"] = rna_adata.var_names

    rna_adata.var_names.name = None

    # --- 3. Fix ATAC .var structure
    atac_adata.var["feature_type"] = "peak"
    atac_adata.var_names.name = None

    # --- 4. Annotate modality labels
    rna_adata.obs["modality"] = "RNA"
    atac_adata.obs["modality"] = "ATAC"

    # --- 5. Create MuData object
    mdata = MuData({"rna": rna_adata, "atac": atac_adata})
    print(f"MuData created with modalities: {list(mdata.mod.keys())}")

    # --- 6. Setup MultiVI
    MULTIVI.setup_mudata(
        mdata,
        modalities={"rna_layer": "rna", "atac_layer": "atac"},
        batch_key="modality"  # optional: passed to each modality's obs
    )


    # --- 7. Train model
    model = MULTIVI(mdata)
    model.train()

    # --- 8. Get latent embeddings
    latents = model.get_latents()
    rna_names = mdata.mod["rna"].obs_names.tolist()
    atac_names = mdata.mod["atac"].obs_names.tolist()
    rna_latents = latents.loc[rna_names].values
    atac_latents = latents.loc[atac_names].values

    return rna_latents, atac_latents, rna_names, atac_names


def parse_args():
    parser = argparse.ArgumentParser(description="Run MultiVI diagonal integration and compute MNN matches.")
    parser.add_argument("--rna_path", required=True, help="Path to RNA .h5ad file (AnnData with .X as counts)")
    parser.add_argument("--atac_fragments", required=True, help="Path to ATAC matrix (.mtx file)")
    parser.add_argument("--atac_peaks", required=True, help="Path to ATAC peak BED file")
    parser.add_argument("--atac_barcodes", required=True, help="Path to ATAC cell barcode .tsv file")
    parser.add_argument("--out_dir", required=True, help="Output directory to store matches.csv")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

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

    # --- Check dimensions with detailed errors
    n_cells_frag, n_peaks_frag = atac_fragments.shape
    n_barcodes = len(atac_barcodes)
    n_peaks = len(atac_peaks)

    assert n_cells_frag == n_barcodes, (
        f"Fragment matrix has {n_cells_frag} rows but {n_barcodes} barcodes."
    )
    assert n_peaks_frag == n_peaks, (
        f"Fragment matrix has {n_peaks_frag} columns but {n_peaks} peaks."
    )

    # --- Run MultiVI integration
    rna_latents, atac_latents, rna_names, atac_names = run_multivi_integration(
        rna_adata, atac_fragments, atac_peaks, atac_barcodes
    )

    # --- Perform MNN matching
    mnn_results = mutual_nearest_neighbors(rna_latents, atac_latents, rna_names, atac_names)

    # --- Save matches to disk
    out_csv = os.path.join(args.out_dir, "matches.csv")
    mnn_results.to_csv(out_csv, index=False)
    print(f"Saved MNN results to: {out_csv}")
