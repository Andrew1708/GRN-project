import os
import pickle
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
import pyranges as pr
import gc
import sys
from pathlib import Path


from pycisTopic.diff_features import impute_accessibility
from pycisTopic.gene_activity import get_gene_activity

sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.scrna_preprocessing import *

def process_sample(
    cistopic_path,
    rna_path,
    chromsizes_path,
    tss_path,
    out_dir,
    n_top_genes=None
):
    os.makedirs(out_dir, exist_ok=True)

    # Load cisTopic object
    with open(cistopic_path, 'rb') as f:
        cistopic_object = pickle.load(f)
    print(f"Loaded cisTopic object from {cistopic_path}")

    # Load chromosome sizes
    chromsizes = pd.read_table(chromsizes_path)
    chromsizes.rename(columns={"# ucsc": "Chromosome", "length": "End"}, inplace=True)
    chromsizes["Start"] = 0
    chromsizes = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])

    # Load TSS annotations
    pr_annotation = pd.read_table(tss_path).rename(columns={"Name": "Gene", "# Chromosome": "Chromosome"})
    pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
    pr_annotation = pr.PyRanges(pr_annotation)

    # Load and preprocess RNA
    adata_rna = sc.read_h5ad(rna_path)
    print(f"✅ Loaded RNA data from {rna_path}")
    #adata_rna = select_hvgs(adata_rna, n_top_genes=n_top_genes)

    # Impute accessibility
    imputed_acc_obj = impute_accessibility(
        cistopic_object,
        selected_cells=None,
        selected_regions=None,
        scale_factor=10**6
    )

    # Calculate gene activity
    gene_act, _ = get_gene_activity(
        imputed_acc_obj,
        pr_annotation,
        chromsizes,
        use_gene_boundaries=True,
        upstream=[1000, 100000],
        downstream=[1000, 100000],
        distance_weight=True,
        decay_rate=1,
        extend_gene_body_upstream=10000,
        extend_gene_body_downstream=500,
        gene_size_weight=False,
        gene_size_scale_factor='median',
        remove_promoters=False,
        average_scores=True,
        scale_factor=1,
        extend_tss=[10, 10],
        gini_weight=True,
        return_weights=True,
        project='Gene_activity'
    )

    # Convert gene activity to AnnData
    X = gene_act.mtx.T
    obs = pd.DataFrame(index=gene_act.cell_names)
    var = pd.DataFrame(index=gene_act.feature_names)
    adata_atac = AnnData(X=X, obs=obs, var=var)
    adata_atac.obs = cistopic_object.cell_data.loc[adata_atac.obs.index].astype("str")

    # === Filter genes (intersection) ===
    common_genes = adata_rna.var_names.intersection(adata_atac.var_names)

    # Filter genes (intersection)
    adata_rna_filtered = adata_rna[:, common_genes].copy()
    adata_atac_filtered = adata_atac[:, common_genes].copy()

    # Remove RNA cells with all-zero counts
    nonzero_rna_cells = np.array((adata_rna_filtered.X.sum(axis=1) > 0)).flatten()
    adata_rna_filtered = adata_rna_filtered[nonzero_rna_cells, :].copy()

    # Optional: Rename metadata columns if needed
    if "Classified_Celltype" in adata_rna_filtered.obs.columns:
        adata_rna_filtered.obs.rename(columns={"Classified_Celltype": "CellType"}, inplace=True)

    # Save filtered data
    atac_filtered_out_path = os.path.join(out_dir, "gene_act_matrix.h5ad")
    rna_filtered_out_path = os.path.join(out_dir, "rna_filtered.h5ad")

    # Replace NaNs with 0s
    adata_atac_filtered.X = np.nan_to_num(adata_atac_filtered.X, nan=0.0)
    adata_rna_filtered.X = np.nan_to_num(adata_rna_filtered.X, nan=0.0)

    adata_atac_filtered.write(atac_filtered_out_path)
    adata_rna_filtered.write(rna_filtered_out_path)
    print(f"Saved filtered ATAC AnnData to: {atac_filtered_out_path}")
    print(f"Saved filtered RNA AnnData to: {rna_filtered_out_path}")


    return adata_rna, adata_atac, adata_rna_filtered, adata_atac_filtered

def parse_args():
    parser = argparse.ArgumentParser(description="Process a single scATAC-seq and scRNA-seq sample.")
    parser.add_argument("--cistopic_path", required=True, help="Path to cisTopic .pkl file")
    parser.add_argument("--rna_path", required=True, help="Path to scRNA-seq .h5ad file")
    parser.add_argument("--chromsizes_path", required=True, help="Path to chromosome sizes .tsv file")
    parser.add_argument("--tss_path", required=True, help="Path to TSS annotation .bed file")
    parser.add_argument("--celltype_col", default="Classified_Celltype", help="Column name for cell types in the metadata")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    return parser.parse_args()

def main():
    args = parse_args()
    adata_rna, adata_atac, adata_rna_filtered, adata_atac_filtered = process_sample(
        args.cistopic_path,
        args.rna_path,
        args.chromsizes_path,
        args.tss_path,
        args.out_dir
    )
    # At this point, all four AnnData objects are available in memory for further analysis if needed.

if __name__ == "__main__":
    main()
