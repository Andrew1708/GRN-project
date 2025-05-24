import os
import pickle
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse, coo_matrix
import pyranges as pr
import gc

from pycisTopic.diff_features import impute_accessibility
from pycisTopic.gene_activity import get_gene_activity

def normalize(adata):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata

def select_hvgs(adata, n_top_genes=2000):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var.highly_variable].copy()
    return adata

def preprocess(adata, n_top_genes=2000):
    adata = normalize(adata)
    adata = select_hvgs(adata, n_top_genes=n_top_genes)
    sc.pp.scale(adata)
    return adata

def build_and_save_matrices(cistopic_path, h5ad_path, chromsizes_path, tss_path, out_dir, threshold=0):
    # Load inputs
    with open(cistopic_path, "rb") as f:
        cistopic_object = pickle.load(f)
    adata_rna = sc.read_h5ad(h5ad_path)
    chromsizes = pd.read_table(chromsizes_path)
    chromsizes.rename(columns={"# ucsc": "Chromosome", "length": "End"}, inplace=True)
    chromsizes["Start"] = 0
    chromsizes = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])
    pr_annotation = pd.read_table(tss_path).rename(columns={"Name": "Gene", "# Chromosome": "Chromosome"})
    pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
    pr_annotation = pr.PyRanges(pr_annotation)

    # Preprocess RNA
    adata_rna = preprocess(adata_rna)

    # Impute accessibility
    imputed_acc_obj = impute_accessibility(cistopic_object, scale_factor=10**6)

    # Gene activity
    _, weights = get_gene_activity(
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

    # Normalize weights
    w_min, w_max = weights["Weight"].min(), weights["Weight"].max()
    weights["norm_weight"] = (weights["Weight"] - w_min) / (w_max - w_min) if w_max > w_min else 1
    weights_filt = weights[weights["norm_weight"] >= threshold]

    # Define features
    filtered_regions = np.array(weights_filt["Name"].unique())
    cistopic_regions = np.array(cistopic_object.region_names)
    unique_regions = np.sort(np.intersect1d(cistopic_regions, filtered_regions))
    filtered_genes = np.array(weights_filt["Gene"].unique())
    rna_genes = np.array(adata_rna.var_names)
    unique_genes = np.sort(np.intersect1d(rna_genes, filtered_genes))

    region_to_index = {region: i for i, region in enumerate(unique_regions)}
    gene_to_index = {gene: j for j, gene in enumerate(unique_genes)}
    positive_entries = weights_filt[
        (weights_filt["Gene"].isin(unique_genes)) &
        (weights_filt["Name"].isin(unique_regions))
    ]

    # Region-to-gene matrix
    row_indices = positive_entries["Name"].map(region_to_index).values
    col_indices = positive_entries["Gene"].map(gene_to_index).values
    data = np.ones(len(positive_entries), dtype=np.int8)
    region2gene = coo_matrix((data, (row_indices, col_indices)), shape=(len(unique_regions), len(unique_genes)))
    region2gene_binary = (region2gene > 0).astype(np.int8).tocoo()
    df_rg_dense = pd.DataFrame(0, index=unique_regions, columns=unique_genes)
    for r, c, v in zip(region2gene_binary.row, region2gene_binary.col, region2gene_binary.data):
        df_rg_dense.iloc[r, c] = v
    df_rg_dense.to_csv(os.path.join(out_dir, f"{cistopic_object.project}_region2gene_dense.csv"))

    # RNA dense
    adata_rna_filtered = adata_rna[:, unique_genes].copy()
    if issparse(adata_rna_filtered.X):
        rna_dense = pd.DataFrame(adata_rna_filtered.X.toarray(), index=adata_rna_filtered.obs_names, columns=unique_genes)
    else:
        rna_dense = pd.DataFrame(adata_rna_filtered.X, index=adata_rna_filtered.obs_names, columns=unique_genes)
    rna_dense.to_csv(os.path.join(out_dir, f"{cistopic_object.project}_RNA_counts_dense.csv"))

    # ATAC dense
    region_indices = [np.where(cistopic_regions == region)[0][0] for region in unique_regions]
    frag_filtered = cistopic_object.fragment_matrix[region_indices, :].T
    if issparse(frag_filtered):
        atac_dense = pd.DataFrame(frag_filtered.toarray(), index=cistopic_object.cell_names, columns=unique_regions)
    else:
        atac_dense = pd.DataFrame(frag_filtered, index=cistopic_object.cell_names, columns=unique_regions)
    atac_dense.to_csv(os.path.join(out_dir, f"{cistopic_object.project}_ATAC_counts_dense.csv"))

    print("Matrices saved successfully.")
    gc.collect()

def parse_args():
    parser = argparse.ArgumentParser(description="Generate dense matrices for diagonal integration.")
    parser.add_argument("--cistopic_path", required=True, help="Path to cistopic object .pkl file")
    parser.add_argument("--h5ad_path", required=True, help="Path to RNA .h5ad file")
    parser.add_argument("--chromsizes_path", required=True, help="Path to chromosome sizes .tsv file")
    parser.add_argument("--tss_path", required=True, help="Path to TSS annotation .bed file")
    parser.add_argument("--out_dir", required=True, help="Output directory for dense matrices")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    build_and_save_matrices(
        args.cistopic_path,
        args.h5ad_path,
        args.chromsizes_path,
        args.tss_path,
        args.out_dir
    )

if __name__ == "__main__":
    main()
