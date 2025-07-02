import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co
import mudata as mu
import scipy.io


def parse_args():
    parser = argparse.ArgumentParser(description="Run CellOracle with optional ATAC-based prior and export MuData.")
    parser.add_argument("--input", type=str, required=True, help="Path to input .h5ad file (scRNA-seq)")
    parser.add_argument("--base_grn", type=str, help="Path to base GRN .parquet file (required if --use_atac is set)")
    parser.add_argument("--use_atac", action="store_true", help="Flag to use ATAC-derived base GRN")
    parser.add_argument("--output", type=str, required=True, help="Output MuData file path (.h5mu)")
    parser.add_argument("--atac_matrix", type=str, help="Path to ATAC peak matrix file (if using ATAC)")
    parser.add_argument("--atac_barcodes", type=str, help="Path to ATAC barcodes file (if using ATAC)")
    parser.add_argument("--atac_peaks", type=str, help="Path to ATAC peaks file in BED format (if using ATAC)")
    parser.add_argument("--downsample", type=int, default=30000, help="Number of cells to downsample")
    return parser.parse_args()


def flatten_links_object(links_obj):
    """Convert CellOracle Links object to a flat DataFrame with cluster labels."""
    records = []
    for cluster_name, df in links_obj.links_dict.items():
        df = df.copy()
        df["cluster"] = cluster_name
        records.append(df)
    return pd.concat(records, ignore_index=True)


def expand_links_with_CREs(links, base_GRN):
    # Prepare base_GRN in long format: one row per (TF, gene)
    tf_columns = [col for col in base_GRN.columns if col not in ["gene_short_name", "chromosome", "start", "end"]]

    melted = base_GRN.melt(
        id_vars=["gene_short_name", "chromosome", "start", "end"],
        value_vars=tf_columns,
        var_name="source",
        value_name="score"
    )

    # Filter non-zero scores (active TF-gene interactions)
    melted = melted[melted["score"] != 0]

    # Merge links with base_GRN annotations
    annotated = pd.merge(
        links,
        melted,
        how="left",
        left_on=["source", "target"],
        right_on=["source", "gene_short_name"]
    )

    # Drop extra column
    annotated = annotated.drop(columns=["score", "gene_short_name"])

    # For rows with no CRE info, fill NaNs as before
    for col in ["chromosome", "start", "end"]:
        if col not in annotated.columns:
            annotated[col] = np.nan

    return annotated


def main():
    args = parse_args()

    # Load scRNA-seq data
    adata = sc.read_h5ad(args.input)

    # Downsample cells if needed
    if adata.shape[0] > args.downsample:
        sc.pp.subsample(adata, n_obs=args.downsample, random_state=123)

    # Load base GRN
    if args.use_atac:
        if args.base_grn is None:
            raise ValueError("Please provide --base_grn when using --use_atac")
        base_GRN = pd.read_parquet(args.base_grn)
    else:
        base_GRN = co.data.load_human_promoter_base_GRN(version="hg38_gimmemotifsv5_fpr2")

    if "peak_id" in base_GRN.columns:
        peak_coords = base_GRN["peak_id"].str.split("_", expand=True)

        is_numeric_start = peak_coords[1].str.isnumeric()
        is_numeric_end = peak_coords[2].str.isnumeric()
        valid_rows = is_numeric_start & is_numeric_end

        if not valid_rows.any():
            raise ValueError("No valid 'peak_id' entries with numeric coordinates found.")

        peak_coords = peak_coords[valid_rows]
        base_GRN = base_GRN.loc[valid_rows].copy()
        base_GRN["chromosome"] = peak_coords[0].values
        base_GRN["start"] = peak_coords[1].astype(int).values
        base_GRN["end"] = peak_coords[2].astype(int).values
    else:
        raise ValueError("Expected 'peak_id' column in base_GRN, but it was not found.")

    print("Creating MuData object...")
    modalities = {"rna": adata}
    if args.use_atac:
        # Load sparse peak matrix (.mtx), shape: (peaks x cells) → transpose to (cells x peaks)
        peak_matrix = scipy.io.mmread(str(args.atac_matrix)).T.tocsr()

        # Load barcodes (1 per line, matching matrix rows after transpose)
        atac_barcodes = pd.read_csv(args.atac_barcodes, header=None)[0].tolist()

        # Load peaks from BED file: no header, 3 columns (chrom, start, end)
        peaks_df = pd.read_csv(
            args.atac_peaks,
            sep="\t",
            header=None,
            names=["seqnames", "start", "end"],
            usecols=[0, 1, 2]
        )
        peaks = (peaks_df["seqnames"] + ":" + peaks_df["start"].astype(str) + "-" + peaks_df["end"].astype(str)).tolist()

        # Sanity checks
        assert len(atac_barcodes) == peak_matrix.shape[0], \
            f"Barcodes count ({len(atac_barcodes)}) doesn't match matrix rows ({peak_matrix.shape[0]})"
        assert len(peaks) == peak_matrix.shape[1], \
            f"Peaks count ({len(peaks)}) doesn't match matrix columns ({peak_matrix.shape[1]})"

        # Build ATAC AnnData
        atac_mod = sc.AnnData(
            X=peak_matrix,
            obs=pd.DataFrame(index=atac_barcodes),
            var=pd.DataFrame(index=peaks)
        )

        modalities["atac"] = atac_mod

    # Instantiate Oracle
    oracle = co.Oracle()
    adata.X = adata.layers["raw_count"].copy()
    print("Metadata columns :", list(adata.obs.columns))
    print("Dimensional reduction: ", list(adata.obsm.keys()))

    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name="louvain",
        embedding_name="X_draw_graph_fa"
    )

    oracle.import_TF_data(TF_info_matrix=base_GRN)

    oracle.perform_PCA()

    # Estimate number of PCs
    plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
    diffs = np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002)
    n_comps = np.where(diffs)[0][0] if np.any(diffs) else 30
    n_comps = min(n_comps, 50)

    n_cell = oracle.adata.shape[0]
    print(f"Cell number: {n_cell}")
    k = int(0.025 * n_cell)
    print(f"Auto-selected k: {k}")

    oracle.knn_imputation(
        n_pca_dims=n_comps,
        k=k,
        balanced=True,
        b_sight=k * 8,
        b_maxl=k * 4,
        n_jobs=4
    )

    print("Inferring GRNs...")
    raw_links = oracle.get_links(
        cluster_name_for_GRN_unit="louvain",
        alpha=10,
        verbose_level=10
    )

    print("Flattening GRN links...")
    links = flatten_links_object(raw_links)

    print("Annotating links with CREs...")
    links = expand_links_with_CREs(links, base_GRN)

    mdata = mu.MuData(modalities)
    mdata.uns["celloracle_links"] = links

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Saving MuData to {args.output}")
    mdata.write(args.output)


if __name__ == "__main__":
    main()
