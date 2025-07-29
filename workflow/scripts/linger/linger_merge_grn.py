import os
import gc
import time
import argparse
import pandas as pd
import numpy as np
import mudata as md
import anndata as ad
import pyarrow.dataset as ds


def parse_args():
    parser = argparse.ArgumentParser(description="TF–TG–RE GRN Inference Pipeline")
    parser.add_argument("--tf_re", required=True, help="Path to TF–RE binding file (TSV)")
    parser.add_argument("--re_tg", required=True, help="Path to RE–TG linkage file (TSV)")
    parser.add_argument("--tf_tg", required=True, help="Path to TF–TG association file (TSV)")
    parser.add_argument("--rna", required=True, help="Path to RNA modality (.h5ad)")
    parser.add_argument("--atac", required=True, help="Path to ATAC modality (.h5ad)")
    parser.add_argument("--outdir", required=True, help="Output directory to store results")
    parser.add_argument("--batch_size", type=int, default=1500, help="Number of genes per TF–TG batch")
    parser.add_argument("--meta", type=str, required=True, help="Path to metadata file containing cell types")
    parser.add_argument("--celltype_col", type=str, default="Classified_Celltype", help="Column name for cell types in RNA/ATAC obs")
    return parser.parse_args()



def transfer_celltypes_from_anndata(
    rna: ad.AnnData,
    meta_adata: ad.AnnData,
    celltype_col: str = "Classified_Celltype"):
    print("Transferring cell type annotations from metadata AnnData...")

    # Normalize metadata barcodes
    normalized_meta_barcodes = meta_adata.obs_names.str.split("-").str[0]
    meta_df = meta_adata.obs.copy()
    meta_df["normalized_barcode"] = normalized_meta_barcodes
    meta_df = meta_df.set_index("normalized_barcode")

    # Extract RNA barcodes from rna.obs["barcode"]
    if "barcode" not in rna.obs.columns:
        raise ValueError("RNA AnnData must contain a 'barcode' column in .obs")

    rna_barcodes = rna.obs["barcode"]
    unmatched = rna_barcodes[~rna_barcodes.isin(meta_df.index)]

    if len(unmatched) > 0:
        raise ValueError(
            f"{len(unmatched)} RNA cell barcodes could not be matched in the metadata AnnData. Halting.\n"
            f"Examples of unmatched barcodes: {unmatched.unique()[:5].tolist()}"
        )

    # Assign in place
    rna.obs[celltype_col] = rna_barcodes.map(meta_df[celltype_col])
    print(f"Cell type column '{celltype_col}' successfully transferred to RNA AnnData.")


def load_tf_re_re_tg(tf_re_path, re_tg_path):
    print("📦 Loading RE–TG and TF–RE matrices...")
    df_re_tg = pd.read_csv(re_tg_path, sep="\t", header=None, names=["Region", "Gene", "tg_re_score"])
    df_tf_re = pd.read_csv(tf_re_path, sep="\t", index_col=0)
    return df_tf_re, df_re_tg


def count_genes(tf_tg_path):
    with open(tf_tg_path) as f:
        return sum(1 for _ in f) - 1


def process_batches(tf_tg_path, df_re_tg, df_tf_re, batch_size, parquet_dir, alpha=0.75):
    os.makedirs(parquet_dir, exist_ok=True)
    n_genes_total = count_genes(tf_tg_path)
    n_batches = (n_genes_total + batch_size - 1) // batch_size
    print(f"📊 Total genes: {n_genes_total} → {n_batches} batches of {batch_size} genes")

    tf_tg_iter = pd.read_csv(tf_tg_path, sep="\t", index_col=0, chunksize=batch_size)

    for i, df_tf_tg_chunk in enumerate(tf_tg_iter):
        print(f"\n🔄 [Batch {i+1}/{n_batches}] Processing {len(df_tf_tg_chunk)} genes...")
        t0 = time.time()

        filtered_tf_tg_rows = []

        for gene, row in df_tf_tg_chunk.iterrows():
            abs_row = row.abs()

            if abs_row.isna().all() or abs_row.max() == 0:
                continue

            avg = abs_row.mean()
            max_val = abs_row.max()
            cutoff = avg + alpha * (max_val - avg)

            selected = abs_row[abs_row > cutoff]

            for tf in selected.index:
                score = row[tf]  # retain original signed score
                filtered_tf_tg_rows.append((gene, tf, score))

        df_tf_tg_long = pd.DataFrame(filtered_tf_tg_rows, columns=["Gene", "TF", "tf_tg_score"])
        print(f"    → {len(df_tf_tg_long)} filtered TF–TG interactions")

        # Merge with RE–TG
        merged = pd.merge(df_re_tg, df_tf_tg_long, on="Gene", how="inner")
        print(f"    → {len(merged)} merged rows")

        # Subset TF–RE
        needed_regions = merged["Region"].unique()
        needed_tfs = merged["TF"].unique()
        df_tf_re_sub = df_tf_re.loc[
            df_tf_re.index.intersection(needed_regions),
            df_tf_re.columns.intersection(needed_tfs)
        ]

        # Melt TF–RE
        df_tf_re_flat = df_tf_re_sub.reset_index().melt(
            id_vars="index", var_name="TF", value_name="tf_re_score"
        ).rename(columns={"index": "Region"})
        print(f"    → {len(df_tf_re_flat)} TF–RE entries")

        # Final merge
        grn = pd.merge(merged, df_tf_re_flat, on=["Region", "TF"], how="left")
        grn["tf_re_score"] = grn["tf_re_score"].fillna(0)
        grn = grn[["TF", "Gene", "Region", "tf_tg_score", "tg_re_score", "tf_re_score"]]
        print(f"    → GRN rows: {len(grn)}")

        # Write to parquet
        grn.to_parquet(os.path.join(parquet_dir, f"grn_batch_{i+1:03d}.parquet"), index=False)

        # Cleanup
        del df_tf_tg_chunk, df_tf_tg_long, merged, df_tf_re_sub, df_tf_re_flat, grn
        gc.collect()
        print(f"  ✅ Batch {i+1} done in {time.time() - t0:.2f}s")


def combine_parquet_batches(parquet_dir):
    print("\n🧬 Loading all GRN batches (via Arrow dataset)...")
    dataset = ds.dataset(parquet_dir, format="parquet")
    table = dataset.to_table()
    grn_all = table.to_pandas()
    print(f"📐 Combined GRN size: {len(grn_all)} rows")
    return grn_all


def build_mudata(grn_df, rna_path, atac_path, meta_path, celltype_col, output_path):
    print("💾 Building MuData container...")
    ad_rna = ad.read_h5ad(rna_path)
    ad_atac = ad.read_h5ad(atac_path)
    metadata = ad.read_h5ad(meta_path)

    ad_atac.obs = ad_atac.obs.reindex(ad_rna.obs.index)

    transfer_celltypes_from_anndata(ad_rna, metadata, celltype_col)

    mdata = md.MuData({"rna": ad_rna, "atac": ad_atac})
    mdata.uns["grn"] = grn_df

    print(f"💾 Saving MuData to {output_path}...")
    mdata.write(output_path)
    print("✅ All done.")


def main():
    args = parse_args()

    parquet_dir = os.path.join(args.outdir, "grn_parquet_batches")
    mudata_output = os.path.join(args.outdir, "mdata.h5mu")
    os.makedirs(args.outdir, exist_ok=True)

    df_tf_re, df_re_tg = load_tf_re_re_tg(args.tf_re, args.re_tg)
    process_batches(args.tf_tg, df_re_tg, df_tf_re, args.batch_size, parquet_dir)
    grn_all = combine_parquet_batches(parquet_dir)
    build_mudata(grn_all, args.rna, args.atac, args.meta, args.celltype_col, mudata_output)


if __name__ == "__main__":
    main()
