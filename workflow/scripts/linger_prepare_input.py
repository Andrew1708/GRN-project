import argparse
import os
import pandas as pd
import scipy.io
import scipy.sparse
import anndata
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare LINGER input files (RNA.txt, ATAC.txt, label.txt)")
    parser.add_argument('--rna_h5ad', type=str, required=True, help='Path to RNA AnnData (.h5ad) file')
    parser.add_argument('--label_column', type=str, required=True, help='Column name in adata.obs to use as label')
    parser.add_argument('--atac_mtx', type=str, required=True, help='Path to ATAC fragment matrix (.mtx)')
    parser.add_argument('--atac_barcodes', type=str, required=True, help='Path to ATAC barcodes.tsv')
    parser.add_argument('--atac_peaks', type=str, required=True, help='Path to ATAC peaks.tsv or BED file')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory')
    return parser.parse_args()

def fast_write_atac_txt(mtx, peaks, barcodes, output_path, chunk_size=100000):
    n_rows = mtx.shape[0]
    with open(output_path, 'w') as f:
        f.write('\t'.join(barcodes) + '\n')
        for start in range(0, n_rows, chunk_size):
            end = min(start + chunk_size, n_rows)
            chunk = mtx[start:end].toarray()
            lines = [
                f"{peaks[i]}\t" + '\t'.join(map(str, chunk[i - start]))
                for i in range(start, end)
            ]
            f.write('\n'.join(lines) + '\n')

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # --- RNA ---
    print("[1/3] Reading RNA AnnData...")
    adata = anndata.read_h5ad(args.rna_h5ad)
    if adata.raw:
        adata = adata.raw.to_adata()
    rna_df = adata.to_df().T  # genes x barcodes
    rna_path = os.path.join(args.output_dir, "RNA.txt")
    rna_df.to_csv(rna_path, sep="\t")
    print(f"✅ Saved RNA.txt: {rna_df.shape}")

    # --- label.txt ---
    print("[2/3] Generating label.txt...")
    if args.label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{args.label_column}' not in adata.obs.")
    labels = adata.obs[[args.label_column]].copy()
    labels['barcode'] = labels.index
    labels.insert(0, 'index', range(len(labels)))
    labels.columns = ['index', 'label', 'barcode']
    label_path = os.path.join(args.output_dir, "label.txt")
    labels.to_csv(label_path, sep="\t", index=False)
    print(f"✅ Saved label.txt: {labels.shape}")

    # --- ATAC ---
    print("[3/3] Reading ATAC matrix...")
    mtx = scipy.io.mmread(args.atac_mtx).tocsc()

    print("Reading barcodes and peaks...")
    barcodes = pd.read_csv(args.atac_barcodes, header=None)[0].tolist()
    peaks_df = pd.read_csv(args.atac_peaks, sep='\t', header=None)
    if peaks_df.shape[1] < 3:
        raise ValueError("Peaks file must have ≥3 columns (chrom, start, end)")
    peaks = (peaks_df[0] + ':' + peaks_df[1].astype(str) + '-' + peaks_df[2].astype(str)).tolist()

    print(f"Writing ATAC.txt to disk (this may take a few minutes)...")
    atac_path = os.path.join(args.output_dir, "ATAC.txt")
    fast_write_atac_txt(mtx, peaks, barcodes, atac_path)
    print(f"✅ Saved ATAC.txt: ({len(peaks)} peaks x {len(barcodes)} cells)")

if __name__ == "__main__":
    main()
