import argparse
import os
import pandas as pd
from scipy.sparse import csc_matrix
import scanpy as sc

from LingerGRN.preprocess import get_adata, preprocess
from LingerGRN.pseudo_bulk import pseudo_bulk
import LingerGRN.LINGER_tr as LINGER_tr
import LingerGRN.LL_net as LL_net

def clean_barcode_columns(columns):
    return columns.str.split('___').str[0].str.replace('-1$', '', regex=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Run LINGER pipeline and save results in MuData")
    parser.add_argument('--data_dir', type=str, required=True, help='Directory with RNA.txt, ATAC.txt, label.txt')
    parser.add_argument('--datadir', type=str, required=True, help='Directory with bulk GRN data (not used by LINGER)')
    parser.add_argument('--outdir', type=str, required=True, help='Unique output directory for this run')
    parser.add_argument('--method', type=str, default='LINGER', choices=['LINGER', 'baseline'])
    parser.add_argument('--genome', type=str, default='hg38')
    parser.add_argument('--network', type=str, default='cell population')
    parser.add_argument('--celltype', type=str, default='all')
    parser.add_argument('--activation', type=str, default='ReLU', choices=['ReLU', 'sigmoid', 'tanh'])
    return parser.parse_args()

def main():
    args = parse_args()

    # Resolve absolute paths so that chdir doesn't break them
    args.data_dir = f"{os.path.abspath(args.data_dir)}/"
    args.datadir = f"{os.path.abspath(args.datadir)}/"
    args.outdir = f"{os.path.abspath(args.outdir)}/"
    

    # Load input matrices and metadata
    label = pd.read_csv(os.path.join(args.data_dir, 'label.txt'), sep='\t')
    RNA = pd.read_csv(os.path.join(args.data_dir, 'RNA.txt'), sep='\t', index_col=0)
    ATAC = pd.read_csv(os.path.join(args.data_dir, 'ATAC.txt'), sep='\t', index_col=0)

    # Setup per-run working directory and data directory
    run_workdir = os.path.join(args.outdir, "linger_wrk")
    data_dir = os.path.join(run_workdir, "data")
    os.makedirs(data_dir, exist_ok=True)


    #TODO meter isto correcto no preparate input => dps de verficar se isto é o correto
    RNA.columns = clean_barcode_columns(RNA.columns)
    ATAC.columns = clean_barcode_columns(ATAC.columns)
    # Clean the "barcode" column and rename it to "barcode_use"
    label['barcode_use'] = label['barcode'].str.split('___').str[0].str.replace('-1$', '', regex=True)
    # Drop the original "index" column if it exists
    label = label.drop(columns=['index'], errors='ignore')
    # Reorder columns: "label", "barcode_use"
    label = label[['label', 'barcode_use']]

    # Combine and label features
    matrix = csc_matrix(pd.concat([RNA, ATAC], axis=0).values)
    features = pd.DataFrame(RNA.index.tolist() + ATAC.index.tolist(), columns=[1])
    features[2] = ['Gene Expression'] * RNA.shape[0] + ['Peaks'] * ATAC.shape[0]
    barcodes = pd.DataFrame(RNA.columns.values, columns=[0])

    # Create AnnData objects
    adata_RNA, adata_ATAC = get_adata(matrix, features, barcodes, label)

    # Filtering
    sc.pp.filter_cells(adata_RNA, min_genes=200)
    sc.pp.filter_genes(adata_RNA, min_cells=3)
    sc.pp.filter_cells(adata_ATAC, min_genes=200)
    sc.pp.filter_genes(adata_ATAC, min_cells=3)

    # Barcode intersection
    selected = list(set(adata_RNA.obs['barcode']) & set(adata_ATAC.obs['barcode']))
    adata_RNA = adata_RNA[[b in selected for b in adata_RNA.obs['barcode']]]
    adata_ATAC = adata_ATAC[[b in selected for b in adata_ATAC.obs['barcode']]]

    # Pseudo-bulk
    samplelist = list(set(adata_ATAC.obs['sample']))
    singlepseudobulk = adata_RNA.obs['sample'].nunique() ** 2 > 100
    TG_pseudobulk = pd.DataFrame()
    RE_pseudobulk = pd.DataFrame()

    for s in samplelist:
        r_temp = adata_RNA[adata_RNA.obs['sample'] == s]
        a_temp = adata_ATAC[adata_ATAC.obs['sample'] == s]
        tg, re = pseudo_bulk(r_temp, a_temp, singlepseudobulk)
        re[re > 100] = 100
        TG_pseudobulk = pd.concat([TG_pseudobulk, tg], axis=1)
        RE_pseudobulk = pd.concat([RE_pseudobulk, re], axis=1)

    TG_pseudobulk = TG_pseudobulk.fillna(0)
    RE_pseudobulk = RE_pseudobulk.fillna(0)

    adata_ATAC.write(os.path.join(data_dir, 'adata_ATAC.h5ad'))
    adata_RNA.write(os.path.join(data_dir, 'adata_RNA.h5ad'))

    # Save LINGER input files into isolated `data/` folder
    pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv(os.path.join(data_dir, "Peaks.txt"), header=None, index=None)
    TG_pseudobulk.to_csv(os.path.join(data_dir, "TG_pseudobulk.tsv"), sep='\t')
    RE_pseudobulk.to_csv(os.path.join(data_dir, "RE_pseudobulk.tsv"), sep='\t')
    print(f"✅ Pseudo-bulk saved: {data_dir}")

    # Run LINGER inside this directory to isolate all data/ paths
    prev_cwd = os.getcwd()
    os.chdir(run_workdir)
    try:
        preprocess(TG_pseudobulk, RE_pseudobulk, args.datadir, args.genome, args.method, args.outdir)
        LINGER_tr.training(args.datadir, args.method, args.outdir, args.activation, "human")

        # Cell population GRN
        LL_net.TF_RE_binding(args.datadir, adata_RNA, adata_ATAC, args.genome, args.method, args.outdir)
        LL_net.cis_reg(args.datadir, adata_RNA, adata_ATAC, args.genome, args.method, args.outdir)
        LL_net.trans_reg(args.datadir, args.method, args.outdir)

        # Cell-type specific GRN
        LL_net.cell_type_specific_TF_RE_binding(args.datadir, adata_RNA, adata_ATAC, args.genome, args.celltype, args.outdir)
        LL_net.cell_type_specific_cis_reg(args.datadir, adata_RNA, adata_ATAC, args.genome, args.celltype, args.outdir)
        LL_net.cell_type_specific_trans_reg(args.datadir, adata_RNA, args.celltype, args.outdir)

    finally:
        os.chdir(prev_cwd)
        print(f"✅ Finished LINGER run in: {run_workdir}")

    # Load results into MuData


if __name__ == '__main__':
    main()
