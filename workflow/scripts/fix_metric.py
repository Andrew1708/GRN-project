import os
import sys
import time
import argparse
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import muon as mu

# Add the parent directory (assuming utils is in ../utils)
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.benchmarking_metrics import *

GRN_TOOLS = Literal["scenicplus", "celloracle"]
SCORE_COL = "score_rank"
TF2GENE_W_COL = "TF2Gene_weight"

def modality_names(grn_tool: GRN_TOOLS, cell_type_col, sample):
    if grn_tool == "scenicplus":
        if "metacell_0" in sample:
            return "scRNA_counts", "scATAC_counts", f"scRNA_counts:{cell_type_col}"
        else:
            return "scRNA_counts", "scATAC_counts", "scRNA_counts:Metacell_Key"
    if grn_tool == "celloracle":
        if "without_atac" in sample:
            return "rna", "", f"rna:{cell_type_col}"
        else:
            return "rna", "atac", f"rna:{cell_type_col}"

def get_adata(mudata, grn_tool: GRN_TOOLS):
    return mudata["scRNA_counts"] if grn_tool == "scenicplus" else mudata["rna"]

def get_mudata(path: str, grn_tool: GRN_TOOLS):
    return mu.read_h5mu(path)

def preprocess_scenicplus(scplus_mdata):
    direct_df = pd.DataFrame(scplus_mdata.uns['direct_e_regulon_metadata'])
    extended_df = pd.DataFrame(scplus_mdata.uns['extended_e_regulon_metadata'])
    grn = pd.concat([direct_df, extended_df], ignore_index=True)
    grn = grn[['Region', 'Gene', 'TF', 'importance_TF2G', 'importance_R2G', 'regulation', 'rho_TF2G', 'triplet_rank']]
    region_split = grn['Region'].str.extract(r'(chr[\w]+):(\d+)-(\d+)')
    region_split.columns = ['Chromosome', 'Start', 'End']
    region_split['Start'] = region_split['Start'].astype(int)
    region_split['End'] = region_split['End'].astype(int)
    grn = pd.concat([region_split, grn], axis=1)
    grn[SCORE_COL] = (grn["triplet_rank"].max() - grn["triplet_rank"]) / (grn["triplet_rank"].max() - grn["triplet_rank"].min())
    grn[TF2GENE_W_COL] = np.tanh(3 * grn["importance_TF2G"] * grn["rho_TF2G"])
    return grn

def preprocess_celloracle(mudata):
    grn = pd.DataFrame(mudata.uns['celloracle_links'])
    grn[SCORE_COL] = (grn["-logp"] - grn["-logp"].min()) / (grn["-logp"].max() - grn["-logp"].min())
    grn[TF2GENE_W_COL] = np.tanh(3 * grn["coef_mean"] * grn["-logp"])
    return grn.rename(columns={"source": "TF", "target": "Gene", "chromosome": "Chromosome", "start": "Start", "end": "End"})

def get_grn_matrix(mudata, grn_tool: GRN_TOOLS):
    return preprocess_scenicplus(mudata) if grn_tool == "scenicplus" else preprocess_celloracle(mudata)

def get_benchmark_matrix(tfb_path, prt_path, frc_path, gst_path, tfm_path):
    tfb_matrix = pd.read_csv(tfb_path, sep="\t", header=None, names=["Chromosome", "Start", "End", "TF"])
    prt_matrix = pd.read_csv(prt_path, sep=',')
    frc_matrix = pd.read_csv(frc_path, index_col=0, sep=',')
    gst_matrix = pd.read_csv(gst_path)
    tfm_matrix = pd.read_csv(tfm_path, header=None, names=['gene'])
    return tfb_matrix, prt_matrix, frc_matrix, gst_matrix, tfm_matrix

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark_table", type=str, required=True, help="Existing benchmark results file to update")
    parser.add_argument("--tests", type=str, required=True, help="Comma-separated list of tests to rerun (e.g. tfb,gst,omics_r2g)")
    parser.add_argument("--grn_path", type=str, required=True)
    parser.add_argument("--grn_tool", type=str, required=True)
    parser.add_argument("--tfb_golden", type=str, required=True)
    parser.add_argument("--prt_golden", type=str, required=True)
    parser.add_argument("--frc_golden", type=str, required=True)
    parser.add_argument("--gst_golden", type=str, required=True)
    parser.add_argument("--tfm_golden", type=str, required=True)
    parser.add_argument("--project_name", type=str, required=True)
    parser.add_argument("--celltype_col", type=str, default="Classified_Celltype")
    parser.add_argument("--out_file", type=str, default=None, help="Optional: path to write corrected output table")
    return parser.parse_args()

if __name__ == "__main__":
    args = argparser()
    start_time = time.time()

    test_names = args.tests.split(",")

    mudata = get_mudata(args.grn_path, args.grn_tool)
    adata = get_adata(mudata, args.grn_tool)
    grn = get_grn_matrix(mudata, args.grn_tool)

    tfb_golden, prt_golden, frc_golden, gst_matrix, tfm_matrix = get_benchmark_matrix(
        args.tfb_golden, args.prt_golden, args.frc_golden, args.gst_golden, args.tfm_golden
    )

    rna_mod, atac_mod, celltype = modality_names(args.grn_tool, args.celltype_col, args.project_name)

    test_dispatch = {
        "tfb": lambda: tfb_test(grn, tfb_golden, SCORE_COL),
        "prt": lambda: prt_test(grn, prt_golden, SCORE_COL, TF2GENE_W_COL, step=0.01),
        "frc": lambda: frc_test(grn, adata.copy(), frc_golden, SCORE_COL, step=0.05),
        "gst": lambda: gst_test(grn, gst_matrix, adata.copy(), SCORE_COL, step=0.05),
        "tfm": lambda: tfm_test(grn, tfm_matrix, adata, SCORE_COL, step=0.01),
        "omics_tf2g": lambda: omic_test(grn, mudata.copy(), "TF", "Gene", SCORE_COL, rna_mod, rna_mod, celltype, step=0.3),
        "omics_r2g": lambda: omic_test(grn, mudata.copy(), "Region", "Gene", SCORE_COL, atac_mod, rna_mod, celltype, step=0.3),
        "omics_r2tf": lambda: omic_test(grn, mudata.copy(), "Region", "TF", SCORE_COL, atac_mod, rna_mod, celltype, step=0.3),
    }

    new_results = []
    for test in test_names:
        if test not in test_dispatch:
            print(f"[WARN] Unknown test: {test}. Skipping.")
            continue
        print(f"[INFO] Rerunning benchmark test: {test}")
        res = test_dispatch[test]()
        res.update({"benchmark_name": test, "project_name": args.project_name})
        new_results.append(res)

    updated_df = pd.read_csv(args.benchmark_table, sep="\t")
    updated_df = updated_df[~updated_df["benchmark_name"].isin(test_names) | (updated_df["project_name"] != args.project_name)]

    new_df = pd.DataFrame(new_results)
    final_df = pd.concat([updated_df, new_df], ignore_index=True)

    # Choose where to save
    save_path = args.out_file if args.out_file else args.benchmark_table
    # Round numeric columns to 3 decimals
    final_df = final_df.round(3)
    final_df.to_csv(save_path, sep="\t", index=False)

    print(f"[INFO] Updated benchmark table saved to {save_path}")
    print(f"Total time: {time.time() - start_time:.2f} sec")
