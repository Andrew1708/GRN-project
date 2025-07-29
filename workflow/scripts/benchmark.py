import os
import sys
import time
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import muon as mu
import pyranges as pr

import matplotlib.pyplot as plt

from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    roc_curve,
    precision_recall_curve
)

from typing import Literal


# Add the parent directory (assuming utils is in ../utils)
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.benchmarking_metrics import *

GRN_TOOLS = Literal["scenicplus", "celloracle", "linger"]
SCORE_COL = "score_rank"
TF2GENE_W_COL = "TF2Gene_weight"

def modality_names(grn_tool:GRN_TOOLS, cell_type_col, sample):
    if grn_tool == "scenicplus":
        if "metacell_0" in sample:
            return "scRNA_counts", "scATAC_counts", f"scRNA_counts:{cell_type_col}"
        else:
            return "scRNA_counts", "scATAC_counts", "scRNA_counts:Metacell_Key"
    if grn_tool == "celloracle":
        if "wihtout_atac" in sample:
            return "rna", "", f"rna:{cell_type_col}"
        else:
            return "rna", "atac", f"rna:{cell_type_col}"
    if grn_tool == "linger":
        return "rna", "atac", f"rna:{cell_type_col}"

def get_adata(mudata, grn_tool:GRN_TOOLS):
    if grn_tool == "scenicplus":
        return mudata["scRNA_counts"]
    elif grn_tool == "celloracle" or grn_tool == "linger":
        return mudata["rna"]
    else:
        raise ValueError(f"Unsupported GRN inference tool: {grn_tool}")


def get_mudata(path:str, grn_tool:GRN_TOOLS):
    if grn_tool == "scenicplus" or grn_tool == "celloracle" or grn_tool == "linger":
        mdata = mu.read_h5mu(path)
        return mdata
    else:
        raise ValueError(f"Unsupported GRN inference tool: {grn_tool}")

# GRNs need to have the following columns: "TF", "Gene", "Region"
def preprocess_scenicplus(scplus_mdata):
    # Extract metadata
    direct_df = pd.DataFrame(scplus_mdata.uns['direct_e_regulon_metadata'])
    extended_df = pd.DataFrame(scplus_mdata.uns['extended_e_regulon_metadata'])

    # Combine into one DataFrame
    grn = pd.concat([direct_df, extended_df], ignore_index=True)

    # Filter the relevant columns
    grn_filtered = grn[['Region', 'Gene', 'TF', 'importance_TF2G', 'importance_R2G','regulation', 'rho_TF2G', 'triplet_rank']].copy()

    # Split the 'Region' column into Chromosome, Start, End
    region_split = grn_filtered['Region'].str.extract(r'(chr[\w]+):(\d+)-(\d+)')
    region_split.columns = ['Chromosome', 'Start', 'End']

    # Convert Start and End to integers
    region_split['Start'] = region_split['Start'].astype(int)
    region_split['End'] = region_split['End'].astype(int)

    grn = pd.concat([region_split, grn_filtered], axis=1)

    max_rank = grn["triplet_rank"].max()
    min_rank = grn["triplet_rank"].min()

    grn[SCORE_COL] = (max_rank - grn["triplet_rank"]) / (max_rank - min_rank)

    raw = grn["importance_TF2G"] * grn["rho_TF2G"]
    grn[TF2GENE_W_COL] = np.tanh(3 * raw) 

    return grn

def preprocess_celloracle(mudata):
    """
    Preprocess CellOracle MuData to extract GRN matrix.
    Assumes the GRN is stored in the 'links' attribute of the MuData object.
    """

    grn = pd.DataFrame(mudata.uns['celloracle_links'])

    # Clean the -logp column
    grn["-logp"] = grn["-logp"].replace([np.inf, -np.inf], np.nan)
    grn["-logp"] = grn["-logp"].fillna(0)

    # Optional: coerce -0.0 → 0.0 for cosmetic + numeric sanity
    grn["-logp"] = grn["-logp"].apply(lambda x: 0.0 if x == 0 else x)

    # Apply ranking
    ranks = grn["-logp"].rank(method="average")  # you can also try "dense" or "first"

    # Normalize ranks to [0, 1]
    min_rank, max_rank = ranks.min(), ranks.max()
    if max_rank > min_rank:
        grn[SCORE_COL] = (ranks - min_rank) / (max_rank - min_rank)
    else:
        grn[SCORE_COL] = 0.0


    raw_score = grn["coef_mean"] * grn["-logp"]
    grn[TF2GENE_W_COL] = np.tanh(3 * raw_score)

    grn["Region"] = grn["chromosome"] + ":" + grn["start"].astype(str) + "-" + grn["end"].astype(str)

    grn = grn.rename(columns={
    "source": "TF",
    "target": "Gene",
    "chromosome": "Chromosome",
    "start": "Start",
    "end": "End"
    })

    return grn

def preprocess_linger(mudata):
    """
    Preprocesses a MuData object for Linger-based GRN inference.

    0) Extracts GRN from mudata.uns["grn"]
    1) Parses 'Region' into 'Chromosome', 'Start', 'End'
    2) Computes 'triplet_score' = |tf_re_score * tg_re_score| (min-max normalized)
    3) Computes 'tf2gene_weight' = tanh(3 * tf_tg_score)
    """
    
    # Step 0: Extract the GRN DataFrame
    grn = mudata.uns["grn"]

    # Step 1: Parse 'Region' into 'Chromosome', 'Start', 'End'
    region_split = grn["Region"].str.extract(r"^(chr[^\:]+):(\d+)-(\d+)$")
    grn["Chromosome"] = region_split[0]
    grn["Start"] = region_split[1].astype(int)
    grn["End"] = region_split[2].astype(int)

    # Step 2: Compute 'triplet_score' = |tf_re_score * tg_re_score|, normalized
    if "tf_re_score" not in grn.columns or "tg_re_score" not in grn.columns:
        raise ValueError("Expected columns 'tf_re_score' and 'tg_re_score' not found in GRN.")
    
    # Compute raw score
    grn["raw_score"] = grn["tf_re_score"] * grn["tg_re_score"] * grn["tf_tg_score"]

    # Use absolute values for ranking (if that’s your criterion)
    abs_score = grn["raw_score"].abs()

    # Create rank (method="average" handles ties gracefully)
    ranks = abs_score.rank(method="average")

    # Normalize ranks to [0,1]
    min_rank, max_rank = ranks.min(), ranks.max()
    if max_rank > min_rank:
        grn[SCORE_COL] = (ranks - min_rank) / (max_rank - min_rank)
    else:
        grn[SCORE_COL] = 0.0

    # Step 3: Compute 'tf2gene_weight' = tanh(3 * tf_tg_score)
    if "tf_tg_score" not in grn.columns:
        raise ValueError("Expected column 'tf_tg_score' not found in GRN.")
    
    grn[TF2GENE_W_COL] = np.tanh(3 * grn["tf_tg_score"])

    # fix on mudata
    adata_rna = mudata.mod["rna"]
    adata_atac = mudata.mod["atac"]

    adata_rna.obs_names = adata_rna.obs['barcode'].astype(str)
    adata_atac.obs_names = adata_atac.obs['barcode'].astype(str)
    adata_atac.var_names = adata_atac.var['gene_ids'].astype(str)
    
    return grn


def get_grn_matrix(mudata, grn_tool:GRN_TOOLS):
    if grn_tool == "scenicplus":
        return preprocess_scenicplus(mudata)
    if grn_tool == "celloracle":
        return preprocess_celloracle(mudata)
    if grn_tool == "linger":
        return preprocess_linger(mudata)


def get_benchmark_matrix(tfb_path, prt_path, frc_path, gst_path, tfm_path):
    tfb_matrix = pd.read_csv(tfb_path, sep="\t", header=None)
    tfb_matrix.columns = ["Chromosome", "Start", "End", "TF"]

    prt_matrix = pd.read_csv(prt_path, sep=',')

    frc_matrix = pd.read_csv(frc_path, index_col=0, sep=',')

    gst_matrix = pd.read_csv(gst_path)

    tfm_matrix = df = pd.read_csv(tfm_path, header=None, names=['gene'])

    return tfb_matrix, prt_matrix, frc_matrix, gst_matrix, tfm_matrix



def argparser():
    parser = argparse.ArgumentParser(description="Benchmark GRN inference tools")
    parser.add_argument("--grn_path", type=str, required=True, help="Path to the GRN inference results")
    parser.add_argument("--grn_tool", type=str, required=True, help="GRN inference tool used")
    parser.add_argument("--tfb_golden", type=str, required=True, help="Path to the golden standard TFB matrix")
    parser.add_argument("--prt_golden", type=str, required=True, help="Path to the golden standard PRT matrix")
    parser.add_argument("--frc_golden", type=str, required=True, help="Path to the golden standard FRC matrix")
    parser.add_argument("--gst_golden", type=str, required=True, help="Path to the golden standard GST matrix")
    parser.add_argument("--tfm_golden", type=str, required=True, help="Path to the golden standard TFM matrix")
    parser.add_argument("--project_name", type=str, help="Name of the project for saving results")
    parser.add_argument("--celltype_col", type=str, default="Classified_Celltype", help="Column name for cell types in the metadata")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory to save the benchmark results")
    return parser.parse_args()




def argparser():
    parser = argparse.ArgumentParser(description="Benchmark GRN inference tools")
    parser.add_argument("--grn_path", type=str, required=True, help="Path to the GRN inference results")
    parser.add_argument("--grn_tool", type=str, required=True, help="GRN inference tool used")
    parser.add_argument("--tfb_golden", type=str, required=True, help="Path to the golden standard TFB matrix")
    parser.add_argument("--prt_golden", type=str, required=True, help="Path to the golden standard PRT matrix")
    parser.add_argument("--frc_golden", type=str, required=True, help="Path to the golden standard FRC matrix")
    parser.add_argument("--gst_golden", type=str, required=True, help="Path to the golden standard GST matrix")
    parser.add_argument("--tfm_golden", type=str, required=True, help="Path to the golden standard TFM matrix")
    parser.add_argument("--project_name", type=str, help="Name of the project for saving results")
    parser.add_argument("--celltype_col", type=str, default="Classified_Celltype", help="Column name for cell types in the metadata")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory to save the benchmark results")
    return parser.parse_args()


if __name__ == "__main__":
    start_time = time.time()

    args = argparser()

    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    mudata = get_mudata(args.grn_path, args.grn_tool)
    grn_inferred = get_grn_matrix(mudata, args.grn_tool)
    adata = get_adata(mudata, args.grn_tool)

    tfb_golden, prt_golden, frc_golden, gst_matrix, tfm_matrix = get_benchmark_matrix(
        args.tfb_golden,
        args.prt_golden,
        args.frc_golden,
        args.gst_golden,
        args.tfm_golden
    )

    # Benchmark PRT
    prt_benchmark = prt_test(
        grn_inferred=grn_inferred, 
        prt_matrix=prt_golden, 
        score_column=SCORE_COL,
        weight_column=TF2GENE_W_COL,
        step=0.01 
    )

    # Benchmark FRC
    frc_benchmark = frc_test(
        grn_inferred=grn_inferred, 
        adata=adata.copy(),
        frc_matrix=frc_golden, 
        score_column=SCORE_COL,
        step=0.05
    )

    # Benchmark TFB
    tfb_benchmark = tfb_test(
        grn_inferred=grn_inferred, 
        tf_binding_matrix=tfb_golden, 
        score_column=SCORE_COL
    )

    gst_benchmark = gst_test(
        grn_inferred=grn_inferred, 
        ptw=gst_matrix, 
        rna=adata.copy(),
        score_column=SCORE_COL,
        step=0.05
    )
   



    tfm_benchmark = tfm_test(
        grn_inferred=grn_inferred, 
        db=tfm_matrix, 
        adata=adata,
        score_column=SCORE_COL,
        step=0.01
    )

    # Benchmark OMIC
    rna_mod_name, atac_mod_name, celltype_col = modality_names(args.grn_tool, args.celltype_col, args.project_name)
    # Omics TF-Gene´
    omics_tf2g = omic_test(
        grn_inferred = grn_inferred,
        mdata = mudata.copy(),
        score_column = SCORE_COL,
        source_column = "TF",
        target_column = "Gene",
        mod_source = rna_mod_name,
        mod_target= rna_mod_name,
        celltype_column= celltype_col,
        step = 0.3
    )

    if "scenicplus" == args.grn_tool or ("celloracle" == args.grn_tool and "with_atac" in args.project_name) or args.grn_tool == "linger":
        # Omics CRE-Gene
        omics_r2g = omic_test(
            grn_inferred = grn_inferred,
            mdata = mudata.copy(),
            score_column = SCORE_COL,
            source_column = "Region",
            target_column = "Gene",
            mod_source = atac_mod_name,
            mod_target= rna_mod_name,
            celltype_column= celltype_col,
            step = 0.3
        )

        omics_r2tf = omic_test(
            grn_inferred = grn_inferred,
            mdata = mudata.copy(),
            score_column = SCORE_COL,
            source_column = "Region",
            target_column = "TF",
            mod_source = atac_mod_name,
            mod_target= rna_mod_name,
            celltype_column= celltype_col,
            step = 0.3
        )



    # Assume project_name is defined
    project_name = args.project_name if args.project_name else "GRN"
    output_file = os.path.join(args.output_dir, f"{project_name}_benchmark_results.csv")

    benchmarks = []

    # Always available
    benchmarks.append(("tfb", tfb_benchmark))
    benchmarks.append(("prt", prt_benchmark))
    benchmarks.append(("frc", frc_benchmark))
    benchmarks.append(("omics_tf2g", omics_tf2g))
    benchmarks.append(("gst", gst_benchmark))
    benchmarks.append(("tfm", tfm_benchmark))

    # Optional benchmarks
    if "scenicplus" == args.grn_tool or ("celloracle" == args.grn_tool and "with_atac" in args.project_name) or args.grn_tool == "linger":
        benchmarks.append(("omics_r2g", omics_r2g))
        benchmarks.append(("omics_r2tf", omics_r2tf))

    # Define the order of columns
    columns = [
        "project_name",
        "benchmark_name",
        "tp",
        "fp",
        "fn",
        "precision",
        "recall",
        "fbeta",
        "auroc",
        "auprc",
        "best_threshold",
        "best_fbeta",
        "best_precision",
        "best_recall"
    ]

    # Initialize a list to collect rows
    rows = []
    # Iterate through each benchmark and build rows
    for name, result in benchmarks:
        row = [
            project_name,
            name,
            result["tp"],
            result["fp"],
            result["fn"],
            round(result["precision"], 3) if pd.notna(result["precision"]) else None,
            round(result["recall"], 3) if pd.notna(result["recall"]) else None,
            round(result["fbeta"], 3) if pd.notna(result["fbeta"]) else None,
            round(result["auroc"], 3) if pd.notna(result["auroc"]) else None,
            round(result["auprc"], 3) if pd.notna(result["auprc"]) else None,
            round(result["best_threshold"], 3) if pd.notna(result["best_threshold"]) else None,
            round(result["best_fbeta"], 3) if pd.notna(result["best_fbeta"]) else None,
            round(result["best_precision"], 3) if pd.notna(result["best_precision"]) else None,
            round(result["best_recall"], 3) if pd.notna(result["best_recall"]) else None
        ]
        rows.append(row)

    # Convert to a DataFrame
    benchmark_df = pd.DataFrame(rows, columns=columns)

    # Save the DataFrame to a CSV file
    benchmark_df.to_csv(output_file, sep='\t', index=False)
    print(f"Benchmark results saved to {output_file}")

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"Elapsed time: {elapsed_time:.2f} seconds")


