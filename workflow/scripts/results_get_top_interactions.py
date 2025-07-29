import argparse
import os
import re
import pickle
import anndata
import mudata
import numpy as np
import pandas as pd

SCORE_COL = "Score"
FINAL_COLUMNS = ["Tool", "Patient", "Time_Point", "TF", "Gene", "CRE", "Rank"]

# ------------------- GRN Preprocessors -------------------

def preprocess_scenicplus(scplus_mdata):
    direct_df = pd.DataFrame(scplus_mdata.uns['direct_e_regulon_metadata'])
    extended_df = pd.DataFrame(scplus_mdata.uns['extended_e_regulon_metadata'])
    grn = pd.concat([direct_df, extended_df], ignore_index=True)
    grn_filtered = grn[['Region', 'Gene', 'TF', 'importance_TF2G', 'importance_R2G', 'regulation', 'rho_TF2G', 'triplet_rank']].copy()
    region_split = grn_filtered['Region'].str.extract(r'(chr[\w]+):(\d+)-(\d+)')
    region_split.columns = ['Chromosome', 'Start', 'End']
    region_split['Start'] = region_split['Start'].astype(int)
    region_split['End'] = region_split['End'].astype(int)
    grn = pd.concat([region_split, grn_filtered], axis=1)
    max_rank = grn["triplet_rank"].max()
    min_rank = grn["triplet_rank"].min()
    grn[SCORE_COL] = (max_rank - grn["triplet_rank"]) / (max_rank - min_rank)
    return grn

def preprocess_celloracle(mudata):
    grn = pd.DataFrame(mudata.uns['celloracle_links'])
    grn["-logp"] = grn["-logp"].replace([np.inf, -np.inf], np.nan).fillna(0)
    grn["-logp"] = grn["-logp"].apply(lambda x: 0.0 if x == 0 else x)
    min_logp = grn["-logp"].min()
    max_logp = grn["-logp"].max()
    grn[SCORE_COL] = (grn["-logp"] - min_logp) / (max_logp - min_logp) if max_logp != min_logp else 0.0
    grn["Region"] = grn["chromosome"] + ":" + grn["start"].astype(str) + "-" + grn["end"].astype(str)
    grn = grn.rename(columns={"source": "TF", "target": "Gene", "chromosome": "Chromosome", "start": "Start", "end": "End"})
    return grn

def preprocess_linger(mudata):
    grn = mudata.uns["grn"]
    region_split = grn["Region"].str.extract(r"^(chr[^\:]+):(\d+)-(\d+)$")
    grn["Chromosome"] = region_split[0]
    grn["Start"] = region_split[1].astype(int)
    grn["End"] = region_split[2].astype(int)
    if "tf_re_score" not in grn.columns or "tg_re_score" not in grn.columns:
        raise ValueError("Missing 'tf_re_score' or 'tg_re_score'")
    abs_score = (grn["tf_re_score"] * grn["tg_re_score"]).abs()
    grn[SCORE_COL] = (abs_score - abs_score.min()) / (abs_score.max() - abs_score.min()) if abs_score.max() > abs_score.min() else 0.0
    return grn

def get_grn_matrix(mdata, grn_tool):
    if grn_tool == "scenicplus":
        return preprocess_scenicplus(mdata)
    elif grn_tool == "celloracle":
        return preprocess_celloracle(mdata)
    elif grn_tool == "linger":
        return preprocess_linger(mdata)
    else:
        raise ValueError(f"Unsupported GRN tool: {grn_tool}")

# ------------------- I/O Helpers -------------------

def load_grn_file(path):
    if path.endswith(".h5mu"):
        return mudata.read(path)
    elif path.endswith(".h5ad"):
        return anndata.read_h5ad(path)
    elif path.endswith(".pkl"):
        with open(path, "rb") as f:
            return pickle.load(f)
    else:
        raise ValueError(f"Unsupported file format: {path}")

def extract_patient_time(folder_name):
    match = re.match(r"([^_]+)_([^_]+)_.*", folder_name)
    if not match:
        raise ValueError(f"Folder name {folder_name} doesn't match expected pattern <Patient>_<TimePoint>_...")
    return match.group(1), match.group(2)

def process_tool(tool_name, base_dir, grn_filename, folder_filter, top_n):
    results = []
    for folder in os.listdir(base_dir):
        if folder_filter and folder_filter not in folder:
            continue
        full_path = os.path.join(base_dir, folder, grn_filename)
        if not os.path.exists(full_path):
            print(f"  [!] Skipping '{folder}': File not found at {full_path}")
            continue
        try:
            patient, time_point = extract_patient_time(folder)
            mdata = load_grn_file(full_path)
            grn_df = get_grn_matrix(mdata, tool_name)
            top_df = grn_df.nlargest(top_n, SCORE_COL).copy()
            top_df["Rank"] = range(1, len(top_df)+1)
            subset = top_df[["TF", "Gene", "Region", "Rank"]].copy()
            subset.columns = ["TF", "Gene", "CRE", "Rank"]
            subset["Tool"] = tool_name
            subset["Patient"] = patient
            subset["Time_Point"] = time_point
            results.append(subset[FINAL_COLUMNS])
        except Exception as e:
            print(f"  [x] Error processing '{folder}': {e}")
    return results

# ------------------- CLI -------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate top-N TF–Gene interactions from GRNs.")
    parser.add_argument("--scenicplus_dir", type=str)
    parser.add_argument("--scenicplus_filename", type=str)
    parser.add_argument("--scenicplus_filter", type=str, default="")
    parser.add_argument("--celloracle_dir", type=str)
    parser.add_argument("--celloracle_filename", type=str)
    parser.add_argument("--celloracle_filter", type=str, default="")
    parser.add_argument("--linger_dir", type=str)
    parser.add_argument("--linger_filename", type=str)
    parser.add_argument("--linger_filter", type=str, default="")
    parser.add_argument("--top_n", type=int, default=100)
    parser.add_argument("--output_csv", type=str, help="Optional path to save the final CSV.")
    return parser.parse_args()

def main():
    args = parse_args()
    all_results = []

    if args.scenicplus_dir and args.scenicplus_filename:
        all_results += process_tool("scenicplus", args.scenicplus_dir, args.scenicplus_filename, args.scenicplus_filter, args.top_n)
    if args.celloracle_dir and args.celloracle_filename:
        all_results += process_tool("celloracle", args.celloracle_dir, args.celloracle_filename, args.celloracle_filter, args.top_n)
    if args.linger_dir and args.linger_filename:
        all_results += process_tool("linger", args.linger_dir, args.linger_filename, args.linger_filter, args.top_n)

    final_df = pd.concat(all_results, ignore_index=True)
    print(f"\n✅ Total interactions collected: {len(final_df)}")

    if args.output_csv:
        final_df.to_csv(args.output_csv, index=False)
        print(f"💾 Saved to {args.output_csv}")

if __name__ == "__main__":
    main()
