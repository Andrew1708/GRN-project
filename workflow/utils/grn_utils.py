import os
import re
import pickle
import pandas as pd
import numpy as np
import anndata
import mudata

SCORE_COL = "Score"

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
