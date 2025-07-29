import argparse
import os
import sys
import gc
from pathlib import Path
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore")

from tqdm import tqdm

# Add the parent directory (assuming utils is in ../utils)
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.grn_utils import load_grn_file, extract_patient_time, get_grn_matrix


def list_grn_paths_by_dir(tool_dir, grn_filename, tool_name):
    grn_map = {}
    for dir_name in os.listdir(tool_dir):
        full_path = os.path.join(tool_dir, dir_name, grn_filename)
        if not os.path.exists(full_path):
            continue
        
        # if tool_name == "linger" and "scdart" not in dir_name:
        #     continue
        
        # if tool_name == "celloracle" and "MO" in dir_name:
        #     continue

        # if tool_name == "scenicplus" and "MO" in dir_name:
        #     continue

        try:
            sample_id, _ = extract_patient_time(dir_name)
        except Exception:
            print(f"[!] Skipping {dir_name}: Could not parse sample ID.")
            continue
        key = (tool_name, dir_name)
        grn_map[key] = {sample_id: full_path}
    return grn_map


def extract_all_interactions_from_df(grn_df, use_region):
    if grn_df is None or grn_df.empty:
        return set()
    if use_region:
        keys = grn_df[["Region", "TF", "Gene"]].dropna().astype(str).apply(
            lambda row: f"{row['Region']}|{row['TF']}|{row['Gene']}", axis=1)
    else:
        keys = grn_df[["TF", "Gene"]].dropna().astype(str).apply(
            lambda row: f"{row['TF']}|{row['Gene']}", axis=1)
    return set(keys)


def compute_overlap(set1, set2):
    if not set1 or not set2:
        return np.nan
    intersection = set1 & set2
    denom = min(len(set1), len(set2))
    return len(intersection) / denom if denom > 0 else 0.0


def compute_overlap_longform(grn_index, use_region, max_cache_size=5, output_csv=None, flush_every=500):
    keys = sorted(grn_index.keys())  # (tool, dir_name)
    labels = [f"{tool}+{cfg}" for tool, cfg in keys]
    n = len(keys)

    # Load existing partial results (e.g. from _temp.csv)
    done_pairs = set()
    if output_csv and os.path.exists(output_csv):
        existing_df = pd.read_csv(output_csv)
        done_pairs = set(zip(existing_df["GRN_A"], existing_df["GRN_B"]))
        print(f"Resuming from {len(done_pairs)} previously computed overlaps\n")

    total_overlaps = (n * (n + 1)) // 2
    print(f"Total overlaps to compute: {total_overlaps} (some may be cached)\n")

    future_use_count = {label: n - i - 1 for i, label in enumerate(labels)}
    grn_cache = {}
    buffer = []

    try:
        with tqdm(total=total_overlaps, desc="Computing overlaps", unit="pair") as pbar:
            for i in range(n):
                a_key = keys[i]
                a_label = labels[i]
                path_a = next(iter(grn_index[a_key].values()))
                grn_a_df = get_grn_matrix(load_grn_file(path_a), a_key[0])
                set_a = extract_all_interactions_from_df(grn_a_df, use_region)
                del grn_a_df
                gc.collect()

                for j in range(i, n):
                    b_key = keys[j]
                    b_label = labels[j]

                    if (a_label, b_label) in done_pairs:
                        pbar.update(1)
                        continue

                    if i == j:
                        coef = 1.0
                    else:
                        if b_label in grn_cache:
                            grn_b_df = grn_cache[b_label]
                        else:
                            path_b = next(iter(grn_index[b_key].values()))
                            grn_b_df = get_grn_matrix(load_grn_file(path_b), b_key[0])
                            if b_label != a_label and len(grn_cache) < max_cache_size:
                                grn_cache[b_label] = grn_b_df

                        set_b = extract_all_interactions_from_df(grn_b_df, use_region)
                        coef = compute_overlap(set_a, set_b)

                        if b_label in future_use_count:
                            future_use_count[b_label] -= 1
                            if future_use_count[b_label] <= 0:
                                grn_cache.pop(b_label, None)
                                del future_use_count[b_label]

                    buffer.append({
                        "GRN_A": a_label,
                        "GRN_B": b_label,
                        "Overlap": coef
                    })

                    if len(buffer) >= flush_every:
                        flush_df = pd.DataFrame(buffer)
                        flush_df.to_csv(output_csv, mode='a', header=not os.path.exists(output_csv), index=False)
                        print(f"Flushed {len(buffer)} rows to {output_csv}")
                        buffer.clear()

                    pbar.update(1)
    finally:
        if buffer:
            flush_df = pd.DataFrame(buffer)
            flush_df.to_csv(output_csv, mode='a', header=not os.path.exists(output_csv), index=False)
            print(f"Final flush: {len(buffer)} rows written to {output_csv}")

    return pd.read_csv(output_csv)


def parse_args():
    parser = argparse.ArgumentParser(description="Compute GRN overlap matrix from tool+config directory layout.")
    parser.add_argument("--scenicplus_dir", type=str)
    parser.add_argument("--scenicplus_filename", type=str)
    parser.add_argument("--celloracle_dir", type=str)
    parser.add_argument("--celloracle_filename", type=str)
    parser.add_argument("--linger_dir", type=str)
    parser.add_argument("--linger_filename", type=str)
    parser.add_argument("--use_region", action="store_true")
    parser.add_argument("--output_csv", type=str, required=True)
    parser.add_argument("--max_cache_size", type=int, default=25,
                        help="Max number of GRNs to cache in memory (excluding the current row GRN)")
    parser.add_argument("--flush_every", type=int, default=500,
                        help="How many rows to accumulate before flushing to disk")
    return parser.parse_args()


def main():
    args = parse_args()

    output_csv_final = Path(args.output_csv).resolve()
    output_csv_temp = output_csv_final.with_name(output_csv_final.stem + "_temp.csv")

    grn_index = {}

    if args.scenicplus_dir and args.scenicplus_filename:
        grn_index.update(
            list_grn_paths_by_dir(args.scenicplus_dir, args.scenicplus_filename, "scenicplus")
        )
    if args.celloracle_dir and args.celloracle_filename:
        grn_index.update(
            list_grn_paths_by_dir(args.celloracle_dir, args.celloracle_filename, "celloracle")
        )
    if args.linger_dir and args.linger_filename:
        grn_index.update(
            list_grn_paths_by_dir(args.linger_dir, args.linger_filename, "linger")
        )

    overlap_df = compute_overlap_longform(
        grn_index,
        args.use_region,
        args.max_cache_size,
        output_csv=output_csv_temp,
        flush_every=args.flush_every
    )

    output_csv_temp.rename(output_csv_final)
    print(f"Overlap results finalized at: {output_csv_final}")


if __name__ == "__main__":
    main()
