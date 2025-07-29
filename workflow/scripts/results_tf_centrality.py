import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.preprocessing import MinMaxScaler

# Add parent directory to path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.grn_utils import load_grn_file, get_grn_matrix

def load_binary_tf_gene_edges(grn_df):
    """Returns edge list as list of (TF, Gene) tuples."""
    grn_df = grn_df.dropna(subset=["TF", "Gene"]).astype(str)
    return list(zip(grn_df["TF"], grn_df["Gene"]))

def build_directed_graph(edges):
    """Returns a directed NetworkX graph from (TF, Gene) edges."""
    G = nx.DiGraph()
    G.add_edges_from(edges)
    return G

def compute_tf_centralities(G):
    """Computes centralities for nodes with out-degree > 0 (i.e., TFs)."""
    tf_nodes = [n for n in G.nodes if G.out_degree(n) > 0]

    degree_out = dict(G.out_degree(tf_nodes))
    degree_in = dict(G.in_degree(tf_nodes))
    closeness = nx.closeness_centrality(G)
    betweenness = nx.betweenness_centrality(G, normalized=True)
    try:
        eigen = nx.eigenvector_centrality(G, max_iter=500)
    except nx.PowerIterationFailedConvergence:
        eigen = {n: 0.0 for n in G.nodes}
    pagerank = nx.pagerank(G)

    df = pd.DataFrame({
        "TF": tf_nodes,
        "out_degree": [degree_out.get(n, 0) for n in tf_nodes],
        "in_degree": [degree_in.get(n, 0) for n in tf_nodes],
        "closeness": [closeness.get(n, 0) for n in tf_nodes],
        "betweenness": [betweenness.get(n, 0) for n in tf_nodes],
        "eigenvector": [eigen.get(n, 0) for n in tf_nodes],
        "pagerank": [pagerank.get(n, 0) for n in tf_nodes],
    })

    return normalize_and_score(df)

def normalize_and_score(df):
    """Normalizes columns and computes average score per TF."""
    metric_cols = df.columns.difference(["TF"])
    scaler = MinMaxScaler()
    df[metric_cols] = scaler.fit_transform(df[metric_cols])
    df["average_score"] = df[metric_cols].mean(axis=1)
    return df

def process_grn(grn_path, tool, top_n=None, score_col=None):
    try:
        mdata = load_grn_file(grn_path)
        grn_df = get_grn_matrix(mdata, tool)

        if top_n and score_col:
            if score_col in grn_df.columns:
                grn_df = grn_df.sort_values(by=score_col, ascending=False).head(top_n)
            else:
                print(f"[!] Column '{score_col}' not found in GRN. Skipping top-N filter.")

        edges = load_binary_tf_gene_edges(grn_df)
        G = build_directed_graph(edges)
        return compute_tf_centralities(G)
    except Exception as e:
        print(f"[!] Failed to process {grn_path}: {e}")
        return None

def parse_args():
    parser = argparse.ArgumentParser(description="Compute TF centrality metrics from a single GRN.")
    parser.add_argument("--grn_path", type=str, required=True,
                        help="Path to a single GRN file (CSV inside a config folder).")
    parser.add_argument("--tool", type=str, required=True,
                        help="Name of tool that produced the GRN (e.g., scenicplus, celloracle, linger).")
    parser.add_argument("--output_csv", type=str, required=True,
                        help="Where to save the centrality scores CSV.")
    parser.add_argument("--top_n", type=int, default=None,
                        help="If set, limit GRN to top-N rows based on --score_col.")
    parser.add_argument("--score_col", type=str, default=None,
                        help="Column to use for top-N filtering (e.g. 'score', 'importance').")
    return parser.parse_args()

def main():
    args = parse_args()
    print(f"🔍 Processing: {args.grn_path}")
    df = process_grn(args.grn_path, args.tool, top_n=args.top_n, score_col=args.score_col)
    if df is not None:
        df.to_csv(args.output_csv, index=False)
        print(f"✅ Saved centrality scores to: {args.output_csv}")

if __name__ == "__main__":
    main()
