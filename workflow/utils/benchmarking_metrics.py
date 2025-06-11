import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

import scipy
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    roc_curve,
    precision_recall_curve
)
from sklearn.model_selection import train_test_split
from xgboost import XGBRegressor

import pyranges as pr
import anndata as ad
import mudata as mu
import muon as mu

import decoupler as dc
from celloracle import Oracle
import celloracle.trajectory.oracle_utility as co

from typing import Callable



########################################
#               GENERAL                #
########################################
def format_benchmarking_metrics(
    tp: int,
    fp: int,
    fn: int,
    precision: float,
    recall: float,
    fbeta: float,
    auroc: float,
    auprc: float,
    best_threshold: float,
    best_fbeta: float,
    best_precision: float,
    best_recall: float
) -> dict:
    """
    Package benchmarking metrics into a standardized dictionary.

    Parameters
    ----------
    tp, fp, fn : int
        Counts of true positives, false positives, and false negatives.
    precision, recall, fbeta : float
        Base precision, recall, and F-beta score.
    auroc, auprc : float
        Area under ROC and precision-recall curves.
    best_threshold : float
        Threshold yielding the best F-beta from PR curve.
    best_fbeta, best_precision, best_recall : float
        Best metrics at optimal threshold.

    Returns
    -------
    dict
        Standardized results dictionary.
    """
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "fbeta": fbeta,
        "auroc": auroc,
        "auprc": auprc,
        "best_threshold": best_threshold,
        "best_fbeta": best_fbeta,
        "best_precision": best_precision,
        "best_recall": best_recall
    }


def compute_auprc(precisions, recalls):
    """
    Compute the Area Under the Precision-Recall Curve (AUPRC)
    using the trapezoidal rule.

    Parameters
    ----------
    precisions : list or np.ndarray
        Sequence of precision values corresponding to recall points.
        Must be ordered according to increasing recall.
    recalls : list or np.ndarray
        Sequence of recall values, in increasing order.

    Returns
    -------
    float
        The estimated area under the precision-recall curve (AUPRC),
        computed using the trapezoidal approximation.

    Notes
    -----
    This function approximates the AUPRC by summing the area of trapezoids
    defined by successive points on the precision-recall curve. It assumes
    that both `precisions` and `recalls` are already sorted by recall in
    ascending order and are of the same length.

    This implementation does not perform any interpolation or sorting.
    """
    auprc = 0.0
    for i in range(1, len(recalls)):
        delta_recall = recalls[i] - recalls[i - 1]
        avg_precision = (precisions[i] + precisions[i - 1]) / 2
        auprc += delta_recall * avg_precision
    return auprc

def threshold_benchmarking(
    grn: pd.DataFrame,
    benchmark_func: Callable,
    benchmark_kwargs: dict,
    score_column: str,
    step: float = 0.01
) -> dict:
    """
    Evaluate an inferred GRN across thresholds using a provided benchmarking function.

    Parameters
    ----------
    grn : pd.DataFrame
        Inferred gene regulatory network.
    benchmark_func : Callable
        Benchmark function that returns either 6 or 3 values (with or without TP/FP/FN).
    benchmark_kwargs : dict
        Keyword arguments to pass to the benchmark function.
    score_column : str
        Name of the column to threshold (e.g., 'weight').
    step : float, optional
        Threshold step size (default is 0.01).

    Returns
    -------
    dict
        Benchmark metrics including AUPRC and best-threshold F-beta.
    """

    # Full network benchmark
    benchmark_result = benchmark_func(**benchmark_kwargs)
    tp0, fp0, fn0, precision0, recall0, fbeta_score0 = benchmark_result


    precisions = []
    recalls = []
    fbeta_scores = []
    thresholds = []

    for threshold in np.arange(step, 1.01, step):
        print(f"Threshold: {threshold:.2f}")
        filtered_grn = grn[abs(grn[score_column]) >= threshold].copy()

        if filtered_grn.empty:
            print("Filtered GRN is empty. Skipping...")
            continue

        try:
            benchmark_kwargs['grn'] = filtered_grn
            result = benchmark_func(**benchmark_kwargs)
            _, _, _, precision, recall, fbeta_score = result

        except ValueError as e:
            print(f"Threshold {threshold:.2f} skipped: {e}")
            continue
        
        print(f"Precision: {precision:.4f}, Recall: {recall:.4f}, F-beta: {fbeta_score:.4f}")
        precisions.append(precision)
        recalls.append(recall)
        fbeta_scores.append(fbeta_score)
        thresholds.append(threshold)

    if len(precisions) == 0 or len(recalls) == 0:
        best_precision = precision0
        best_recall = recall0
        best_fbeta = fbeta_score0
        best_threshold = 0.0
        auprc = float("nan")
    else:
        recall_precision_pairs = sorted(zip(recalls, precisions))
        recalls_sorted, precisions_sorted = zip(*recall_precision_pairs)
        auprc = compute_auprc(precisions_sorted, recalls_sorted)

        best_idx = np.argmax(fbeta_scores)
        best_threshold = thresholds[best_idx]
        best_precision = precisions[best_idx]
        best_recall = recalls[best_idx]
        best_fbeta = fbeta_scores[best_idx]

    return format_benchmarking_metrics(
        tp=tp0,
        fp=fp0,
        fn=fn0,
        precision=precision0,
        recall=recall0,
        fbeta=fbeta_score0,
        auroc=float('nan'),
        auprc=auprc,
        best_threshold=best_threshold,
        best_fbeta=best_fbeta,
        best_precision=best_precision,
        best_recall=best_recall
    )

def f_beta_score(prc, rcl, beta=0.1):
    if prc + rcl == 0:
        return 0
    return (1 + beta**2) * (prc * rcl) / ((prc * beta**2) + rcl)

#########################################
#               TF-BINDING              #
#########################################

def tfb_test(
    grn_inferred: pd.DataFrame,
    tf_binding_matrix: pd.DataFrame,
    score_column: str,
    beta: float = 0.1,
    slack: int = 100
) -> dict:
    """
    Evaluate an inferred TF–region regulatory network against a gold-standard
    TF–region binding matrix using overlap, F-beta, and ROC/PR metrics.

    Parameters
    ----------
    grn_inferred : pd.DataFrame
        DataFrame with inferred TF–region interactions. Must contain columns:
        'Chromosome', 'Start', 'End', 'TF', and a score column.
    tf_binding_matrix : pd.DataFrame
        Gold-standard TF–region binding matrix. Must contain:
        'Chromosome', 'Start', 'End', 'TF'.
    score_column : str
        Name of the column in `grn_inferred` to use as prediction score.
    beta : float, optional
        Beta parameter for F-beta score (default is 0.1).
    slack : int, optional
        Allowed slack (bp) in overlap between predicted and true regions (default is 100).

    Returns
    -------
    dict
        Dictionary with precision, recall, F-beta, AUROC, AUPRC,
        and best threshold/metrics from the PR curve.
    """

    # Step 0: Filter TFs to shared ones
    shared_tfs = set(grn_inferred["TF"]).intersection(set(tf_binding_matrix["TF"]))
    tf_binding_matrix_filtered = tf_binding_matrix[tf_binding_matrix["TF"].isin(shared_tfs)]

    # Step 1: Convert to PyRanges
    grn = pr.PyRanges(grn_inferred[["Chromosome", "Start", "End", "TF"]])
    golden = pr.PyRanges(tf_binding_matrix_filtered[["Chromosome", "Start", "End", "TF"]])

    # Step 2: Overlap predicted and true regions with slack
    overlap = grn.join(golden, strandedness=False, slack=slack, suffix="_golden")

    # Step 3: True Positives = matching TF identity
    overlap_df = overlap.df
    match_df = overlap_df[overlap_df["TF"] == overlap_df["TF_golden"]]
    tp_df = match_df[["Chromosome", "Start", "End", "TF"]].drop_duplicates()

    # Step 4: Count TP, FP, FN
    grn_df = grn.df.drop_duplicates(subset=["Chromosome", "Start", "End", "TF"])
    golden_df = golden.df.drop_duplicates(subset=["Chromosome", "Start", "End", "TF"])

    tp = len(tp_df)
    fp = len(grn_df) - tp
    fn = len(golden_df) - tp

    # Step 5: Base metrics
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    fbeta_score = (1 + beta**2) * (precision * recall) / (beta**2 * precision + recall + 1e-10)

    # Step 6: Label inferred GRN regions for ROC/PR evaluation
    df = grn_inferred[["Chromosome", "Start", "End", "TF", score_column]].copy()
    df["_merge"] = df.merge(
        tp_df,
        on=["Chromosome", "Start", "End", "TF"],
        how="left",
        indicator=True
    )["_merge"]
    df["label"] = (df["_merge"] == "both").astype(int)

    y_true = df["label"]
    y_score = df[score_column]

    # Step 7: Remove NaNs from score
    valid = y_score.notna()
    y_true = y_true[valid]
    y_score = y_score[valid]

    # Step 8: Compute AUROC and AUPRC
    if y_true.nunique() < 2:
        auroc = float("nan")
        auprc = float("nan")
    else:
        auroc = roc_auc_score(y_true, y_score)
        auprc = average_precision_score(y_true, y_score)

    # Step 9: Precision-recall curve and optimal F-beta
    precision_curve, recall_curve, thresholds = precision_recall_curve(y_true, y_score)
    precision_curve = precision_curve[1:]
    recall_curve = recall_curve[1:]
    thresholds = thresholds[:len(precision_curve)]

    f_beta = (1 + beta**2) * (precision_curve * recall_curve) / (
        beta**2 * precision_curve + recall_curve + 1e-10
    )

    best_idx = f_beta.argmax()
    best_thresh = thresholds[best_idx] if len(thresholds) > 0 else None
    best_fbeta = f_beta[best_idx] if len(f_beta) > 0 else 0.0
    best_precision = precision_curve[best_idx] if len(precision_curve) > 0 else 0.0
    best_recall = recall_curve[best_idx] if len(recall_curve) > 0 else 0.0

    # Step 10: Return results
    return format_benchmarking_metrics(
        tp=tp,
        fp=fp,
        fn=fn,
        precision=precision,
        recall=recall,
        fbeta=fbeta_score,
        auroc=auroc,
        auprc=auprc,
        best_threshold=best_thresh,
        best_fbeta=best_fbeta,
        best_precision=best_precision,
        best_recall=best_recall
    )


#########################################
#               TF-ACTIVITY             #
#########################################

def tf_activity_benchmark(
    grn: pd.DataFrame,
    golden: pd.DataFrame,
    weight_column: str,
    beta: float = 0.1,
    q_fdr: float = 0.05
) -> tuple:
    """
    Benchmark an inferred GRN using TF perturbation data and decoupler's ULM method.

    This function evaluates whether the inferred regulatory interactions in `grn`
    are able to recover known transcription factor (TF) perturbation effects from
    a matrix of gene expression changes following TF perturbation experiments.

    The ULM method (from decoupler-py) is used to score TF activity per experiment.
    A TF perturbation is considered a true positive if:
      - the perturbed TF is present in the GRN (i.e., used as a source),
      - the enrichment score from ULM is positive, and
      - the Benjamini-Hochberg corrected FDR is below `q_fdr`.

    Parameters
    ----------
    grn : pd.DataFrame
        Inferred gene regulatory network with at least the following columns:
        'source' (TF), 'target' (gene), and 'weight' (interaction score).
    golden : pd.DataFrame
        TF perturbation matrix. Must include:
        - 'Dataset ID': experiment identifier,
        - 'TF_name': name of the perturbed TF,
        - gene expression logFC values as columns (one per gene).
    beta : float, optional
        Beta value used for computing the F-beta score (default is 0.1).
    q_fdr : float, optional
        FDR threshold used to call an enrichment as significant (default is 0.1).

    Returns
    -------
    tuple
        A tuple containing:
        - tp : int
            Number of true positives (perturbed TFs recovered with positive enrichment and significant FDR).
        - fp : int
            Number of false positives (perturbed TFs tested but not meeting TP criteria).
        - fn : int
            Number of false negatives (perturbed TFs absent from the GRN).
        - precision : float
            Precision of recovered TFs.
        - recall : float
            Recall of recovered TFs.
        - fbeta_score : float
            F-beta score summarizing precision and recall.
    """
    # Sort so that the lowest triplet_rank comes first
    grn = grn.copy()
    grn = grn.rename(columns={'TF': 'source', 'Gene': 'target', weight_column: "weight"})
    grn = grn.sort_values('triplet_rank', ascending=True)

    # Drop duplicates, keeping the best (lowest-ranked) edge per source-target pair
    grn = grn.drop_duplicates(subset=['source', 'target'], keep='first')

    # Prepare infered GRN
    grn['source'] = grn['source'].str.upper()
    grn['target'] = grn['target'].str.upper()

    expr_df = golden.copy()
    # === Step 1: Prepare expression matrix ===
    expr_df = expr_df.set_index("Dataset ID")
    expr_df.columns = expr_df.columns.str.upper()

    # Prepare TF map with uppercase TF names
    expr_tf_map = golden[['Dataset ID', 'TF_name']].copy()
    expr_tf_map.columns = expr_tf_map.columns.str.upper()

    # === Step 2: Filter expression matrix to TFs in GRN ===
    network_tfs = set(grn['source'])
    network_genes = set(grn['target'])

    expr_filtered = expr_df[expr_df['TF_NAME'].isin(network_tfs)].copy()
    expr_matrix = expr_filtered.drop(columns=['TF_NAME'])

    expr_genes = set(expr_matrix.columns)
    shared_targets = expr_genes & network_genes
    shared_sources = expr_genes & network_tfs
    common_genes = sorted(shared_targets | shared_sources)

    expr_matrix = expr_filtered[common_genes]
    expr_tf_map_filtered = expr_tf_map[expr_tf_map['DATASET ID'].isin(expr_matrix.index)].copy()

    try:
        estimate, pvals = dc.mt.ulm(
            data=expr_matrix,
            net=grn,
            tmin=0,
            verbose=False
        )
    except AssertionError as e:
        print(f"[WARNING] Skipping TF-activity benchmarking for this dataset. Reason: {e}")
        return (0, 0, 0, 0, 0, 0)

    # === Step 5: BH-correct p-values row-wise ===
    adj_pvals = pd.DataFrame(index=pvals.index, columns=pvals.columns)
    for idx in pvals.index:
        raw_p = pvals.loc[idx].values.astype(float)
        _, fdr_vals, _, _ = multipletests(raw_p, method='fdr_bh')
        adj_pvals.loc[idx] = fdr_vals

    # === Step 6: Count true positives (score > 0 and FDR < threshold) ===
    tp_results = []
    fp = 0
    for _, row in expr_tf_map_filtered.iterrows():
        exp = row['DATASET ID']
        tf = row['TF_NAME']
        if tf not in estimate.columns:
            continue
        score = estimate.loc[exp, tf]
        fdr = adj_pvals.loc[exp, tf]
        if score < 0 and fdr <= q_fdr:
            tp_results.append({
                "experiment_id": exp,
                "TF": tf,
                "score": score,
                "FDR": fdr
            })
        else:
            fp += 1

    tp_df = pd.DataFrame(tp_results)

    # === Step 7: Compute precision/recall ===
    tp = len(tp_df)

    perturbed_tfs = set(expr_tf_map['TF_NAME'].str.upper())
    diff_tf = perturbed_tfs - network_tfs
    fn = len(diff_tf)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    fbeta_score = (1 + beta**2) * (precision * recall) / (beta**2 * precision + recall) if (precision + recall) > 0 else 0

    return tp, fp, fn, precision, recall, fbeta_score



def prt_test(
    grn_inferred: pd.DataFrame,
    prt_matrix: pd.DataFrame,
    weight_column: str,
    score_column: str,
    beta: float = 0.1,
    q_fdr: float = 0.05,
    step: float = 0.01,
) -> dict:
    """
    Benchmark a GRN's ability to recover TF perturbation effects using enrichment scoring.

    This function evaluates whether the inferred gene regulatory network (GRN) correctly
    identifies perturbed transcription factors (TFs) based on perturbation expression
    assays. It applies decoupler-py's ULM method via `tf_activity_benchmark` and 
    systematically filters the GRN across thresholds based on interaction strength 
    (absolute value of `score_column`).

    For each threshold, the GRN is filtered and evaluated against the perturbation
    matrix. A true positive (TP) is defined as a TF that:
      - is present in the GRN,
      - has a positive enrichment score,
      - and passes a Benjamini-Hochberg adjusted FDR cutoff (`q_fdr`).

    The function summarizes overall performance (on the full GRN), and computes:
      - precision, recall, F-beta
      - best threshold by F-beta
      - area under the precision-recall curve (AUPRC)

    Parameters
    ----------
    grn_inferred : pd.DataFrame
        Inferred gene regulatory network with at least the columns:
        ['source', 'target', score_column].
    prt_matrix : pd.DataFrame
        TF perturbation expression matrix. Must include 'Dataset ID', 'TF_name',
        and gene expression values.
    score_column : str
        Column in `grn_inferred` representing interaction strength.
    beta : float, optional
        Beta value for F-beta score (default is 0.1).
    q_fdr : float, optional
        FDR threshold for considering a TF activity significant (default is 0.05).
    step : float, optional
        Step size for thresholding (default is 0.01, yields ~100 thresholds).

    Returns
    -------
    dict
        Dictionary of benchmarking metrics containing:
        - tp, fp, fn
        - precision, recall, fbeta
        - auprc (area under PR curve)
        - best_threshold, best_fbeta, best_precision, best_recall
    """

    return threshold_benchmarking(
        grn=grn_inferred,
        benchmark_func=tf_activity_benchmark,
        benchmark_kwargs={
            "grn":grn_inferred,
            "weight_column":weight_column,
            "golden": prt_matrix,
            "beta": beta,
            "q_fdr": q_fdr
        },
        score_column=score_column,
        step=step
    )


##########################################
#               FORECASTING              #
##########################################

def init_celloracle(adata, grn, fit_grn):
    oracle = Oracle()
    oracle.adata = adata
    oracle.adata.obsm['X_umap'] = np.zeros((adata.shape[0], 2))
    oracle.adata.layers['imputed_count'] = oracle.adata.X
    oracle.adata.obs['cluster'] = 'cluster'
    oracle.cluster_column_name = None
    oracle.embedding_name = 'X_umap'
    oracle.pcs = np.zeros((oracle.adata.shape[0], 2))
    oracle.knn = True
    oracle.k_knn_imputation = True
    oracle.high_var_genes = list(oracle.adata.var_names)
    oracle.adata.obs['cluster'] = oracle.adata.obs['cluster'].astype('category')
    oracle.adata.uns['cluster_colors'] = ['#1f77b4']
    col_dict = co._get_clustercolor_from_anndata(adata=oracle.adata,
                                                cluster_name='cluster',
                                                return_as="dict")
    oracle.colorandum = np.array([col_dict[i] for i in oracle.adata.obs['cluster']])
    tf_dict = grn.groupby(['target'])['source'].apply(lambda x: sorted(list(x))).to_dict()   # Refit GRN
    oracle.addTFinfo_dictionary(tf_dict)
    # Add grn
    if fit_grn:
        oracle.fit_GRN_for_simulation(alpha=0, GRN_unit="whole")
    return oracle


def simulate_delta(oracle, tfs, n_steps=3):
    perturb_condition = dict()

    # Knockdown the specified transcription factors
    for g in tfs:
        perturb_condition[g] = 0  # Use zero for knockdown

    # Perform the simulation
    oracle.simulate_shift(perturb_condition=perturb_condition, n_propagation=n_steps)

    # Retrieve and process the delta values
    delta = pd.DataFrame(oracle.adata.layers['delta_X'], index=oracle.adata.obs_names, columns=oracle.adata.var_names)

    # Remove the columns that are 0s.
    delta = delta.loc[:, delta.abs().sum(0) != 0]

    return delta


def forecasting_benchmark(grn, rna, benchmark_matrix):
    grn = grn.copy()
    df_grn = grn[["TF", "Gene"]].copy()
    df_grn = df_grn.rename(columns={"TF": "source", "Gene": "target"})

    # Subset data to grn
    genes = set(df_grn['source']) | set(df_grn['target'])
    rna = rna[:, rna.var_names.isin(genes)]

    # Init object
    oracle = init_celloracle(rna, df_grn, fit_grn=True)
    coef_mat = oracle.coef_matrix
    oracle = ad.AnnData(oracle.adata.to_df().mean(0).to_frame().T)
    oracle = init_celloracle(oracle, df_grn, fit_grn=False)
    oracle.coef_matrix = coef_mat
    tf_n_trgs = (coef_mat != 0).sum(0)
    tf_n_trgs = set(tf_n_trgs[tf_n_trgs >= 3].index)
    oracle.all_regulatory_genes_in_TFdict = [t for t in oracle.all_regulatory_genes_in_TFdict if t in tf_n_trgs]

    # Read benchmark data
    obs = benchmark_matrix[["TF_name"]].copy()
    mat = benchmark_matrix.drop(columns=["TF_name"]).copy()

    # Subset bench data to dataset
    msk = obs['TF_name'].isin(rna.var_names)
    obs = obs.loc[msk, :]
    mat = mat.loc[msk, :]

    # Subset by overlap with rna
    genes = list(genes & set(mat.columns))
    mat = mat.loc[:, genes].copy()
    
    coefs = []
    pvals = []
    for dataset in tqdm(obs.index):
        # Extract
        tf = obs.loc[dataset, 'TF_name']
        tf_mat = mat.loc[[dataset], :]
        tf_mat = tf_mat[tf_mat != 0].dropna(axis=1)
        if tf in oracle.all_regulatory_genes_in_TFdict:
            # Run the simulation for the current TF
            x = simulate_delta(oracle, [tf], n_steps=3)
            
            # Intersect
            y = tf_mat
            inter = np.intersect1d(x.columns, y.columns)
            x, y = x.loc[:, inter].values[0], y.loc[:, inter].values[0]

            # Compute correlation
            if x.size >= 10:
                r, p = scipy.stats.spearmanr(x, y)
            else:
                r, p = 0., 1.
            coefs.append(r)
            pvals.append(p)
    
    # Compute recall
    coefs = np.array(coefs)
    pvals = np.array(pvals)
    padj = scipy.stats.false_discovery_control(pvals, method='bh')
    tp = np.sum((coefs > 0.05) & (padj < 0.05))
    fp = coefs.size - tp
    fn = obs.shape[0] - tp
    if tp > 0:
        prc = tp / coefs.size
        rcl = tp / obs.shape[0]
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0., 0., 0.
    
    return tp, fp, fn, prc, rcl, f01

def frc_test(
    grn_inferred: pd.DataFrame,
    adata,
    frc_matrix,
    score_column,
    step: float = 0.05,
    beta: float = 0.1
) -> dict:
    """
    Benchmark a GRN's ability to predict perturbation outcomes using a forecasting benchmark.

    This function applies threshold-based filtering on the GRN's interaction weights
    and evaluates each filtered network using a user-defined forecasting benchmark.
    The benchmarking function must return (precision, recall, fbeta).

    Parameters
    ----------
    grn : pd.DataFrame
        Inferred gene regulatory network with interaction scores.
    adata : AnnData
        Single-cell RNA-seq dataset used in the forecasting benchmark.
    frc_golden : str or pd.DataFrame
        Path or object representing the benchmark target outcomes.
    score_column : str, optional
        Column name containing GRN edge weights (default is "weight").
    step : float, optional
        Threshold step size (default is 0.05).
    beta : float, optional
        Beta value used to compute F-beta score (default is 0.1).

    Returns
    -------
    dict
        Dictionary of benchmarking metrics including precision, recall, F-beta,
        AUPRC, and best-threshold scores.
    """
    return threshold_benchmarking(
        grn=grn_inferred,
        benchmark_func=forecasting_benchmark,
        benchmark_kwargs={
            "grn": grn_inferred,
            "rna": adata,
            "benchmark_matrix": frc_matrix
        },
        score_column=score_column,
        step=step
    )


###########################################
#                 OMICS                   #
###########################################  


def test_predictability(mdata, train, test, grn, col_source='source', col_target='target', mod_source='rna', mod_target='rna', ntop=5):
    def remove_zeros(X, y):
        msk = y != 0.
        y = y[msk]
        X = X[msk, :]
        return X, y
    net = grn.iloc[np.argsort(-abs(grn['score'])), :].drop_duplicates([col_source, col_target])
    net = net.groupby(col_target)[col_source].apply(lambda x: list(x) if ntop is None else list(x)[:ntop])
    cor = []
    for target in tqdm(net.index):
        sources = net[target]
        sources = [s for s in sources if s != target]
        if len(sources) > 0:
            train_X = mdata.mod[mod_source][train, sources].X
            train_y = mdata.mod[mod_target][train, target].X
            train_y = train_y.toarray().ravel() if scipy.sparse.issparse(train_y) else np.ravel(train_y)

            test_X = mdata.mod[mod_source][test, sources].X
            test_y = mdata.mod[mod_target][test, target].X
            test_y = test_y.toarray().ravel() if scipy.sparse.issparse(test_y) else np.ravel(test_y)
            if test_y.size >= 10:
                reg = XGBRegressor(random_state=0).fit(train_X, train_y)
                pred_y = reg.predict(test_X)
                if np.any(pred_y != pred_y[0]):
                    s, p = scipy.stats.spearmanr(pred_y, test_y)  # Spearman to control for outliers
                    cor.append([target, pred_y.size, len(sources), s, p])
    cor = pd.DataFrame(cor, columns=['target', 'n_obs', 'n_vars', 'coef', 'pval'])
    if cor.shape[0] > 0:
        valid = (cor['pval'] >= 0) & (cor['pval'] <= 1) & np.isfinite(cor['pval'])
        cor = cor[valid]
        if cor.shape[0] > 0:
            cor['padj'] = scipy.stats.false_discovery_control(cor['pval'], method='bh')
        else:
            cor['padj'] = pd.Series(dtype=float)
    else:
        cor['padj'] = pd.Series(dtype=float)
    return cor



def omics_benchmark(
    grn: pd.DataFrame,
    mdata: mu.MuData,
    mod_source: str,
    mod_target: str,
    celltype_column: str,
) -> pd.DataFrame:
    """
    Evaluate the predictability of a gene regulatory network (GRN) using a test set and compute key performance metrics.

    Parameters
    ----------
    grn : pd.DataFrame
        The inferred gene regulatory network as a pandas DataFrame. If empty, performance metrics will be set as NaN.
    grn_name : str
        Name of the GRN being evaluated.
    mdata : muon.MuData
        Multi-omics data object containing cell annotations and modality-specific gene names.
    mod_target : str
        The target modality name used to extract variable gene names (e.g., 'rna').
    test_predictability : function
        A function that computes predictability correlations between regulator activity and target gene expression.
        It should accept the arguments `mdata`, `train`, `test`, and `grn`.
    f_beta_score : function
        A function that computes the F0.1 score given precision and recall.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns ['name', 'prc', 'rcl', 'f01'], containing the GRN name and its precision, recall,
        and F0.1 score, respectively.
    """
    train, test = train_test_split(
        mdata.obs_names,
        test_size=0.33,
        random_state=42,
        stratify=mdata.obs[celltype_column]
    )
    cor = test_predictability(mdata=mdata, train=train, test=test, grn=grn, mod_source=mod_source, mod_target=mod_target)
    sig_cor = cor[(cor['padj'] < 0.05) & (cor['coef'] > 0.05)]
    n_hits = sig_cor.shape[0]

    if n_hits > 0:
        universe_size = mdata.mod[mod_target].var_names.size
        rcl = n_hits / universe_size
        prc = n_hits / cor.shape[0]
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0.0, 0.0, 0.0

    tp = n_hits
    fp = universe_size - tp
    fn = cor.shape[0] - tp

    return tp, fp, fn, prc, rcl, f01


def omic_test(
    grn_inferred,
    mdata,
    source_column,
    target_column,
    score_column,
    mod_source,
    mod_target,
    celltype_column,
    step = 0.2
):
    """
    Evaluate the performance of a gene regulatory network (GRN) using score-based thresholding 
    and compute key performance metrics through benchmarking.

    This function takes an inferred GRN and:
    1. Renames the GRN columns to standardize the source, target, and score fields.
    2. Calls a threshold benchmarking function to systematically evaluate the GRN's performance 
       at different score thresholds.

    Parameters
    ----------
    grn_inferred : pd.DataFrame
        The inferred gene regulatory network. Must include columns for transcription factor (TF), 
        target gene (Gene), and the score column specified by `score_col`.
    mdata : muon.MuData
        Multi-omics data object containing the expression/accessibility matrices used for predictability testing.
    score_col : str
        The column name in `grn_inferred` containing the score that ranks the regulatory edges.
    mod_target : str
        The target modality (e.g., 'scRNA_counts') to use as the target gene expression data.
    celltype_column : str
        The column in `mdata.obs` that contains cell type annotations for stratified splitting.


    Returns
    -------
    pd.DataFrame
        A DataFrame containing benchmarking results across different score thresholds. Typically includes
        columns such as ['name', 'prc', 'rcl', 'f01'] per threshold, summarizing precision, recall, and F0.1 score.
    """
    grn = grn_inferred.copy()
    grn = grn.rename(columns={source_column: "source", target_column: "target", score_column: "score"})

    return threshold_benchmarking(
        grn=grn,
        benchmark_func=omics_benchmark,
        benchmark_kwargs={
            "grn": grn,
            "mdata": mdata,
            "mod_source": mod_source, 
            "mod_target": mod_target,
            "celltype_column": celltype_column
        },
        score_column="score",
        step=step
    )


###########################################
#               GENE SETS                 #
###########################################  

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

def get_sig_pws(grn: pd.DataFrame, db: pd.DataFrame, thr_pval: float = 0.01) -> np.ndarray:
    """
    For each TF in the GRN, perform one-sided Fisher's exact test
    to test overrepresentation of its targets in each gene set.
    Collect gene sets with BH-adjusted p-values < thr_pval.
    
    Parameters:
    - grn: pd.DataFrame with columns ['source', 'target']
    - db: pd.DataFrame with columns ['source', 'target'] representing gene sets
    - thr_pval: float, FDR threshold (default 0.01)
    
    Returns:
    - sig_pws: np.ndarray of significant gene sets
    """
    sig_pws = set()
    
    # Preprocess gene sets from db
    all_genes = set(grn['target'].unique())
    gene_sets = db.groupby('source')['target'].apply(set)
    
    for tf in grn['source'].unique():
        tf_targets = set(grn.loc[grn['source'] == tf, 'target'])
        if len(tf_targets) == 0:
            continue
        
        pvals = []
        terms = []
        
        for gene_set_name, gene_set in gene_sets.items():
            # Build contingency table:
            #            in gene set     not in gene set
            # in TF targets     A                  B
            # not in TF targets C                  D
            A = len(tf_targets & gene_set)
            B = len(tf_targets - gene_set)
            C = len(gene_set - tf_targets)
            D = len(all_genes - (tf_targets | gene_set))
            
            # Skip gene sets with no overlap
            if A == 0:
                continue
            
            table = np.array([[A, B], [C, D]])
            _, pval = fisher_exact(table, alternative='greater')
            pvals.append(pval)
            terms.append(gene_set_name)
        
        # Multiple testing correction
        if pvals:
            padj = multipletests(pvals, method='fdr_bh')[1]
            for term, adj_p in zip(terms, padj):
                if adj_p < thr_pval:
                    sig_pws.add(term)
    
    return np.array(list(sig_pws))


def eval_metrics(y_pred, y):
    tp = np.intersect1d(y_pred, y).size
    if tp > 0.:
        fp = np.setdiff1d(y_pred, y).size
        fn = np.setdiff1d(y, y_pred).size
        prc = tp / (tp + fp)
        rcl = tp / (tp + fn)
        f1 = f_beta_score(prc, rcl)
    else:
        tp, fp, fn, prc, rcl, f1 = 0., 0., 0.
    return tp, fp, fn, prc, rcl, f1


def eval_grn(data, grn, db, thr_pval=0.01, thr_prop=0.2):
    hits = data.uns['ulm_hits']
    sig_pws = get_sig_pws(grn, db, thr_pval)
    tp, fp, fn, prc, rcl, f1 = eval_metrics(y_pred=sig_pws, y=hits)
    return tp, fp, fn, prc, rcl, f1


def get_pw_hits(data, thr_pval, thr_prop):
    pvals = data.obsm['padj_ulm'].copy()  # decoupler already applied BH
    acts = data.obsm['score_ulm'].copy()

    hits = ((pvals < thr_pval) & (acts > 0)).sum(0).sort_values(ascending=False) / pvals.shape[0]
    hits = hits[hits > thr_prop].index.values.astype('U')
    return hits


def gst_benchmark(
    ptw : pd.DataFrame,
    rna : mu.MuData,
    grn : pd.DataFrame,
    thr_pval: float = 0.01,
    thr_prop: float = 0.2
):

    # Infer pathway activities
    if 'ulm_hits' not in rna.uns:
        dc.mt.ulm(
            data=rna,
            net=ptw,
            verbose=True
        )

        hits = get_pw_hits(rna, thr_pval, thr_prop)
        rna.uns['ulm_hits'] = hits

    return eval_grn(rna, grn, ptw, thr_pval=0.01, thr_prop=0.2)

def gst_test(
    grn_inferred: pd.DataFrame,
    ptw: pd.DataFrame,
    rna: mu.MuData,
    score_column,
    thr_pval: float = 0.01,
    thr_prop: float = 0.2,
    step: float = 0.01
):

    grn = grn_inferred.copy()
    grn = grn.rename(columns={"TF": "source", "Gene": "target"})
    
    return threshold_benchmarking(
        grn=grn,
        benchmark_func=gst_benchmark,
        benchmark_kwargs={
            "ptw": ptw,
            "rna": rna,
            "grn": grn,
            "thr_pval": thr_pval,
            "thr_prop": thr_prop
        },
        score_column=score_column,
        step=step
    )


###########################################
#              TF MARKERS                 #
########################################### 


def tf_marker_benchmark(grn, db, adata):
    """
    Evaluate a GRN against a database of known TF markers.

    Parameters
    ----------
    grn : pd.DataFrame
        DataFrame with at least a 'source' column listing the predicted regulators.
    db : pd.DataFrame
        DataFrame with a single column 'gene' containing known TF markers.
    adata : AnnData or MuData
        Object containing the measured gene names in .var_names.

    Returns
    -------
    tp : int
        Number of true positives.
    fp : int
        Number of false positives.
    fn : int
        Number of false negatives.
    prc : float
        Precision.
    rcl : float
        Recall.
    f01 : float
        F0.1 score.
    """
    
    # Filter the resource database by measured genes
    db_filtered = db[db['gene'].astype(str).isin(adata.var_names.astype(str))]
    
    # Compute evaluation metrics
    y_pred = grn['source'].unique().astype(str)
    y_true = db_filtered['gene'].unique().astype(str)
    
    tp = np.intersect1d(y_pred, y_true).size
    fp = np.setdiff1d(y_pred, y_true).size
    fn = np.setdiff1d(y_true, y_pred).size
    
    if tp > 0:
        prc = tp / (tp + fp)
        rcl = tp / (tp + fn)
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0.0, 0.0, 0.0
    
    return tp, fp, fn, prc, rcl, f01


def tfm_test(
    grn_inferred: pd.DataFrame,
    db: pd.DataFrame,
    adata: mu.MuData,
    score_column,
    step: float = 0.01,
):
    """
    Benchmark a GRN's ability to recover known TF markers.

    This function applies threshold-based filtering on the GRN's interaction weights
    and evaluates each filtered network using a user-defined TF marker benchmark.
    The benchmarking function must return (precision, recall, F-beta).

    Parameters
    ----------
    grn : pd.DataFrame
        Inferred gene regulatory network with interaction scores.
    db : pd.DataFrame
        Database of known TF markers.
    adata : AnnData or MuData
        Object containing the measured gene names in .var_names.
    score_column : str, optional
        Column name containing GRN edge weights (default is "weight").
    step : float, optional
        Threshold step size (default is 0.01).
    beta : float, optional
        Beta value used to compute F-beta score (default is 0.1).

    Returns
    -------
    dict
        Dictionary of benchmarking metrics including precision, recall, F-beta,
        AUPRC, and best-threshold scores.
    """

    grn = grn_inferred.copy()
    grn = grn.rename(columns={"TF": "source"})
    
    return threshold_benchmarking(
        grn=grn,
        benchmark_func=tf_marker_benchmark,
        benchmark_kwargs={
            "grn": grn,
            "db": db,
            "adata": adata
        },
        score_column=score_column,
        step=step
    )