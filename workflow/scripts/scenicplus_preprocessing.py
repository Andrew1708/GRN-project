import argparse
import pickle
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates
import os
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np
import scanpy as sc


cell_type_col_name = "Classified_Celltype"
malignant_col_name = "Malignant"
metacell_key = "Metacell_Key"

def topic_binarization(cistopic_object, out_dir):
    region_bin_topics_top_3k = binarize_topics(
        cistopic_object, method='ntop', ntop = 3_000,
        plot=True, num_columns=5
    )

    region_bin_topics_otsu = binarize_topics(
        cistopic_object, method='otsu',
        plot=True, num_columns=5
    )

    binarized_cell_topic = binarize_topics(
        cistopic_object,
        target='cell',
        method='li',
        plot=True,
        num_columns=5, nbins=100
    )

    for topic in region_bin_topics_otsu:
        region_names_to_coordinates(
            region_bin_topics_otsu[topic].index
        ).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
            sep = "\t",
            header = False, index = False
        )

    for topic in region_bin_topics_top_3k:
        region_names_to_coordinates(
            region_bin_topics_top_3k[topic].index
        ).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
            sep = "\t",
            header = False, index = False
        )


def compute_DAR(cistopic_object, out_dir, temp_dir):
    imputed_acc_obj = impute_accessibility(
        cistopic_object,
        selected_cells=None,
        selected_regions=None,
        scale_factor=10**6
    )

    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj,
        min_disp = 0.05,
        min_mean = 0.0125,
        max_mean = 3,
        max_disp = np.inf,
        n_bins=20,
        n_top_features=None,
        plot=True
    )   

    markers_dict= find_diff_features(
        cistopic_object,
        imputed_acc_obj,
        variable= metacell_key,
        var_features=variable_regions,
        contrasts=None,
        adjpval_thr=0.05,
        log2fc_thr=np.log2(1.5),
        n_cpu=5,
        split_pattern = '___'
    )
    
    for cell_type in markers_dict:
        region_names_to_coordinates(
            markers_dict[cell_type].index
        ).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
            sep = "\t",
            header = False, index = False
        )
    
        # Define the directory containing .bed files
    directory = os.path.join(out_dir, "region_sets/DARs_cell_type")
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".bed"):  # Check if the file has a .bed extension
            file_path = os.path.join(directory, filename)

            # Check if the file is empty
            if os.path.isfile(file_path) and os.stat(file_path).st_size == 0:
                os.remove(file_path)  # Delete the empty file
                print(f"Deleted empty file: {filename}")

def normalize(adata):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata

def select_hvgs(adata, n_top_genes=None):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var.highly_variable].copy()
    return adata

def preprocess(adata, n_top_genes=None):
    adata.raw = adata
    adata = normalize(adata)
    adata = select_hvgs(adata, n_top_genes=n_top_genes)
    sc.pp.scale(adata)
    return adata

def treat_adata(adata, out_dir):
    preprocess(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    adata.write(os.path.join(out_dir, "rna_adata.h5ad"))


def parse_arguments():
    parser = argparse.ArgumentParser(description="Preprocess Cistopic object for SCENIC.")
    parser.add_argument(
        "--cistopic_object",
        type=str,
        required=True,
        help="Path to the Cistopic object."
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Output directory for results."
    )
    parser.add_argument(
        "--temp_dir",
        type=str,
        required=True,
        help="Temporary directory for intermediate files."
    )

    parser.add_argument(
        "--rna_adata",
        type=str,
        required=True,
        help="Path to the RNA AnnData object."
    )

    parser.add_argument(
        "--sample_name",
        type=str,
        required=True,
        help="Sample name for the output files."
    )

    return parser.parse_args()

def create_metacell_keys(cistopic_object, adata_rna):
    # Create metacell key from ATAC
    cistopic_object.cell_data[metacell_key] = (
        cistopic_object.cell_data[cell_type_col_name].astype(str) + "." +
        cistopic_object.cell_data[malignant_col_name].astype(str)
    )

    # Create metacell key from RNA
    adata_rna.obs[metacell_key] = (
        adata_rna.obs[cell_type_col_name].astype(str) + "." +
        adata_rna.obs[malignant_col_name].astype(str)
    )


def main():
    args = parse_arguments()

    cistopic_object = pickle.load(open(args.cistopic_object, "rb"))
    # Load the RNA AnnData object
    rna_adata = sc.read_h5ad(args.rna_adata)

    # Create metacell keys
    create_metacell_keys(cistopic_object, rna_adata)


    out_dir = args.out_dir

    # Create output directories if they don't exist
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok=True)

    cistopic_object.region_data.to_csv(f"{out_dir}/regions.bed", sep='\t', header=False, index=False)

    # Perform topic binarization
    topic_binarization(cistopic_object, out_dir)
    print("Topic binarization completed.")

    # Compute differential accessibility regions (DARs)
    compute_DAR(cistopic_object, out_dir, args.temp_dir)
    print("Differential accessibility regions computed.")

    # Preprocess the RNA data
    treat_adata(rna_adata, args.out_dir)
    print("RNA data preprocessing completed.")

    # Save pycistopic obj
    cistopic_out_dir = os.path.join(out_dir, f"cistopic_object.pkl")
    with open(cistopic_out_dir, "wb") as f:
        pickle.dump(cistopic_object, f)

if __name__ == "__main__":
    main()