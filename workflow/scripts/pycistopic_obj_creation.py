import pandas as pd
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl
import os
import pickle
import os
import argparse

N_CPU = 10

def load_tsv_as_dataframe(file_path):
    """
    Load a .tsv file into a pandas DataFrame with standardized column names.
    
    Parameters:
    file_path (str): Path to the .tsv file.

    Returns:
    pd.DataFrame: DataFrame with renamed columns:
        - "Chromosome"
        - "Start"
        - "End"
        - "Name"
        - "Strand"
    """
    # Define expected column names (assuming order in the file is correct)
    column_names = ["Chromosome", "Start", "End", "Name", "Strand"]
    
    # Load the file
    df = pd.read_csv(file_path, sep="\t", header=None, names=column_names)

    return df


def create_cistopic_objects(fragment_path, region_path, barcode_path, path_to_blacklist, project, temp_dir):
    """
    Create cistopic objects from fragments files and save them as pickle files.

    Parameters:
    fragments_files (list): List of fragment files.
    regions_dict (dict): Dictionary mapping fragment file names to region file paths.
    barcode_dict (dict): Dictionary mapping fragment file names to barcode file paths.
    path_to_blacklist (str): Path to the blacklist file.
    temp_dir (str): Directory to save the temporary files.
    source_dir (str): Source directory containing the fragment files.
    """

    df = load_tsv_as_dataframe(fragment_path)

    cistopic_object = create_cistopic_object_from_fragments(
        path_to_fragments = fragment_path,
        path_to_regions= region_path,
        path_to_blacklist = path_to_blacklist,
        valid_bc = barcode_path,
        fragments_df = df,
        project = project, #VERY SPECIFIC TO THIS DATASET
        n_cpu = N_CPU,
    )

    return cistopic_object

def add_metadata_to_cistopic_object_list(cistopic_obj, metadata_path, out_dir):
    # Load metadata
    df = pd.read_csv(metadata_path, sep="\t")
    df["Modified_Index"] = df["Cell_Barcode"] + "___" + df["Library_ID"]
    df.set_index("Modified_Index", inplace=True)

    # Add metadata
    cistopic_obj.add_cell_data(df)

    # Save updated object
    file_name = f"{cistopic_obj.project}_pycistopic_obj.pkl"
    with open(os.path.join(out_dir, file_name), "wb") as f:
        pickle.dump(cistopic_obj, f)

    print(f"Saved: {file_name}")



def parse_args():
    parser = argparse.ArgumentParser(description="Process directories for barcode and region files.")
    parser.add_argument("--temp_dir", required=True, help="Temporary working directory")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--path_to_fragments", required=True, help="Path to the fragments file")
    parser.add_argument("--path_to_regions", required=True, help="Path to the regions file")
    parser.add_argument("--path_to_barcodes", required=True, help="Path to the barcodes file")
    parser.add_argument("--path_to_blacklist", required=True, help="Path to the blacklist file")
    parser.add_argument("--project", required=True, help="Project name for the cistopic object")
    parser.add_argument("--metadata_path", required=True, help="Path to the metadata file")
    return parser.parse_args()

def main():
    args = parse_args()

    temp_dir = args.temp_dir
    out_dir = args.out_dir
    path_to_fragments = args.path_to_fragments
    path_to_regions = args.path_to_regions
    path_to_barcodes = args.path_to_barcodes
    path_to_blacklist = args.path_to_blacklist
    project = args.project
    metadata_path = args.metadata_path

    # (Optional) Create directories if needed
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    # Create cistopic objects
    cistopic_obj = create_cistopic_objects(
        fragment_path=path_to_fragments,
        region_path=path_to_regions, 
        barcode_path=path_to_barcodes, 
        path_to_blacklist=path_to_blacklist, 
        project=project, 
        temp_dir=temp_dir
    )

    # Add metadata to cistopic objects
    add_metadata_to_cistopic_object_list(
        cistopic_obj=cistopic_obj, 
        metadata_path=metadata_path,
        out_dir=out_dir
    )

if __name__ == "__main__":
    main()

