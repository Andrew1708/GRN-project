import argparse
import scanpy as sc
import pandas as pd
import pickle
import os
import scipy.io
import scipy.sparse
import numpy as np
import ast


def argparser():
    parser = argparse.ArgumentParser(description="Create a multiome dataset from single-cell RNA and ATAC data.")
    parser.add_argument("--match_file", type=str, required=True, help="Path to the file containing matched cell IDs for RNA and ATAC data.")
    parser.add_argument("--rna_file", type=str, required=True, help="Path to the RNA data file.")
    parser.add_argument("--cistopic_file", type=str, required=True, help="Path to the ATAC data file.")
    parser.add_argument("--input_rna_dir", type=str, required=True, help="Directory containing the RNA data files.")
    parser.add_argument("--input_atac_dir", type=str, required=True, help="Directory containing the ATAC data files.")
    parser.add_argument("--cistopic_dir", type=str, required=True, help="Directory containing the Cistopic data files.")
    parser.add_argument("--project_name", type=str, required=True, help="Name of the project for saving the multiome dataset.")
    return parser.parse_args()

def main():
    args = argparser()

    # Load the matched cell IDs
    match_data = pd.read_csv(args.match_file, sep=",")

    # Filter only matches with a weight greater than 0.8 (to mitigate distant neighbors)
    match_data = match_data[match_data['Weight'] > 0.8]

    # Treat RNA data
    rna_adata = sc.read_h5ad(args.rna_file)

    # Filter RNA data to only include matched cells
    rna_adata = rna_adata[rna_adata.obs_names.isin(match_data['RNA']), :]

    # Save embeddings
    # Parse embedding strings (no commas)
    embedding_list = match_data['RNA_embedding'].apply(
        lambda x: np.fromstring(x.strip('[]'), sep=' ')
    )

    # Stack into 2D array
    embedding_array = np.stack(embedding_list)

    # Align rows if needed (ensure match_data['RNA'] is aligned with rna_adata.obs_names)
    embedding_df = pd.DataFrame(embedding_array, index=match_data['RNA'])
    embedding_df = embedding_df.loc[rna_adata.obs_names]

    # Store in .obsm
    rna_adata.obsm['integration_embedding'] = embedding_df.values

    # Rename cell barcodes to match ATAC data
    rna_to_atac = dict(zip(match_data['RNA'], match_data['ATAC']))
    rna_adata.obs_names = rna_adata.obs_names.to_series().map(rna_to_atac)

    # Load ATAC cistopic data
    atac_cistopic = pickle.load(open(args.cistopic_file, 'rb'))
    atac_cistopic = atac_cistopic.subset(match_data['ATAC'].values, copy=True)

    # Save embeddings
    embedding_df = pd.DataFrame(match_data['ATAC_embedding'], index=match_data['ATAC'])
    embedding_df = embedding_df.loc[atac_cistopic.cell_data.index]
    atac_cistopic.cell_data['integration_embedding'] = embedding_df.values

    # Save old cell type annotation
    atac_cistopic.cell_data['Old_Classified_Celltype'] = atac_cistopic.cell_data['Classified_Celltype'].copy()

    # Transfer classified cell type from RNA to ATAC
    atac_cistopic.cell_data['Classified_Celltype'] = (
        rna_adata.obs['Classified_Celltype']
        .reindex(atac_cistopic.cell_data.index)
    )

    atac_cistopic.project_name = args.project_name

    # Save the multiome dataset
    # RNA data
    rna_adata.write_h5ad(os.path.join(args.input_rna_dir, f"{args.project_name}_filtered.h5ad"))
    print(f"RNA data saved to {os.path.join(args.input_rna_dir, f'{args.project_name}_filtered.h5ad')}")
    # Cistopic 
    pickle.dump(
        atac_cistopic,
        open(os.path.join(args.cistopic_dir, f"{args.project_name}_pycistopic_obj.pkl"), 'wb')
    )
    print(f"ATAC data saved to {os.path.join(args.cistopic_dir, f'{args.project_name}_pycistopic_obj.pkl')}")
    # ATAC data: fragment matrix, .bed file and barcodes
    # 1) Save cell names as a TSV file
    cell_names = atac_cistopic.cell_names
    with open(os.path.join(args.input_atac_dir, f"{args.project_name}_filtered_barcodes.tsv"), 'w') as f:
        for name in cell_names:
            f.write(f"{name}\n")

    print(f"Cell names saved to {os.path.join(args.input_atac_dir, f'{args.project_name}_filtered_barcodes.tsv')}")
    # 2) Save fragment_matrix as an .mtx file
    # Assuming fragment_matrix is a scipy sparse matrix
    fragment_matrix = atac_cistopic.fragment_matrix
    scipy.io.mmwrite(os.path.join(args.input_atac_dir, f"{args.project_name}_filtered_peak_matrix.mtx"), fragment_matrix.astype(np.int32))
    print(f"Fragment matrix saved to {os.path.join(args.input_atac_dir, f'{args.project_name}_filtered_peak_matrix.mtx')}")

    # 3) Save region_data as a BED file
    # Assuming region_data is a pandas DataFrame
    region_data = atac_cistopic.region_data
    bed_data = region_data[['Chromosome', 'Start', 'End']]
    bed_data.to_csv(os.path.join(args.input_atac_dir, f"{args.project_name}_peaks.bed"), sep='\t', header=False, index=False)
    print(f"Region data saved to {os.path.join(args.input_atac_dir, f'{args.project_name}_peaks.bed')}")

if __name__ == "__main__":
    main()



