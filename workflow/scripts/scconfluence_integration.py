import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import mudata as md
import scanpy as sc
import scconfluence
import torch
from scipy.spatial.distance import cdist

# Local imports
sys.path.append(str(Path(__file__).resolve().parent.parent))
from utils.mutual_nearest_neighbors import *

def main(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load MuData
    print("Reading input MuData file...")
    mdata = md.read_h5mu(args.input)

    # --- PATCH: Make barcodes unique by adding modality suffix ---
    print("Applying modality suffix to barcodes...")
    original_rna_names = mdata["rna"].obs_names.copy()
    original_atac_names = mdata["atac"].obs_names.copy()

    mdata["rna"].obs_names = [f"{name}#RNA" for name in original_rna_names]
    mdata["atac"].obs_names = [f"{name}#ATAC" for name in original_atac_names]

    # Compute cross-modality distance
    print("Computing cross-modality distances...")
    key = "cross_rna+atac"
    mdata.uns[key] = cdist(mdata["rna"].obsm["cm_features"],
                           mdata["atac"].obsm["cm_features"])
    mdata.uns["cross_keys"] = [key]

    # Initialize autoencoders
    print("Initializing autoencoders...")
    autoencoders = {
        "rna": scconfluence.model.AutoEncoder(mdata["rna"], modality="rna"),
        "atac": scconfluence.model.AutoEncoder(mdata["atac"], modality="atac")
    }

    # Initialize and optionally move model to GPU
    print("Initializing ScConfluence model...")
    model = scconfluence.model.ScConfluence(mdata, unimodal_aes=autoencoders)
    if hasattr(model, "to"):
        model = model.to(device)

    # Train model
    print("Training model...")
    model.fit(save_path=args.output)

    # Extract latent embeddings
    print("Extracting latent embeddings...")
    z_all = model.get_latent()

    # Extract modality-specific embeddings and names (safe)
    rna_obs_names = mdata["rna"].obs_names
    atac_obs_names = mdata["atac"].obs_names

    z_rna = z_all.loc[rna_obs_names].to_numpy()
    z_atac = z_all.loc[atac_obs_names].to_numpy()

    # Strip suffixes to restore original names in same order
    rna_names = [name.replace("#RNA", "") for name in rna_obs_names]
    atac_names = [name.replace("#ATAC", "") for name in atac_obs_names]

    # Final guarantee
    assert len(rna_names) == z_rna.shape[0], "Mismatch in RNA name and embedding count"
    assert len(atac_names) == z_atac.shape[0], "Mismatch in ATAC name and embedding count"

    # Save embeddings
    print(f"Saving latent embeddings to: {args.output}")
    os.makedirs(args.output, exist_ok=True)
    np.save(os.path.join(args.output, "z_rna.npy"), z_rna)
    np.save(os.path.join(args.output, "z_atac.npy"), z_atac)

    # Perform mutual nearest neighbor matching
    print("Matching cells using mutual nearest neighbors...")
    matches = mutual_nearest_neighbors(z_rna, z_atac, rna_names, atac_names)

    matches.to_csv(os.path.join(args.output, "matches.csv"), index=False)
    print(f"Saved matches for {args.project_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ScConfluence on RNA + ATAC MuData")
    parser.add_argument("--input", type=str, required=True, help="Path to input .h5mu file")
    parser.add_argument("--output", type=str, default="results", help="Directory to save model outputs")
    parser.add_argument("--project_name", type=str, default="scconfluence_run", help="Project name for logs")
    args = parser.parse_args()

    main(args)
