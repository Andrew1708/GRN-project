#!/usr/bin/env python3

import argparse
import pandas as pd
from celloracle import motif_analysis as ma
import celloracle as co

def parse_args():
    parser = argparse.ArgumentParser(
        description="Integrate scATAC-seq peaks with Cicero coaccessibility scores "
                    "and annotate with TSS information."
    )
    parser.add_argument(
        "--peaks",
        type=str,
        required=True,
        help="Path to the all_peaks.csv file."
    )
    parser.add_argument(
        "--cicero",
        type=str,
        required=True,
        help="Path to the cicero_connections.csv file."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="processed_peak_file.csv",
        help="Path to the output processed peaks file."
    )
    parser.add_argument(
        "--ref-genome",
        type=str,
        default="hg38",
        help="Reference genome to use (e.g., 'hg38' or 'mm10')."
    )
    return parser.parse_args()

def main():
    args = parse_args()

    # Load scATAC-seq peak list
    print("Reading peaks file...")
    peaks = pd.read_csv(args.peaks)
    peaks = peaks.x.values

    # Load Cicero coaccessibility scores
    print("Reading Cicero connections file...")
    cicero_connections = pd.read_csv(args.cicero, index_col=0)

    # Annotate TSS information
    print(f"Annotating peaks with TSS for reference genome '{args.ref_genome}'...")
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=args.ref_genome)

    # Integrate TSS with Cicero connections
    print("Integrating TSS and Cicero connections...")
    integrated = ma.integrate_tss_peak_with_cicero(
        tss_peak=tss_annotated,
        cicero_connections=cicero_connections
    )

    # Filter high coaccessibility peaks
    print("Filtering high coaccessibility peaks (threshold >= 0.8)...")
    filtered_peaks = integrated[integrated.coaccess >= 0.8]
    filtered_peaks = filtered_peaks[["peak_id", "gene_short_name"]].reset_index(drop=True)

    # Save the processed peak file
    print(f"Saving processed peak file to '{args.output}'...")
    filtered_peaks.to_csv(args.output, index=False)

    print("Processing completed successfully.")

if __name__ == "__main__":
    main()
