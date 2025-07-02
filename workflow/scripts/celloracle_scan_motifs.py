#!/usr/bin/env python3

import argparse
import pandas as pd
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
import genomepy


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scan scATAC-seq peaks for transcription factor motifs using CellOracle."
    )
    parser.add_argument(
        "--peaks",
        type=str,
        required=True,
        help="Path to the input peak file (CSV)."
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="output",
        help="Prefix for output files (TFinfo object and DataFrame)."
    )
    parser.add_argument(
        "--ref-genome",
        type=str,
        default="mm10",
        help="Reference genome to use (e.g., 'mm10' or 'hg38')."
    )
    parser.add_argument(
        "--fpr",
        type=float,
        default=0.02,
        help="False positive rate for motif scanning (default: 0.02)."
    )
    parser.add_argument(
        "--score-threshold",
        type=float,
        default=10.0,
        help="Threshold for filtering motif scores (default: 10)."
    )
    return parser.parse_args()


def ensure_genome_installed(ref_genome, genomes_dir=None):
    if not ma.is_genome_installed(ref_genome=ref_genome, genomes_dir=genomes_dir):
        print(f"Installing reference genome '{ref_genome}'...")
        genomepy.install_genome(name=ref_genome, provider="UCSC", genomes_dir=genomes_dir)
    else:
        print(f"Reference genome '{ref_genome}' is already installed.")


def run_motif_scan(peaks_path, ref_genome, fpr, score_threshold, output_prefix):
    print("Loading peaks...")
    peaks = pd.read_csv(peaks_path)
    peaks = ma.check_peak_format(peaks, ref_genome)

    print("Creating TFinfo object...")
    tfi = ma.TFinfo(peak_data_frame=peaks, ref_genome=ref_genome)

    print("Scanning motifs...")
    tfi.scan(fpr=fpr, motifs=None, verbose=True)

    tfinfo_path = f"{output_prefix}/info.celloracle.tfinfo"
    print(f"Saving TFinfo object to '{tfinfo_path}'...")
    tfi.to_hdf5(file_path=tfinfo_path)

    print(f"Filtering motifs with score >= {score_threshold}...")
    tfi.reset_filtering()
    tfi.filter_motifs_by_score(threshold=score_threshold)

    print("Formatting TFinfo DataFrame...")
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    df = tfi.to_dataframe()

    df_path = f"{output_prefix}/base_GRN_dataframe.parquet"
    print(f"Saving TFinfo DataFrame to '{df_path}'...")
    df.to_parquet(df_path)


def main():
    args = parse_args()
    ensure_genome_installed(args.ref_genome)
    run_motif_scan(
        peaks_path=args.peaks,
        ref_genome=args.ref_genome,
        fpr=args.fpr,
        score_threshold=args.score_threshold,
        output_prefix=args.output_prefix
    )
    print("Motif scanning completed.")


if __name__ == "__main__":
    main()
