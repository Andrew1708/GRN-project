# summarise_grn_streaming.py
import os
import csv
import gc
from pathlib import Path
from typing import Dict, Iterable

import pandas as pd
import mudata


TF_COL   = "TF"        # regulator column in every GRN table
GENE_COL = "Gene"      # target-gene column in every GRN table
OUT_COLS = ("sample", "n_triplets", "n_tfs", "n_genes")


# ------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------
def _summary_record(sample: str, grn: pd.DataFrame) -> Dict[str, object]:
    """Return basic GRN statistics for one sample."""
    if TF_COL not in grn.columns or GENE_COL not in grn.columns:
        missing = {TF_COL, GENE_COL}.difference(grn.columns)
        raise KeyError(f"Missing expected column(s): {missing}")

    return {
        "sample":    sample,
        "n_triplets": int(grn.shape[0]),
        "n_tfs":     int(grn[TF_COL].nunique()),
        "n_genes":   int(grn[GENE_COL].nunique()),
    }


def _completed_samples(csv_path: Path) -> set[str]:
    """Return the set of sample identifiers already stored in *csv_path*."""
    if not csv_path.exists():
        return set()
    return set(pd.read_csv(csv_path, usecols=["sample"])["sample"])


def _append_rows(csv_path: Path, rows: Iterable[Dict[str, object]]) -> None:
    """Append *rows* to *csv_path*, writing a header only when the file is new."""
    write_header = not csv_path.exists()
    with csv_path.open("a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_COLS)
        if write_header:
            writer.writeheader()
        for r in rows:
            writer.writerow(r)


# ------------------------------------------------------------------
# Per-tool processing
# ------------------------------------------------------------------
def _process_tool(
    tool_name: str,
    root_dir: Path,
    h5mu_name: str,
    extract_grn,
    out_path: Path,
) -> None:
    """
    Iterate over *root_dir/*sample/h5mu_name*, compute summary statistics,
    and append them to *out_path* for samples that are not yet recorded.

    *extract_grn* is a callable that receives *mudata.MuData* and returns
    a ``pd.DataFrame`` with columns ``TF_COL`` and ``GENE_COL``.
    """
    done = _completed_samples(out_path)
    new_rows = []

    for entry in root_dir.iterdir():
        if not entry.is_dir():
            continue
        sample = entry.name
        if sample in done:
            continue

        f = entry / h5mu_name
        if not f.exists():
            continue

        try:
            mdata = mudata.read(f)
            grn = extract_grn(mdata)
            new_rows.append(_summary_record(sample, grn))
        except Exception as e:
            print(f"[WARN] {tool_name} ⇒ {sample}: {e}")
        finally:
            # free memory before the next loop iteration
            del mdata, grn
            gc.collect()

        # optionally flush every N samples to keep memory footprint tiny
        if len(new_rows) >= 20:
            _append_rows(out_path, new_rows)
            new_rows.clear()

    if new_rows:                              # leftovers
        _append_rows(out_path, new_rows)


# ------------------------------------------------------------------
# Public entry point
# ------------------------------------------------------------------
def get_and_save_n_triplets(
    scenicplus_dir: Path,
    celloracle_dir: Path,
    linger_dir: Path,
    out_dir: Path,
    scenic_name: str = "scplusmdata.h5mu",
    celloracle_name: str = "mdata.h5mu",
    linger_name: str = "mdata.h5mu",
) -> None:
    """Stream GRN statistics from three tools into per-tool CSV files."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- SCENICPLUS ------------------------------------------------
    # try:
    #     _process_tool(
    #         tool_name="SCENICPLUS",
    #         root_dir=Path(scenicplus_dir),
    #         h5mu_name=scenic_name,
    #         extract_grn=lambda m: pd.concat(
    #             [pd.DataFrame(m.uns["direct_e_regulon_metadata"]),
    #              pd.DataFrame(m.uns["extended_e_regulon_metadata"])],
    #             ignore_index=True,
    #         ),
    #         out_path=out_dir / "scenicplus.csv",
    #     )
    # except Exception as e:
    #     print(f"[FATAL] Failed to process SCENICPLUS: {e}")

    # ---- CellOracle ------------------------------------------------
    try:
        _process_tool(
            tool_name="CellOracle",
            root_dir=Path(celloracle_dir),
            h5mu_name=celloracle_name,
            extract_grn=lambda m: pd.DataFrame(m.uns["celloracle_links"])[["source", "target"]].rename(
                columns={"source": "TF", "target": "Gene"}
            ),
            out_path=out_dir / "celloracle.csv",
        )
    except Exception as e:
        print(f"[FATAL] Failed to process CellOracle: {e}")

    # ---- LINGER ----------------------------------------------------
    # try:
    #     _process_tool(
    #         tool_name="LINGER",
    #         root_dir=Path(linger_dir),
    #         h5mu_name=linger_name,
    #         extract_grn=lambda m: pd.DataFrame(m.uns["grn"]),
    #         out_path=out_dir / "linger.csv",
    #     )
    # except Exception as e:
    #     print(f"[FATAL] Failed to process LINGER: {e}")


# ------------------------------------------------------------------
# Command-line use
# ------------------------------------------------------------------
if __name__ == "__main__":
    SCENICPLUS_DIR = Path("../../data/output/scenicplus/results/")
    CELLORACLE_DIR = Path("../../data/output/celloracle/out")
    LINGER_DIR     = Path("../../data/output/linger/out")
    OUT_DIR        = Path("../results")

    get_and_save_n_triplets(
        scenicplus_dir=SCENICPLUS_DIR,
        celloracle_dir=CELLORACLE_DIR,
        linger_dir=LINGER_DIR,
        out_dir=OUT_DIR,
    )
