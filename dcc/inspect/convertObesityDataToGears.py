

#!/usr/bin/env python3
"""
convert_obesity_h5ad_to_gears.py

Read an obesity AnnData .h5ad and add the GEARS-required obs fields, then save
a new GEARS-ready .h5ad.

GEARS expects (at minimum) these adata.obs columns:
  - condition   : "GENE+ctrl" (single pert) or "GENE1+GENE2" (double pert) or "ctrl"
  - control     : 1 for ctrl else 0
  - cell_type   : a string label (can be constant if you have one cell type)
  - dose_val    : "1+1" for gene+ctrl, "1" for ctrl  (matches GEARS examples)
  - condition_name (optional but very commonly used): "<cell_type>_<condition>_<dose_val>"

This script assumes your obs has:
  - gene column: "gene" (perturbed gene symbol per cell)
  - guide UMI count column: "nCount_guide" (0 => ctrl by default)
  - optional one-hot cell-type columns: adipo, pre_adipo, other, lipo

Usage:
  python convert_obesity_h5ad_to_gears.py \
    --input  /path/to/obesity.h5ad \
    --output /path/to/obesity_gears.h5ad

Optional:
  --gene-col gene
  --guide-count-col nCount_guide
  --guide-threshold 0
  --cell-type-cols adipo,pre_adipo,other,lipo
  --default-cell-type adipocyte
  --ctrl-label ctrl
"""
# imports
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import scanpy as sc
except ImportError:
    print("ERROR: scanpy is required. Install with: pip install scanpy anndata", file=sys.stderr)
    raise


# constants
DIR_DATA = "/home/javaprog/Data/Broad/GeneticsML/Gears202512/Data/broad_obesity"
FILE_OBESITY_ORIGINAL = "{}/obesity_challenge_1.h5ad".format(DIR_DATA)
FILE_OBESITY_GEARS = "{}/perturb_processed.h5ad".format(DIR_DATA)
COUNT_CUTOFF_GUIDE_COUNTS = 50


# methods
def derive_cell_type_from_onehot(obs, onehot_cols: list[str]) -> np.ndarray:
    """
    Derive a categorical cell_type from one-hot columns by taking argmax per row.
    If a row is all zeros, returns "unknown".
    """
    mat = obs[onehot_cols].to_numpy()
    sums = mat.sum(axis=1)
    idx = mat.argmax(axis=1)
    ct = np.array([onehot_cols[i] for i in idx], dtype=object)
    ct[sums == 0] = "unknown"
    return ct


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", "-i", default=FILE_OBESITY_ORIGINAL, help="Input .h5ad")
    ap.add_argument("--output", "-o", default=FILE_OBESITY_GEARS, help="Output .h5ad (GEARS-ready)")
    ap.add_argument("--gene-col", default="gene", help="obs column with perturbed gene symbol (default: gene)")
    ap.add_argument(
        "--guide-count-col",
        default="nCount_guide",
        help="obs column with guide UMI counts used to identify controls (default: nCount_guide)",
    )
    ap.add_argument(
        "--guide-threshold",
        type=float,
        default=COUNT_CUTOFF_GUIDE_COUNTS,
        help="Cells with guide_count_col <= threshold treated as control (default: 0.0)",
    )
    ap.add_argument(
        "--cell-type-cols",
        default="adipo,pre_adipo,other,lipo",
        help="Comma-separated one-hot columns to derive cell_type (default: adipo,pre_adipo,other,lipo)",
    )
    ap.add_argument(
        "--default-cell-type",
        default="adipocyte",
        help="Fallback cell_type if one-hot columns missing (default: adipocyte)",
    )
    ap.add_argument("--ctrl-label", default="ctrl", help="Control label (default: ctrl)")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)

    if not in_path.exists():
        raise SystemExit(f"ERROR: input does not exist: {in_path}")
    if in_path.suffix.lower() != ".h5ad":
        raise SystemExit(f"ERROR: expected .h5ad input, got: {in_path}")

    adata = sc.read_h5ad(str(in_path))
    obs = adata.obs

    # Validate required source columns
    missing = []
    if args.gene_col not in obs.columns:
        missing.append(args.gene_col)
    if args.guide_count_col not in obs.columns:
        missing.append(args.guide_count_col)
    if missing:
        cols_preview = list(obs.columns)
        raise SystemExit(
            "ERROR: missing required adata.obs columns: "
            + ", ".join(missing)
            + "\nAvailable columns include: "
            + ", ".join(cols_preview[:60])
            + (" ..." if len(cols_preview) > 60 else "")
        )

    # Determine controls
    guide_counts = obs[args.guide_count_col].to_numpy()
    is_ctrl = guide_counts <= args.guide_threshold

    # Clean gene symbols
    genes = obs[args.gene_col].astype(str).to_numpy()
    genes_clean = np.where(
        (genes == "") | (genes == "nan") | (genes == "None") | (genes == "NA"),
        "unknown",
        genes,
    )

    # Build GEARS condition column
    #   - ctrl cells => "ctrl"
    #   - pert cells => "<GENE>+ctrl"
    obs["condition"] = np.where(is_ctrl, args.ctrl_label, genes_clean + f"+{args.ctrl_label}")

    # GEARS control flag
    obs["control"] = (obs["condition"] == args.ctrl_label).astype(int)

    # Dose values: match GEARS examples
    obs["dose_val"] = np.where(obs["control"] == 1, "1", "1+1")

    # cell_type: from one-hot if present, else constant
    onehot_cols = [c.strip() for c in args.cell_type_cols.split(",") if c.strip()]
    if onehot_cols and all(c in obs.columns for c in onehot_cols):
        obs["cell_type"] = derive_cell_type_from_onehot(obs, onehot_cols)
    else:
        obs["cell_type"] = args.default_cell_type

    # condition_name: handy unique label used across many GEARS flows
    obs["condition_name"] = (
        obs["cell_type"].astype(str) + "_" + obs["condition"].astype(str) + "_" + obs["dose_val"].astype(str)
    )

    # Sanity checks / summary
    print("GEARS fields added to adata.obs: condition, control, dose_val, cell_type, condition_name")
    print("\nTop conditions:")
    print(obs["condition"].value_counts().head(20))

    n_unknown = int((~is_ctrl & (genes_clean == "unknown")).sum())
    if n_unknown:
        print(
            f"\nWARNING: {n_unknown} non-control cells had missing/blank gene names -> condition='unknown+{args.ctrl_label}'. "
            "GEARS will likely drop these."
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(out_path))
    print(f"\nSaved GEARS-ready file: {out_path}")


if __name__ == "__main__":
    main()


