"""
Section 2B.1 — Assemble per-cell IKZF1(+) AUCell table
======================================================
Joins the aggregated pySCENIC AUCell loom (100K-cell subsample) to per-cell
donor / substate / diffusion-pseudotime labels from the trajectory object.

Inputs:
  results/phase2/GRN/scenic_auc_aggregated.loom   (ca: RegulonsAUC, CellID, state, ...)
  results/phase1/trajectory/microglia_trajectory.h5ad   (obs: state, Donor ID, dpt_pseudotime, Brain Region)

Output:
  results/phase7/downsample/ikzf1_auc_cells.parquet
    columns: cell_id, donor, state, region, ikzf1_auc, dpt

Barcode reconciliation: loom CellID and h5ad obs_names share BARCODE-LIBRARY
but differ in the trailing field (hash vs donor id). Join key = rsplit("-", 1)[0].
Aborts if overlap < 90%.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import loompy

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
LOOM = PROJECT / "results/phase2/GRN/scenic_auc_aggregated.loom"
H5AD = PROJECT / "results/phase1/trajectory/microglia_trajectory.h5ad"
OUT_DIR = PROJECT / "results/phase7/downsample"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "ikzf1_auc_cells.parquet"

REGULON = "IKZF1(+)"


def strip_key(x: str) -> str:
    return x.rsplit("-", 1)[0]


def main():
    print(f"Loading AUCell loom: {LOOM.name}")
    with loompy.connect(str(LOOM), mode="r", validate=False) as ds:
        auc = ds.ca["RegulonsAUC"]
        if REGULON not in auc.dtype.names:
            sys.exit(f"ERROR: {REGULON} not in loom regulons: {auc.dtype.names[:10]}")
        loom_df = pd.DataFrame({
            "cell_id": list(ds.ca["CellID"]),
            "ikzf1_auc": auc[REGULON].astype(float),
            "state_loom": list(ds.ca["state"]),
        })
    loom_df["key"] = loom_df["cell_id"].map(strip_key)
    print(f"  {len(loom_df):,} cells in loom")

    print(f"Loading trajectory obs: {H5AD.name}")
    A = ad.read_h5ad(H5AD, backed="r")
    obs = A.obs.copy()
    A.file.close()
    obs["key"] = pd.Index(obs.index).map(strip_key)
    lab = obs[["key", "state", "Donor ID", "Brain Region", "dpt_pseudotime"]].rename(
        columns={"state": "state", "Donor ID": "donor",
                 "Brain Region": "region", "dpt_pseudotime": "dpt"}
    )

    # Guard: key uniqueness on the label side
    dup = lab["key"].duplicated().sum()
    if dup:
        print(f"  WARNING: {dup} duplicate join keys in trajectory obs — keeping first")
        lab = lab.drop_duplicates("key", keep="first")

    merged = loom_df.merge(lab, on="key", how="left")
    overlap = merged["donor"].notna().mean()
    print(f"  join overlap: {overlap:.1%}")
    if overlap < 0.90:
        sys.exit(f"ERROR: join overlap {overlap:.1%} < 90% — barcode key mismatch")

    # Sanity: loom state vs trajectory state agreement
    agree = (merged["state_loom"].astype(str) == merged["state"].astype(str)).mean()
    print(f"  loom-state vs trajectory-state agreement: {agree:.1%}")

    out = merged.dropna(subset=["donor"])[
        ["cell_id", "donor", "state", "region", "ikzf1_auc", "dpt"]
    ].reset_index(drop=True)
    out["donor"] = out["donor"].astype(str)
    out["state"] = out["state"].astype(str)
    out["region"] = out["region"].astype(str)
    out.to_parquet(OUT, index=False)
    print(f"\nWrote {len(out):,} cells -> {OUT}")
    print("\nPer-state cell / donor counts:")
    summ = out.groupby("state").agg(
        n_cells=("cell_id", "size"),
        n_donors=("donor", "nunique"),
        median_cells_per_donor=("donor", lambda s: s.value_counts().median()),
        mean_auc=("ikzf1_auc", "mean"),
    )
    print(summ.to_string())


if __name__ == "__main__":
    main()
