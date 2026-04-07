"""
Prepare per-state scRNA + ATAC MTX input for scMultiomeGRN.

Memory strategy:
  1. Load RNA h5ad (already microglia-only, ~3 GB)
  2. Open ATAC h5ad in backed (lazy) mode to read obs only
  3. Intersect barcodes → load only microglia ATAC cells into RAM
  4. Iterate states, write MTX; skip states whose output already exists

Outputs per state (e.g. DAM/):
  scrna/matrix.mtx + barcodes.tsv + genes.tsv + var_features.tsv
  atac/matrix.mtx  + barcodes.tsv + peaks.tsv  (names: chr_start_end)
"""

import os
import re
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import scanpy as sc
from pathlib import Path

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
TRAJ    = PROJECT / "results/phase1/trajectory"
OUT     = PROJECT / "results/phase2/scMultiomeGRN/input"

RNA_H5AD  = TRAJ / "microglia_trajectory.h5ad"
ATAC_H5AD = Path(os.environ.get("ATAC_H5AD",
    str(PROJECT / "data/raw/SEA-AD/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad")))

STATES = ["Homeostatic", "DAM", "DAM-IRM", "IRM", "LateAD-DAM", "LAM"]
N_HVG  = 2000


def to_underscore(name: str) -> str:
    """chr1:12345-67890 → chr1_12345_67890 (required by scMultiomeGRN)."""
    return re.sub(r"[:\-]", "_", str(name))


def state_done(state: str) -> bool:
    """Return True if both MTX files already exist for this state."""
    tag = state.replace("-", "_")
    return (
        (OUT / tag / "scrna" / "matrix.mtx").exists() and
        (OUT / tag / "atac"  / "matrix.mtx").exists()
    )


def export_state(rna, atac, state: str):
    tag       = state.replace("-", "_")
    scrna_dir = OUT / tag / "scrna"
    atac_dir  = OUT / tag / "atac"
    scrna_dir.mkdir(parents=True, exist_ok=True)
    atac_dir.mkdir(parents=True, exist_ok=True)

    # ── RNA subset ────────────────────────────────────────────────────────────
    rna_s = rna[rna.obs["state"] == state].copy()
    print(f"  RNA  : {rna_s.n_obs:,} cells × {rna_s.n_vars:,} genes")

    X_rna = rna_s.X.T
    if not sp.issparse(X_rna):
        X_rna = sp.csr_matrix(X_rna)
    scipy.io.mmwrite(str(scrna_dir / "matrix.mtx"), X_rna.astype(np.float32))
    pd.Series(rna_s.obs_names).to_csv(scrna_dir / "barcodes.tsv",    index=False, header=False)
    pd.Series(rna_s.var_names).to_csv(scrna_dir / "genes.tsv",       index=False, header=False)

    if "highly_variable" in rna_s.var.columns:
        hvg = rna_s.var_names[rna_s.var["highly_variable"]].tolist()[:N_HVG]
    else:
        means = np.asarray(rna_s.X.mean(axis=0)).ravel()
        hvg = rna_s.var_names[np.argsort(means)[::-1][:N_HVG]].tolist()
    pd.Series(hvg).to_csv(scrna_dir / "var_features.tsv", index=False, header=False)
    print(f"  HVG  : {len(hvg)}")

    # ── ATAC subset ───────────────────────────────────────────────────────────
    shared = rna_s.obs_names.intersection(atac.obs_names)
    if len(shared) == 0:
        print(f"  WARNING: no shared barcodes — skipping ATAC for {state}")
        return
    atac_s = atac[shared].copy()
    print(f"  ATAC : {atac_s.n_obs:,} cells × {atac_s.n_vars:,} peaks")

    X_atac = atac_s.X.T
    if not sp.issparse(X_atac):
        X_atac = sp.csr_matrix(X_atac)
    scipy.io.mmwrite(str(atac_dir / "matrix.mtx"), X_atac.astype(np.float32))
    pd.Series(atac_s.obs_names).to_csv(atac_dir / "barcodes.tsv", index=False, header=False)
    peak_names = [to_underscore(p) for p in atac_s.var_names]
    pd.Series(peak_names).to_csv(atac_dir / "peaks.tsv", index=False, header=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # ── Check which states still need processing ───────────────────────────────
    todo = [s for s in STATES if not state_done(s)]
    done = [s for s in STATES if     state_done(s)]
    if done:
        print(f"Already done (skipping): {done}")
    if not todo:
        print("All states already exported. Nothing to do.")
        sys.exit(0)
    print(f"To export: {todo}")

    # ── Load RNA ─────────────────────────────────────────────────────────────
    print(f"\nLoading RNA: {RNA_H5AD}")
    rna = sc.read_h5ad(RNA_H5AD)
    print(f"  {rna.n_obs:,} cells × {rna.n_vars:,} genes")

    # ── Load ATAC: backed mode first to get obs, then subset into RAM ─────────
    if not ATAC_H5AD.exists():
        print(f"ERROR: ATAC h5ad not found: {ATAC_H5AD}")
        sys.exit(1)

    print(f"\nOpening ATAC (backed): {ATAC_H5AD}")
    atac_backed = sc.read_h5ad(ATAC_H5AD, backed='r')
    print(f"  Total: {atac_backed.n_obs:,} cells × {atac_backed.n_vars:,} peaks")

    # Intersect with all microglia barcodes in the todo states
    todo_barcodes = rna.obs_names[rna.obs["state"].isin(todo)]
    shared_all = todo_barcodes.intersection(atac_backed.obs_names)
    print(f"  Shared with todo-state microglia: {len(shared_all):,} cells")

    if len(shared_all) == 0:
        print("ERROR: No shared barcodes. Check that RNA and ATAC h5ads come from the same dataset.")
        atac_backed.file.close()
        sys.exit(1)

    # Load only the shared cells into RAM
    print("  Loading shared cells into RAM...")
    atac = atac_backed[shared_all].to_memory()
    atac_backed.file.close()
    print(f"  ATAC in RAM: {atac.n_obs:,} × {atac.n_vars:,}")

    # ── Export per state ───────────────────────────────────────────────────────
    for state in todo:
        if state not in rna.obs["state"].values:
            print(f"\nSkip {state} — not in RNA obs")
            continue
        print(f"\n=== {state} ===")
        export_state(rna, atac, state)

    print(f"\nDone. Input files in: {OUT}")


if __name__ == "__main__":
    main()
