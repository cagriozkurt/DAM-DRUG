"""
DAM-DRUG Phase 2 — pySCENIC Gene Regulatory Network Inference
=============================================================
Uses pySCENIC CLI (not Python API) to avoid dask version incompatibilities.

Pipeline:
  Step 0: Export microglia h5ad → loom
  Step 1: pyscenic grn   (GRNBoost2 co-expression)
  Step 2: pyscenic ctx   (RcisTarget motif pruning → regulons)
  Step 3: pyscenic aucell (AUCell regulon activity scoring)
  Step 4: Summarize regulon activity by microglial state

Resources required (download on login node):
  data/resources/motifs/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  data/resources/motifs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
  data/resources/tf_lists/allTFs_hg38.txt
"""

import os
import sys
import subprocess
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import loompy
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
RAW      = PROJECT / "data/raw/SEA-AD"
RESOURCE = PROJECT / "data/resources"
RES      = PROJECT / "results/phase2/GRN"
RES.mkdir(parents=True, exist_ok=True)

H5AD = Path(os.environ.get(
    "SEA_AD_H5AD",
    str(RAW / "SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad")
))

RANKINGS  = RESOURCE / "motifs/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
MOTIF_TBL = RESOURCE / "motifs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TF_LIST   = RESOURCE / "tf_lists/allTFs_hg38.txt"

# Intermediate / output files
LOOM_IN   = RES / "microglia_raw.loom"
ADJ_CSV   = RES / "adj_matrix.tsv"
REG_CSV   = RES / "regulons.csv"
AUC_LOOM  = RES / "scenic_auc.loom"

N_WORKERS = int(os.environ.get("SLURM_NTASKS", 20))

SUPERTYPE_COL = "Supertype"
STATE_MAP = {
    "Micro-PVM_1":         "Homeostatic",
    "Micro-PVM_2":         "DAM",
    "Micro-PVM_2_3-SEAAD": "DAM-IRM",
    "Micro-PVM_2_1-SEAAD": "LateAD-DAM",
    "Micro-PVM_3-SEAAD":   "IRM",
    "Micro-PVM_4-SEAAD":   "LAM",
}
TF_TARGETS = ["SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
              "RELB", "BHLHE40", "BHLHE41"]


# ==============================================================================
def check_resources():
    missing = [str(f) for f in [RANKINGS, MOTIF_TBL, TF_LIST] if not f.exists()]
    if missing:
        log.error("Missing resource files:\n  " + "\n  ".join(missing))
        sys.exit(1)
    log.info("All resource files present.")


def run_cmd(cmd, step_name):
    log.info(f"Running {step_name}:\n  {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        log.error(f"{step_name} failed with exit code {result.returncode}")
        sys.exit(result.returncode)
    log.info(f"{step_name} complete.")


# ==============================================================================
# Step 0: Load h5ad, filter, export to loom
# ==============================================================================
def prepare_loom():
    if LOOM_IN.exists():
        log.info(f"Loom already exists: {LOOM_IN}")
        return

    log.info(f"Loading {H5AD.name}...")
    mg = sc.read_h5ad(H5AD)
    log.info(f"  {mg.n_obs:,} nuclei × {mg.n_vars:,} genes")

    # Filter to microglia only
    mg = mg[mg.obs[SUPERTYPE_COL].isin(list(STATE_MAP.keys()))].copy()
    log.info(f"  After filtering: {mg.n_obs:,} cells")

    # Use raw UMI counts
    if "UMIs" in mg.layers:
        mg.X = mg.layers["UMIs"].copy()
        log.info("  Using layers['UMIs'] as raw counts")

    # Filter low-expressed genes
    sc.pp.filter_genes(mg, min_cells=10)
    log.info(f"  After gene filter: {mg.n_vars:,} genes")

    # Add state label
    mg.obs["state"] = mg.obs[SUPERTYPE_COL].map(STATE_MAP)

    # Subsample for GRNBoost2 — 100K cells is standard practice; reduces dask graph to ~25 GB
    # Simple random draw without replacement; NOT stratified by state or donor.
    # Substate representation is proportional by chance (see manuscript Methods §GRN inference).
    N_GRN = 100_000
    if mg.n_obs > N_GRN:
        sc.pp.subsample(mg, n_obs=N_GRN, random_state=42)
        log.info(f"  Subsampled to {mg.n_obs:,} cells (uniform random, random_state=42)")

    # Keep only essential obs columns — loompy rejects column names with "/"
    keep_cols = ["state", SUPERTYPE_COL, "Donor ID", "Braak Stage", "APOE Genotype", "Brain Region"]
    mg.obs = mg.obs[[c for c in keep_cols if c in mg.obs.columns]].copy()
    # Sanitize any remaining special characters in column names
    mg.obs.columns = [c.replace("/", "_").replace(" ", "_") for c in mg.obs.columns]

    # pySCENIC CLI requires row attribute "Gene" and column attribute "CellID"
    mg.var.index.name = "Gene"
    mg.obs.index.name = "CellID"

    # Export to loom (pySCENIC CLI expects genes as rows)
    log.info(f"Writing loom: {LOOM_IN}")
    mg.write_loom(str(LOOM_IN), write_obsm_varm=False)
    log.info(f"  Loom written: {mg.n_obs:,} cells × {mg.n_vars:,} genes")

    return mg


# ==============================================================================
# Step 1: GRNBoost2 co-expression (pyscenic grn)
# ==============================================================================
def run_grn():
    if ADJ_CSV.exists():
        log.info(f"Adjacency matrix exists, skipping GRN step: {ADJ_CSV}")
        return

    run_cmd([
        "pyscenic", "grn",
        str(LOOM_IN),
        str(TF_LIST),
        "--output", str(ADJ_CSV),
        "--num_workers", str(N_WORKERS),
        "--seed", "42",
    ], "GRNBoost2")


# ==============================================================================
# Step 2: RcisTarget motif pruning (pyscenic ctx)
# ==============================================================================
def run_ctx():
    if REG_CSV.exists():
        log.info(f"Regulons exist, skipping CTX step: {REG_CSV}")
        return

    run_cmd([
        "pyscenic", "ctx",
        str(ADJ_CSV),
        str(RANKINGS),
        "--annotations_fname", str(MOTIF_TBL),
        "--expression_mtx_fname", str(LOOM_IN),
        "--output", str(REG_CSV),
        "--num_workers", str(N_WORKERS),
        "--mode", "dask_multiprocessing",
        "--mask_dropouts",
    ], "RcisTarget")


# ==============================================================================
# Step 3: AUCell regulon activity scoring (pyscenic aucell)
# ==============================================================================
def run_aucell():
    if AUC_LOOM.exists():
        log.info(f"AUC loom exists, skipping AUCell step: {AUC_LOOM}")
        return

    run_cmd([
        "pyscenic", "aucell",
        str(LOOM_IN),
        str(REG_CSV),
        "--output", str(AUC_LOOM),
        "--num_workers", str(N_WORKERS),
    ], "AUCell")


# ==============================================================================
# Step 4: Summarize regulon activity by state
# ==============================================================================
def summarize():
    if not AUC_LOOM.exists():
        log.error("AUC loom not found — run steps 1-3 first")
        return

    log.info("Summarizing regulon activity by microglial state...")
    with loompy.connect(str(AUC_LOOM), mode="r", validate=False) as ds:
        # pySCENIC writes RegulonsAUC as a structured numpy array (one field per regulon)
        auc_raw = ds.ca["RegulonsAUC"]
        regulon_names = list(auc_raw.dtype.names)
        auc_mtx = pd.DataFrame(
            {name: auc_raw[name] for name in regulon_names},
            index=ds.ca["CellID"],
        )
        cell_states = pd.Series(ds.ca["state"], index=ds.ca["CellID"])

    auc_mtx["state"] = cell_states
    state_mean = auc_mtx.groupby("state").mean()
    state_mean.to_csv(RES / "regulon_auc_by_state.csv")
    log.info(f"  State AUC matrix: {state_mean.shape}")

    # TF target summary
    log.info(f"\n{'='*60}")
    log.info("TF TARGET REGULON SUMMARY")
    log.info(f"{'='*60}")

    summary_rows = []
    for tf in TF_TARGETS:
        matches = [c for c in state_mean.columns if c.startswith(f"{tf}(")]
        if matches:
            for reg in matches:
                aucs = state_mean[reg]
                peak = aucs.idxmax()
                row = {"TF": tf, "regulon": reg,
                       "peak_state": peak, "peak_AUC": aucs.max()}
                for state in state_mean.index:
                    row[f"AUC_{state}"] = aucs[state]
                summary_rows.append(row)
                log.info(f"  {reg:<35} peak={peak} ({aucs.max():.4f})")
        else:
            log.info(f"  {tf:<35} NO REGULON FOUND")
            summary_rows.append({"TF": tf, "regulon": "NOT_FOUND"})

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(RES / "TF_regulon_summary.csv", index=False)
    log.info(f"\n  Saved: {RES / 'TF_regulon_summary.csv'}")
    log.info(f"Done. All results in: {RES}")


# ==============================================================================
if __name__ == "__main__":
    check_resources()
    prepare_loom()
    run_grn()
    run_ctx()
    run_aucell()
    summarize()
    log.info("Next step: sbatch code/slurm/07_grn_multiseed.slurm (5-seed GRNBoost2 on HPC)")
