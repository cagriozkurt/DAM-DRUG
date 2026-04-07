"""
DAM-DRUG Phase 1 — Trajectory Inference: PAGA + Diffusion Pseudotime
=====================================================================
Resumable: each expensive step saves a checkpoint h5ad.
On restart the script loads the most recent checkpoint and skips ahead.

Checkpoints (in results/phase1/trajectory/):
  _ckpt_01_neighbors.h5ad  — after sc.pp.neighbors  (~1.5 min)
  _ckpt_02_paga.h5ad       — after PAGA + graph plot (~10 s)
  _ckpt_03_umap.h5ad       — after PAGA-init UMAP   (~3.5 min)
  _ckpt_04_diffmap.h5ad    — after diffusion map     (~20 s)

Final outputs:
  results/phase1/trajectory/pseudotime.csv
  results/phase1/trajectory/paga_graph.pdf
  results/phase1/trajectory/umap_paga_pseudotime.pdf
  results/phase1/trajectory/microglia_trajectory.h5ad
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
RAW      = PROJECT / "data/raw/SEA-AD"
RES      = PROJECT / "results/phase1/trajectory"
RES.mkdir(parents=True, exist_ok=True)

H5AD = Path(os.environ.get(
    "SEA_AD_H5AD",
    str(RAW / "SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad")
))

# Checkpoints
CKPT = {
    1: RES / "_ckpt_01_neighbors.h5ad",
    2: RES / "_ckpt_02_paga.h5ad",
    3: RES / "_ckpt_03_umap.h5ad",
    4: RES / "_ckpt_04_diffmap.h5ad",
}

OUT_H5AD = RES / "microglia_trajectory.h5ad"
OUT_PT   = RES / "pseudotime.csv"

sc.settings.figdir   = str(RES)
sc.settings.verbosity = 2

SUPERTYPE_COL = "Supertype"
STATE_MAP = {
    "Micro-PVM_1":         "Homeostatic",
    "Micro-PVM_2":         "DAM",
    "Micro-PVM_2_3-SEAAD": "DAM-IRM",
    "Micro-PVM_2_1-SEAAD": "LateAD-DAM",
    "Micro-PVM_3-SEAAD":   "IRM",
    "Micro-PVM_4-SEAAD":   "LAM",
}
STATE_ORDER = ["Homeostatic", "DAM", "DAM-IRM", "IRM", "LateAD-DAM", "LAM"]
STATE_COLORS = {
    "Homeostatic": "#4477AA", "DAM": "#EE6677", "DAM-IRM": "#CCBB44",
    "IRM": "#228833", "LateAD-DAM": "#AA3377", "LAM": "#66CCEE",
}

N_NEIGHBORS = 30
N_DCS       = 15


def save_ckpt(mg, step):
    path = CKPT[step]
    log.info(f"  Saving checkpoint {step}: {path.name}...")
    mg.write_h5ad(path, compression="gzip")
    log.info(f"  Checkpoint {step} saved.")


def load_source():
    log.info(f"Loading {H5AD.name}...")
    mg = sc.read_h5ad(H5AD)
    log.info(f"  {mg.n_obs:,} nuclei × {mg.n_vars:,} genes")

    mg = mg[mg.obs[SUPERTYPE_COL].isin(STATE_MAP.keys())].copy()
    log.info(f"  After filtering: {mg.n_obs:,} microglia")

    mg.obs["state"] = pd.Categorical(
        mg.obs[SUPERTYPE_COL].map(STATE_MAP),
        categories=STATE_ORDER, ordered=True,
    )
    mg.uns["state_colors"] = [STATE_COLORS[s] for s in STATE_ORDER]

    if "UMIs" in mg.layers:
        mg.X = mg.layers["UMIs"].copy()
        log.info("  Using layers['UMIs'] as X")

    if "X_scVI" not in mg.obsm:
        log.error("X_scVI not in obsm. Available: " + str(list(mg.obsm.keys())))
        sys.exit(1)
    log.info(f"  X_scVI shape: {mg.obsm['X_scVI'].shape}")
    return mg


# ==============================================================================
def main():
    # Short-circuit if fully done
    if OUT_PT.exists() and OUT_H5AD.exists():
        log.info(f"Outputs already exist — nothing to do.\n  {OUT_PT}\n  {OUT_H5AD}")
        return

    # ── Determine resume point ─────────────────────────────────────────────────
    resume_from = 0
    for step in sorted(CKPT.keys(), reverse=True):
        if CKPT[step].exists():
            resume_from = step
            break

    if resume_from == 0:
        mg = load_source()
    else:
        log.info(f"Resuming from checkpoint {resume_from}: {CKPT[resume_from].name}")
        mg = sc.read_h5ad(CKPT[resume_from])
        log.info(f"  Loaded: {mg.n_obs:,} cells")

    # ── Step 1: Neighbors ──────────────────────────────────────────────────────
    if resume_from < 1:
        log.info(f"Step 1: Neighbors (k={N_NEIGHBORS}, use_rep=X_scVI)...")
        sc.pp.neighbors(mg, n_neighbors=N_NEIGHBORS, use_rep="X_scVI", random_state=42)
        save_ckpt(mg, 1)
    else:
        log.info("Step 1: Neighbors — skipped (checkpoint exists)")

    # ── Step 2: PAGA ──────────────────────────────────────────────────────────
    if resume_from < 2:
        log.info("Step 2: PAGA...")
        sc.tl.paga(mg, groups="state")
        conn = pd.DataFrame(
            mg.uns["paga"]["connectivities"].toarray(),
            index=STATE_ORDER, columns=STATE_ORDER,
        )
        log.info(f"  PAGA connectivity:\n{conn.round(3)}")
        sc.pl.paga(mg, color=["state"], layout="eq_tree",
                   fontsize=10, save="_graph.pdf", show=False)
        log.info(f"  PAGA graph saved.")
        save_ckpt(mg, 2)
    else:
        log.info("Step 2: PAGA — skipped (checkpoint exists)")

    # ── Step 3: PAGA-initialized UMAP ─────────────────────────────────────────
    if resume_from < 3:
        log.info("Step 3: PAGA-initialized UMAP...")
        sc.tl.umap(mg, init_pos="paga", random_state=42)
        save_ckpt(mg, 3)
    else:
        log.info("Step 3: UMAP — skipped (checkpoint exists)")

    # ── Step 4: Diffusion map ─────────────────────────────────────────────────
    if resume_from < 4:
        log.info("Step 4: Diffusion map...")
        sc.tl.diffmap(mg, n_comps=N_DCS)
        save_ckpt(mg, 4)
    else:
        log.info("Step 4: Diffmap — skipped (checkpoint exists)")

    # ── Step 5: DPT ───────────────────────────────────────────────────────────
    log.info("Step 5: Diffusion pseudotime...")
    hm_mask = np.asarray(mg.obs["state"] == "Homeostatic")

    if "P2RY12" in mg.var_names:
        gene_col = mg.X[:, mg.var_names.get_loc("P2RY12")]
        if hasattr(gene_col, "toarray"):
            gene_col = gene_col.toarray().ravel()
        else:
            gene_col = np.asarray(gene_col).ravel()
        hm_p2ry12  = gene_col[hm_mask]
        root_within = int(np.argmax(hm_p2ry12))
        root_global = int(np.where(hm_mask)[0][root_within])
        log.info(f"  Root: {mg.obs_names[root_global]} "
                 f"(Homeostatic, P2RY12={float(hm_p2ry12[root_within]):.1f})")
    else:
        root_global = int(np.where(hm_mask)[0][0])
        log.warning("P2RY12 not found — using first Homeostatic cell as root")

    mg.uns["iroot"] = root_global
    sc.tl.dpt(mg, n_dcs=N_DCS)

    pt_stats = mg.obs.groupby("state", observed=True)["dpt_pseudotime"].agg(["mean", "median"])
    log.info(f"  Pseudotime by state:\n{pt_stats.round(3)}")

    # ── Step 6: Plots ──────────────────────────────────────────────────────────
    log.info("Step 6: Saving plots...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(mg, color="state", ax=axes[0], show=False,
               title="Microglial states", frameon=False,
               palette=STATE_COLORS)
    sc.pl.umap(mg, color="dpt_pseudotime", ax=axes[1], show=False,
               title="Diffusion pseudotime", frameon=False, color_map="viridis")
    sc.pl.paga(mg, ax=axes[2], show=False, fontsize=9,
               title="PAGA graph", frameon=False)
    plt.tight_layout()
    fig.savefig(RES / "umap_paga_pseudotime.pdf", bbox_inches="tight")
    plt.close()
    log.info(f"  Plot saved: umap_paga_pseudotime.pdf")

    # ── Step 7: Export ────────────────────────────────────────────────────────
    log.info("Step 7: Exporting pseudotime.csv and final h5ad...")
    pt_df = mg.obs[["state", "dpt_pseudotime"]].copy()
    pt_df.index.name = "cell_id"
    pt_df.to_csv(OUT_PT)
    log.info(f"  pseudotime.csv: {len(pt_df):,} cells → {OUT_PT}")

    del mg.layers
    mg.write_h5ad(OUT_H5AD, compression="gzip")
    log.info(f"  trajectory h5ad → {OUT_H5AD}")

    # Clean up checkpoints
    for ckpt in CKPT.values():
        if ckpt.exists():
            ckpt.unlink()
    log.info("  Checkpoints cleaned up.")

    log.info("\nDone. Next: run 2.7 regulon-pseudotime correlation.")


if __name__ == "__main__":
    main()
