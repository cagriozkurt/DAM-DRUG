"""
DAM-DRUG Supplementary Figure S3 — QC & Batch Correction
==========================================================
Panels:
  A. Violin: nGenes per cell before/after QC, by sample
  B. Violin: nUMI per cell before/after QC, by sample
  C. Violin: % mitochondrial reads before/after QC
  D. UMAP before Harmony integration (colored by donor/batch)
  E. UMAP after Harmony integration (colored by donor/batch)
  F. UMAP after integration, colored by microglial state

Requires: microglia_raw.h5ad + microglia_trajectory.h5ad (TRUBA)

Run on TRUBA:
  apptainer exec --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \\
    ~/containers/scanpy-env.sif \\
    conda run -n scenic python code/phase6_figures/supp_fig_s3_qc.py
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
EXPLORE  = PROJECT / "results/phase1/explore"
TRAJ     = PROJECT / "results/phase1/trajectory"
OUT      = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})

STATE_ORDER  = ["Homeostatic", "DAM", "DAM-IRM", "IRM", "LAM", "LateAD-DAM"]
STATE_COLORS = {"Homeostatic": "#0072B2", "DAM": "#E69F00", "DAM-IRM": "#F0E442",
                "IRM": "#56B4E9", "LAM": "#CC79A7", "LateAD-DAM": "#D55E00"}

# ── Load data ──────────────────────────────────────────────────────────────
# Pre-QC: raw counts with QC metrics computed
adata_raw  = sc.read_h5ad(EXPLORE / "microglia_raw_qc.h5ad")
# Post-QC, post-Harmony: final trajectory object
adata_traj = sc.read_h5ad(TRAJ / "microglia_trajectory.h5ad")

# ── Figure layout ─────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 3, hspace=0.4, wspace=0.35, figure=fig)

axes_top = [fig.add_subplot(gs[0, i]) for i in range(3)]
axes_bot = [fig.add_subplot(gs[1, i]) for i in range(3)]

QC_METRICS = [
    ("n_genes_by_counts", "# Genes per cell", "A"),
    ("total_counts",      "UMI per cell",      "B"),
    ("pct_counts_mt",     "% Mitochondrial",   "C"),
]

# ── Panels A–C: QC violins before/after filtering ────────────────────────
for ax, (metric, ylabel, label) in zip(axes_top, QC_METRICS):
    # Raw (before QC) and filtered data
    pre_vals  = adata_raw.obs[metric].values
    post_vals = adata_traj.obs[metric].values if metric in adata_traj.obs else np.array([])

    data_violin = []
    labels_v    = []
    if len(pre_vals) > 0:
        data_violin.append(pre_vals)
        labels_v.append("Pre-QC")
    if len(post_vals) > 0:
        data_violin.append(post_vals)
        labels_v.append("Post-QC")

    vp = ax.violinplot(data_violin, positions=range(len(data_violin)),
                       showmedians=True, widths=0.7)
    for i, (body, color) in enumerate(zip(vp["bodies"], ["#AAAAAA", "#0072B2"])):
        body.set_facecolor(color)
        body.set_alpha(0.7)
    vp["cmedians"].set_colors(["#333333"] * len(data_violin))

    ax.set_xticks(range(len(labels_v)))
    ax.set_xticklabels(labels_v, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.15, 1.05, label, transform=ax.transAxes, fontsize=11, fontweight="bold")

    # Annotate median
    for i, vals in enumerate(data_violin):
        ax.text(i, np.median(vals), f" n={len(vals):,}", va="center",
                fontsize=6, color="#333333")

# ── Panel D: UMAP pre-Harmony (colored by sample/donor) ─────────────────
ax_d = axes_bot[0]
if "X_umap_preharmony" in adata_traj.obsm:
    coords = adata_traj.obsm["X_umap_preharmony"]
    donors = adata_traj.obs["donor_id"] if "donor_id" in adata_traj.obs.columns else adata_traj.obs.iloc[:, 0]
    unique_donors = donors.unique()
    cmap = plt.cm.get_cmap("tab20", len(unique_donors))
    for i, d in enumerate(unique_donors):
        mask = donors == d
        ax_d.scatter(coords[mask, 0], coords[mask, 1],
                     s=0.3, alpha=0.4, c=[cmap(i)], rasterized=True)
else:
    ax_d.text(0.5, 0.5, "X_umap_preharmony\nnot found in h5ad",
              transform=ax_d.transAxes, ha="center", va="center", fontsize=8)
ax_d.set_title("Pre-Harmony UMAP\n(by donor)", fontsize=8, fontweight="bold")
ax_d.axis("off")
ax_d.text(-0.05, 1.04, "D", transform=ax_d.transAxes, fontsize=11, fontweight="bold")

# ── Panel E: UMAP post-Harmony (colored by donor) ────────────────────────
ax_e = axes_bot[1]
coords = adata_traj.obsm["X_umap"]
donors = adata_traj.obs["donor_id"] if "donor_id" in adata_traj.obs.columns else adata_traj.obs.iloc[:, 0]
unique_donors = donors.unique()
cmap = plt.cm.get_cmap("tab20", len(unique_donors))
for i, d in enumerate(unique_donors):
    mask = donors == d
    ax_e.scatter(coords[mask, 0], coords[mask, 1],
                 s=0.3, alpha=0.4, c=[cmap(i)], rasterized=True)
ax_e.set_title("Post-Harmony UMAP\n(by donor)", fontsize=8, fontweight="bold")
ax_e.axis("off")
ax_e.text(-0.05, 1.04, "E", transform=ax_e.transAxes, fontsize=11, fontweight="bold")

# ── Panel F: UMAP by microglial state ────────────────────────────────────
ax_f = axes_bot[2]
state_col = "microglia_state" if "microglia_state" in adata_traj.obs.columns else "leiden"
states = adata_traj.obs[state_col]
for state in STATE_ORDER:
    mask = states == state
    if mask.sum() == 0:
        continue
    ax_f.scatter(coords[mask, 0], coords[mask, 1],
                 s=0.5, alpha=0.6, c=STATE_COLORS.get(state, "grey"),
                 label=state, rasterized=True)
ax_f.legend(markerscale=6, fontsize=6.5, loc="lower left",
            framealpha=0.8, handlelength=0.8)
ax_f.set_title("Post-Harmony UMAP\n(by microglial state)", fontsize=8, fontweight="bold")
ax_f.axis("off")
ax_f.text(-0.05, 1.04, "F", transform=ax_f.transAxes, fontsize=11, fontweight="bold")

fig.suptitle(
    "Supplementary Figure S3 — scRNA-seq Quality Control and Batch Correction",
    fontsize=10
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"supp_fig_S3_qc.{ext}",
                bbox_inches="tight", dpi=300)
    print(f"Saved → {OUT}/supp_fig_S3_qc.{ext}")

plt.close(fig)
