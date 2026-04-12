"""
DAM-DRUG Supplementary Figure S3 — QC & Dataset Overview
==========================================================
Note: Per-dataset QC filtering was skipped (step 1.6) — the SEA-AD
pre-release h5ad is already QC'd by the Allen Institute. Panels show
the post-QC distribution from microglia_trajectory.h5ad.

Panels:
  A. Violin: nGenes per cell by microglial state
  B. Violin: nUMI per cell by microglial state
  C. Violin: % mitochondrial reads by microglial state
  D. UMAP colored by donor_id (batch overview)
  E. UMAP colored by microglial state
  F. Stacked bar: cell counts per state by Braak stage

Requires: microglia_trajectory.h5ad (TRUBA)

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
import matplotlib.patches as mpatches
from pathlib import Path

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
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
adata = sc.read_h5ad(TRAJ / "microglia_trajectory.h5ad")
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Identify obs columns (column names vary across SEA-AD versions)
state_col  = next((c for c in ["state", "microglia_state", "subtype", "leiden"]
                   if c in adata.obs.columns), None)
donor_col  = next((c for c in ["Donor ID", "donor_id", "donor", "batch"]
                   if c in adata.obs.columns), None)
braak_col  = next((c for c in ["Braak", "Braak stage", "braak_stage", "BRAAK"]
                   if c in adata.obs.columns), None)
ngenes_col = next((c for c in ["Genes detected", "n_genes_by_counts", "n_genes", "nFeature_RNA"]
                   if c in adata.obs.columns), None)
numi_col   = next((c for c in ["Number of UMIs", "total_counts", "n_counts", "nCount_RNA"]
                   if c in adata.obs.columns), None)
mito_col   = next((c for c in ["Fraction mitochondrial UMIs", "pct_counts_mt", "pct_counts_mito"]
                   if c in adata.obs.columns), None)

print(f"state_col={state_col}, donor_col={donor_col}, braak_col={braak_col}")
print(f"ngenes_col={ngenes_col}, numi_col={numi_col}, mito_col={mito_col}")

# ── Figure layout ─────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 3, hspace=0.45, wspace=0.35, figure=fig)
axes_top = [fig.add_subplot(gs[0, i]) for i in range(3)]
axes_bot = [fig.add_subplot(gs[1, i]) for i in range(3)]

# ── Panels A–C: QC metric violins by state ────────────────────────────────
QC_DEFS = [
    (ngenes_col, "# Genes per cell",   "A"),
    (numi_col,   "UMI per cell",        "B"),
    (mito_col,   "% Mitochondrial",     "C"),
]

for ax, (metric, ylabel, label) in zip(axes_top, QC_DEFS):
    ax.text(-0.15, 1.06, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold")

    if metric is None or state_col is None:
        ax.text(0.5, 0.5, f"Column not found\nin h5ad obs",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color="#999999", style="italic")
        ax.set_title(ylabel, fontsize=8)
        continue

    data_by_state = []
    positions     = []
    tick_labels   = []
    for i, state in enumerate(STATE_ORDER):
        mask = adata.obs[state_col] == state
        if mask.sum() == 0:
            continue
        vals = adata.obs.loc[mask, metric].dropna().values
        data_by_state.append(vals)
        positions.append(i)
        tick_labels.append(state)

    if not data_by_state:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
        continue

    vp = ax.violinplot(data_by_state, positions=positions,
                       showmedians=True, widths=0.7)
    colors = [STATE_COLORS.get(s, "#AAAAAA") for s in tick_labels]
    for body, color in zip(vp["bodies"], colors):
        body.set_facecolor(color)
        body.set_alpha(0.7)
        body.set_edgecolor("none")
    vp["cmedians"].set_colors(["#333333"] * len(data_by_state))
    for part in ("cbars", "cmins", "cmaxes"):
        vp[part].set_colors(["#666666"] * len(data_by_state))
        vp[part].set_linewidth(0.7)

    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels, rotation=35, ha="right", fontsize=7)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(ylabel, fontsize=8, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # Annotate n per state
    for pos, vals in zip(positions, data_by_state):
        ax.text(pos, ax.get_ylim()[0], f"n={len(vals):,}",
                ha="center", va="bottom", fontsize=5.5, color="#555555", rotation=90)

# ── Panel D: UMAP by donor ────────────────────────────────────────────────
ax_d  = axes_bot[0]
ax_d.text(-0.08, 1.04, "D", transform=ax_d.transAxes, fontsize=11, fontweight="bold")
coords = adata.obsm["X_umap"]

if donor_col is not None:
    donors = adata.obs[donor_col].astype(str)
    unique_donors = sorted(donors.unique())
    cmap = matplotlib.colormaps["tab20"].resampled(len(unique_donors))
    for i, d in enumerate(unique_donors):
        mask = donors == d
        ax_d.scatter(coords[mask, 0], coords[mask, 1],
                     s=0.2, alpha=0.3, c=[cmap(i)], rasterized=True)
    ax_d.set_title(f"UMAP by donor (n={len(unique_donors)} donors)",
                   fontsize=8, fontweight="bold")
else:
    ax_d.scatter(coords[:, 0], coords[:, 1], s=0.2, alpha=0.3,
                 c="#AAAAAA", rasterized=True)
    ax_d.set_title("UMAP (donor column not found)", fontsize=8)
ax_d.axis("off")

# ── Panel E: UMAP by microglial state ────────────────────────────────────
ax_e = axes_bot[1]
ax_e.text(-0.08, 1.04, "E", transform=ax_e.transAxes, fontsize=11, fontweight="bold")

if state_col is not None:
    states = adata.obs[state_col].astype(str)
    for state in STATE_ORDER:
        mask = states == state
        if mask.sum() == 0:
            continue
        ax_e.scatter(coords[mask, 0], coords[mask, 1],
                     s=0.5, alpha=0.6, c=STATE_COLORS.get(state, "grey"),
                     label=state, rasterized=True)
    legend_handles = [mpatches.Patch(color=STATE_COLORS.get(s, "grey"), label=s)
                      for s in STATE_ORDER if (adata.obs[state_col] == s).any()]
    ax_e.legend(handles=legend_handles, markerscale=1, fontsize=6.5,
                loc="lower left", framealpha=0.85, handlelength=1.0)
else:
    ax_e.scatter(coords[:, 0], coords[:, 1], s=0.3, alpha=0.4,
                 c="#0072B2", rasterized=True)
ax_e.set_title("UMAP by microglial state", fontsize=8, fontweight="bold")
ax_e.axis("off")

# ── Panel F: Cell counts by state (± by Braak if available) ──────────────
ax_f = axes_bot[2]
ax_f.text(-0.08, 1.04, "F", transform=ax_f.transAxes, fontsize=11, fontweight="bold")

if state_col is not None and braak_col is not None:
    # Stacked bar: Braak stage × state
    obs = adata.obs[[state_col, braak_col]].copy()
    obs[braak_col] = obs[braak_col].astype(str)
    ct = obs.groupby([braak_col, state_col]).size().unstack(fill_value=0)
    ct = ct.reindex(columns=[s for s in STATE_ORDER if s in ct.columns])
    bottom = np.zeros(len(ct))
    braak_labels = ct.index.tolist()
    for state in ct.columns:
        vals = ct[state].values
        ax_f.bar(braak_labels, vals, bottom=bottom,
                 color=STATE_COLORS.get(state, "grey"), label=state,
                 edgecolor="white", linewidth=0.3)
        bottom += vals
    ax_f.set_xlabel("Braak stage", fontsize=8)
    ax_f.set_ylabel("Cell count", fontsize=8)
    ax_f.set_title("Cell counts by Braak stage", fontsize=8, fontweight="bold")
    ax_f.tick_params(axis="x", labelsize=7)
elif state_col is not None:
    # Simple bar: total cells per state
    counts = adata.obs[state_col].value_counts().reindex(STATE_ORDER, fill_value=0)
    ax_f.bar(counts.index, counts.values,
             color=[STATE_COLORS.get(s, "grey") for s in counts.index],
             edgecolor="white", linewidth=0.4)
    ax_f.set_xlabel("State", fontsize=8)
    ax_f.set_ylabel("Cell count", fontsize=8)
    ax_f.set_title("Cell counts per microglial state", fontsize=8, fontweight="bold")
    ax_f.tick_params(axis="x", rotation=35, labelsize=7)
else:
    ax_f.text(0.5, 0.5, "State column not found",
              transform=ax_f.transAxes, ha="center", va="center", fontsize=8)

ax_f.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    f"Supplementary Figure S3 — scRNA-seq Dataset Overview & QC "
    f"({adata.shape[0]:,} nuclei, SEA-AD pre-release)",
    fontsize=10
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"supp_fig_S3_qc.{ext}",
                bbox_inches="tight", dpi=300)
    print(f"Saved → {OUT}/supp_fig_S3_qc.{ext}")

plt.close(fig)
