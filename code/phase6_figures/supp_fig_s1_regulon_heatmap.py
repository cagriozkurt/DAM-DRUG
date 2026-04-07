"""
DAM-DRUG Supplementary Figure S1 — Full 46-Regulon AUCell Heatmap
==================================================================
Panel: Heatmap of mean AUCell scores for all 46 pySCENIC regulons
       across 6 microglial states (rows=states, cols=regulons).
       Z-score normalised across states per regulon; CVD palette.

Runs locally (no h5ad needed):
  python code/phase6_figures/supp_fig_s1_regulon_heatmap.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from pathlib import Path

PROJECT = Path("/Volumes/PortableSSD/untitled folder/DAM-DRUG")
GRN     = PROJECT / "results/phase2/GRN"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.linewidth": 0.6, "pdf.fonttype": 42,
})

STATE_ORDER  = ["Homeostatic", "DAM", "DAM-IRM", "IRM", "LAM", "LateAD-DAM"]
STATE_COLORS = {"Homeostatic": "#0072B2", "DAM": "#E69F00", "DAM-IRM": "#F0E442",
                "IRM": "#56B4E9", "LAM": "#CC79A7", "LateAD-DAM": "#D55E00"}

# ── Load AUCell scores ─────────────────────────────────────────────────────
auc = pd.read_csv(GRN / "regulon_auc_by_state_aggregated.csv", index_col=0)
auc = auc.reindex(STATE_ORDER)           # enforce state order
auc.columns = [c.replace("(+)", "")      # strip (+) suffix for cleaner labels
               for c in auc.columns]

# Z-score per regulon (column-wise)
from scipy.stats import zscore as scipy_zscore
auc_z = auc.apply(scipy_zscore, axis=0)  # each column normalised to μ=0, σ=1

# ── Cluster regulons by similarity ────────────────────────────────────────
# Use scipy linkage on regulon columns (similar states will co-cluster)
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

col_linkage = linkage(pdist(auc_z.T, "euclidean"), method="average")
col_order   = leaves_list(col_linkage)
auc_z_sorted = auc_z.iloc[:, col_order]

# ── Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    2, 1,
    figsize=(14, 4.5),
    gridspec_kw={"height_ratios": [0.12, 1], "hspace": 0.02}
)

# Top strip: state colour bars (one coloured cell per state)
strip_ax = axes[0]
state_rgb = np.array([[mcolors.to_rgb(STATE_COLORS[s]) for s in STATE_ORDER]])
strip_ax.imshow(state_rgb, aspect="auto")
strip_ax.set_xticks([])
strip_ax.set_yticks([0])
strip_ax.set_yticklabels(["State"], fontsize=7)
for spine in strip_ax.spines.values():
    spine.set_visible(False)

# Heatmap
heat_ax = axes[1]
im = heat_ax.imshow(
    auc_z_sorted.values,
    aspect="auto",
    cmap="RdBu_r",
    vmin=-2, vmax=2,
    interpolation="nearest",
)

# Axes labels
heat_ax.set_yticks(range(len(STATE_ORDER)))
heat_ax.set_yticklabels(STATE_ORDER, fontsize=8)
heat_ax.set_xticks(range(auc_z_sorted.shape[1]))
heat_ax.set_xticklabels(auc_z_sorted.columns, rotation=90, fontsize=6.5, ha="center")
heat_ax.set_xlabel("Regulon (TF)", fontsize=9, labelpad=4)

# Highlight key TFs in x-tick labels
key_tfs = {"IKZF1"}
for lbl in heat_ax.get_xticklabels():
    if lbl.get_text() in key_tfs:
        lbl.set_fontweight("bold")
        lbl.set_color("#D55E00")

# Colorbar
cbar = fig.colorbar(im, ax=heat_ax, orientation="vertical",
                    fraction=0.015, pad=0.01, shrink=0.8)
cbar.set_label("AUCell z-score", fontsize=8)
cbar.ax.tick_params(labelsize=7)

fig.suptitle(
    "Supplementary Figure S1 — pySCENIC 46-Regulon AUCell Scores Across Microglial States",
    fontsize=10, y=1.01
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"supp_fig_S1_regulon_heatmap.{ext}",
                bbox_inches="tight", dpi=300)
    print(f"Saved → {OUT}/supp_fig_S1_regulon_heatmap.{ext}")

plt.close(fig)
