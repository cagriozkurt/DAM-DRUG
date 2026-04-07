"""
DAM-DRUG Supplementary Figure S4 — Docking Validation ROC Curves
=================================================================
Enrichment/ROC analysis for each docking target, using known
actives (DUD-E / ChEMBL actives for the target) vs the
1,677-compound CNS approved-drug library as decoys.

Panels (one per dockable target, 3×3 or similar grid):
  PPARG, IRF8, BHLHE41, SPI1, MAF, RUNX1, IKZF1

Requires outputs from docking enrichment analysis on TRUBA
(results/phase4/roc/):
  {target}_roc.csv  — columns: fpr, tpr, auc, label

Run on TRUBA (after running docking enrichment step 4.x):
  apptainer exec --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \\
    ~/containers/scanpy-env.sif \\
    conda run -n scenic python code/phase6_figures/supp_fig_s4_docking_roc.py

NOTE: This script is TRUBA-dependent. The ROC data files are generated
by the docking enrichment analysis (roc/ subdirectory under phase4).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
ROC_DIR  = PROJECT / "results/phase4/roc"
OUT      = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})

# Okabe-Ito palette per target
TARGET_COLORS = {
    "PPARG":   "#56B4E9",
    "IRF8":    "#E69F00",
    "BHLHE41": "#009E73",
    "SPI1":    "#E91E63",
    "MAF":     "#F0E442",
    "RUNX1":   "#0072B2",
    "IKZF1":   "#D55E00",
}

TARGETS = list(TARGET_COLORS.keys())

# ── Plot ───────────────────────────────────────────────────────────────────
ncols = 3
nrows = (len(TARGETS) + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(11, nrows * 3.5),
                          constrained_layout=True)
axes_flat = axes.flatten()

for ax, target in zip(axes_flat, TARGETS):
    roc_file = ROC_DIR / f"{target}_roc.csv"
    if not roc_file.exists():
        ax.text(0.5, 0.5, f"ROC data not found\n{roc_file.name}",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color="#999999", style="italic")
        ax.set_title(target, fontsize=9, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        continue

    df  = pd.read_csv(roc_file)
    auc = df["auc"].iloc[0] if "auc" in df.columns else np.trapz(df["tpr"], df["fpr"])

    ax.plot(df["fpr"], df["tpr"], color=TARGET_COLORS[target],
            linewidth=2, label=f"AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], color="grey", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.fill_between(df["fpr"], df["tpr"], alpha=0.12, color=TARGET_COLORS[target])

    # EF1% annotation: enrichment factor at 1% FPR
    ef_idx = np.searchsorted(df["fpr"].values, 0.01, side="right")
    if ef_idx < len(df):
        tpr_at_1pct = df["tpr"].iloc[min(ef_idx, len(df) - 1)]
        ax.axvline(0.01, color="grey", linewidth=0.6, linestyle=":")
        ax.text(0.03, tpr_at_1pct, f"EF1%={tpr_at_1pct/0.01:.1f}×",
                fontsize=6.5, color="#333333", va="bottom")

    ax.set_xlabel("FPR (false-positive rate)", fontsize=7.5)
    ax.set_ylabel("TPR (recall)", fontsize=7.5)
    ax.set_title(target, fontsize=9, fontweight="bold")
    ax.legend(fontsize=7.5, loc="lower right", framealpha=0.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.spines[["top", "right"]].set_visible(False)

# Hide unused axes
for ax in axes_flat[len(TARGETS):]:
    ax.set_visible(False)

fig.suptitle(
    "Supplementary Figure S4 — Docking Enrichment ROC Curves (Known Actives vs CNS Library)",
    fontsize=10
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"supp_fig_S4_docking_roc.{ext}",
                bbox_inches="tight", dpi=300)
    print(f"Saved → {OUT}/supp_fig_S4_docking_roc.{ext}")

plt.close(fig)
