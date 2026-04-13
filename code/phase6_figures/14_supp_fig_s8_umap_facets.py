"""
DAM-DRUG Supplementary Figure S8 — UMAP Faceted by Microglial Substate
=======================================================================
Small-multiple UMAPs (one per substate) as accessibility companion to
fig1 Panel A.  Removes colour-overlap masking for dense overlapping states.

Requires: microglia_trajectory.h5ad (TRUBA)

Run on TRUBA:
  apptainer exec containers/scenic.sif \
    python code/phase6_figures/14_supp_fig_s8_umap_facets.py
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
TRAJ    = PROJECT / "results/phase1/trajectory"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})

STATE_ORDER  = ["Homeostatic", "DAM", "DAM-IRM", "IRM", "LateAD-DAM", "LAM"]
STATE_COLORS = {"Homeostatic": "#0072B2", "DAM": "#E69F00", "DAM-IRM": "#F0E442",
                "IRM": "#56B4E9", "LateAD-DAM": "#D55E00", "LAM": "#CC79A7"}


def main():
    print(f"Loading microglia_trajectory.h5ad ...")
    mg = sc.read_h5ad(TRAJ / "microglia_trajectory.h5ad")
    print(f"  {mg.n_obs:,} cells × {mg.n_vars:,} genes")

    coords = mg.obsm["X_umap"]
    states = mg.obs["state"].astype(str)

    n_states = len(STATE_ORDER)
    ncols = 3
    nrows = (n_states + ncols - 1) // ncols   # = 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 8),
                             gridspec_kw={"hspace": 0.38, "wspace": 0.18})
    axes = axes.flat

    for ax, state in zip(axes, STATE_ORDER):
        # Background: all cells in light grey
        ax.scatter(coords[:, 0], coords[:, 1],
                   s=0.15, alpha=0.15, c="#CCCCCC", rasterized=True,
                   linewidths=0)
        # Foreground: highlighted state
        mask = states == state
        n = mask.sum()
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   s=0.8, alpha=0.7,
                   c=STATE_COLORS.get(state, "grey"),
                   rasterized=True,
                   linewidths=0.15, edgecolors=(0, 0, 0, 0.20))
        ax.set_title(f"{state}\n(n={n:,})", fontsize=8, fontweight="bold",
                     color=STATE_COLORS.get(state, "black"), pad=3)
        ax.axis("off")

    # Hide any spare axes (none here for 6 states / 3 cols)
    for ax in list(axes)[n_states:]:
        ax.set_visible(False)

    fig.suptitle(
        "Supplementary Figure S8 — Microglial substate UMAPs (faceted)\n"
        "Grey points = all other nuclei; coloured = highlighted state",
        fontsize=9, y=1.01
    )

    for ext in ("pdf", "png"):
        out = OUT / f"supp_fig_S8_umap_facets.{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=300)
        print(f"Saved → {out}")

    plt.close(fig)


if __name__ == "__main__":
    main()
