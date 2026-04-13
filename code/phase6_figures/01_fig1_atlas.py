"""
DAM-DRUG Figure 1 — Microglial Atlas & Substate Landscape
==========================================================
Panels:
  A. UMAP colored by microglial state
  B. Dot plot — marker genes × 6 states
  C. PAGA graph
  D. UMAP colored by diffusion pseudotime
  E. Violin: pseudotime distribution by state
  F. Stacked bar: substate composition by Braak stage

Runs on TRUBA (needs microglia_trajectory.h5ad ~3 GB).
Submit via 10_fig1_atlas.slurm, or run interactively:
  apptainer exec containers/scenic.sif python code/phase6_figures/01_fig1_atlas.py
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
from matplotlib.patheffects import withStroke
import seaborn as sns
from pathlib import Path

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
TRAJ    = PROJECT / "results/phase1/trajectory"
EXPLORE = PROJECT / "results/phase1/explore"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

H5AD = TRAJ / "microglia_trajectory.h5ad"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9,
    "axes.linewidth": 0.8, "pdf.fonttype": 42,
})

STATE_ORDER  = ["Homeostatic", "DAM", "DAM-IRM", "IRM", "LateAD-DAM", "LAM"]
# Okabe-Ito CVD-safe palette
STATE_COLORS = {"Homeostatic": "#0072B2", "DAM": "#E69F00", "DAM-IRM": "#F0E442",
                "IRM": "#56B4E9", "LateAD-DAM": "#D55E00", "LAM": "#CC79A7"}

MARKER_GENES = [
    # Homeostatic
    "P2RY12", "CX3CR1", "TMEM119",
    # DAM / disease
    "TREM2", "APOE", "SPP1", "LGALS3",
    # IRM
    "IFIT1", "MX1", "ISG15",
    # LAM
    "CD9", "GPNMB", "LPL",
    # TF targets
    "IKZF1", "IRF8", "SPI1", "BHLHE41", "PPARG",
]


def main():
    print(f"Loading {H5AD.name}...")
    mg = sc.read_h5ad(H5AD)
    print(f"  {mg.n_obs:,} cells × {mg.n_vars:,} genes")

    # Ensure state is categorical with correct order
    mg.obs["state"] = pd.Categorical(
        mg.obs["state"], categories=STATE_ORDER, ordered=True
    )
    mg.uns["state_colors"] = [STATE_COLORS[s] for s in STATE_ORDER]

    # ── Figure layout ──────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(20, 16))
    gs_outer = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.32)

    ax_a = fig.add_subplot(gs_outer[0, 0])
    ax_b = fig.add_subplot(gs_outer[0, 1])
    ax_c = fig.add_subplot(gs_outer[0, 2])
    ax_d = fig.add_subplot(gs_outer[1, 0])
    ax_e = fig.add_subplot(gs_outer[1, 1])
    ax_f = fig.add_subplot(gs_outer[1, 2])

    # ── Panel A: UMAP by state ─────────────────────────────────────────────────
    print("Panel A: UMAP by state...")
    sc.pl.umap(mg, color="state", ax=ax_a, show=False, frameon=False,
               palette=STATE_COLORS, legend_loc=None,
               title="", s=3, alpha=0.4)
    # Add thin semi-transparent border to each point collection to reduce
    # colour-overlap masking (OUP accessibility requirement)
    for coll in ax_a.collections:
        coll.set_linewidths(0.15)
        coll.set_edgecolors((0, 0, 0, 0.18))
        coll.set_rasterized(True)   # keep file size manageable for 236K points
    # Manual legend outside plot area
    handles = [mpatches.Patch(color=STATE_COLORS[s], label=s) for s in STATE_ORDER]
    ax_a.legend(handles=handles, fontsize=7, frameon=False,
                bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)
    ax_a.set_title("A  Microglial states (n=236,002)",
                   fontweight="bold", loc="left", fontsize=10)

    # ── Panel B: Dot plot ─────────────────────────────────────────────────────
    print("Panel B: Dot plot...")
    genes_present = [g for g in MARKER_GENES if g in mg.var_names]
    sc.pl.dotplot(mg, var_names=genes_present, groupby="state",
                  ax=ax_b, show=False, swap_axes=False,
                  categories_order=STATE_ORDER,
                  standard_scale="var", colorbar_title="Scaled expr",
                  size_title="% cells")
    ax_b.set_title("B  Marker gene expression",
                   fontweight="bold", loc="left", fontsize=10)
    ax_b.tick_params(axis="x", rotation=45)

    # ── Panel C: PAGA ─────────────────────────────────────────────────────────
    print("Panel C: PAGA graph...")
    if "paga" in mg.uns:
        # node_size_power=0 → all nodes identical size (n^0=1); scale=4 makes them large
        # color="state" explicit to guarantee correct palette mapping
        sc.pl.paga(mg, ax=ax_c, show=False, frameon=False,
                   color="state", fontsize=9,
                   node_size_scale=4.0, node_size_power=0,
                   edge_width_scale=1.2)
        # Add white halos to all node labels for contrast against edges
        for txt in ax_c.texts:
            txt.set_path_effects([withStroke(linewidth=2.5, foreground="white")])
        ax_c.set_title("C  PAGA connectivity graph",
                       fontweight="bold", loc="left", fontsize=10)
    else:
        ax_c.text(0.5, 0.5, "PAGA not in object\n(run sc.tl.paga first)",
                  ha="center", va="center", transform=ax_c.transAxes)

    # ── Panel D: UMAP by pseudotime ───────────────────────────────────────────
    print("Panel D: UMAP by pseudotime...")
    sc.pl.umap(mg, color="dpt_pseudotime", ax=ax_d, show=False, frameon=False,
               color_map="viridis", title="", s=3, alpha=0.6,
               colorbar_loc="right")
    for coll in ax_d.collections:
        coll.set_rasterized(True)
    ax_d.set_title("D  Diffusion pseudotime",
                   fontweight="bold", loc="left", fontsize=10)

    # ── Panel E: Violin pseudotime by state ───────────────────────────────────
    print("Panel E: Pseudotime violin...")
    pt_data = mg.obs[["state", "dpt_pseudotime"]].copy()
    sns.violinplot(
        data=pt_data, x="state", y="dpt_pseudotime",
        order=STATE_ORDER, hue="state", hue_order=STATE_ORDER,
        palette=STATE_COLORS, ax=ax_e, legend=False,
        inner="box", linewidth=0.8, density_norm="width", cut=0,
    )
    ax_e.set_xticks(range(len(STATE_ORDER)))
    ax_e.set_xticklabels(STATE_ORDER, rotation=30, ha="right", fontsize=8)
    ax_e.set_xlabel("")
    ax_e.set_ylabel("Diffusion pseudotime", fontsize=8)
    ax_e.set_ylim(0, 0.5)
    ax_e.set_title("E  Pseudotime distribution by state",
                   fontweight="bold", loc="left", fontsize=10)

    # ── Panel F: Braak composition stacked bar ────────────────────────────────
    print("Panel F: Braak composition...")
    braak_col = "Braak" if "Braak" in mg.obs.columns else next(
        (c for c in mg.obs.columns if "braak" in c.lower()), None
    )
    if braak_col:
        print(f"  Braak column: '{braak_col}' ({mg.obs[braak_col].dtype})")
        # Already string categories (e.g. "Braak I", "Braak V"); just fill NAs
        mg.obs["Braak_stage"] = mg.obs[braak_col].astype(object).fillna("Unknown")
        braak_col = "Braak_stage"
        braak_df = (
            mg.obs.groupby([braak_col, "state"], observed=True)
            .size()
            .reset_index(name="n")
        )
        braak_df["pct"] = braak_df.groupby(braak_col)["n"].transform(
            lambda x: x / x.sum() * 100
        )
        pivot = braak_df.pivot_table(
            index=braak_col, columns="state", values="pct", fill_value=0
        )
        # Sort Braak stages; exclude Unknown
        pivot = pivot[pivot.index != "Unknown"]
        stage_order = sorted(pivot.index, key=lambda s: (int("".join(filter(str.isdigit, s))) if any(c.isdigit() for c in s) else 99))
        pivot = pivot.reindex(stage_order)
        pivot = pivot[[s for s in STATE_ORDER if s in pivot.columns]]
        # Line plot per state — shows trajectories clearly for rare states too
        MARKERS = {"Homeostatic": "o", "DAM": "s", "DAM-IRM": "^",
                   "IRM": "D", "LateAD-DAM": "v", "LAM": "P"}
        x = range(len(pivot))
        for state in pivot.columns:
            light = state == "DAM-IRM"   # yellow needs visible edge
            ax_f.plot(x, pivot[state],
                      marker=MARKERS.get(state, "o"), markersize=5,
                      markeredgecolor="grey" if light else STATE_COLORS[state],
                      markeredgewidth=0.8,
                      color=STATE_COLORS[state], label=state, linewidth=1.5)
        ax_f.set_xticks(list(x))
        ax_f.set_xticklabels(pivot.index, rotation=30, ha="right", fontsize=8)
        ax_f.set_xlabel("Braak Stage", fontsize=8)
        ax_f.set_ylabel("% of cells", fontsize=8)
        ax_f.legend(title="State", fontsize=7, title_fontsize=7,
                    bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
    else:
        braak_candidates = [c for c in mg.obs.columns if "braak" in c.lower() or "stage" in c.lower()]
        msg = f"Braak column not found.\nobs columns with 'braak'/'stage':\n{braak_candidates}\nAll obs: {list(mg.obs.columns[:10])}..."
        print(f"  WARNING: {msg}")
        ax_f.text(0.5, 0.5, "Braak stage column\nnot found in obs",
                  ha="center", va="center", transform=ax_f.transAxes, fontsize=8, color="grey")
    ax_f.set_title("F  Substate composition by Braak stage",
                   fontweight="bold", loc="left", fontsize=10)

    # ── Save ──────────────────────────────────────────────────────────────────
    out_pdf = OUT / "fig1_atlas.pdf"
    out_png = OUT / "fig1_atlas.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    print(f"Saved:\n  {out_pdf}\n  {out_png}")


if __name__ == "__main__":
    main()
