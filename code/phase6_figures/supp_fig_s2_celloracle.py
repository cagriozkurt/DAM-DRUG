"""
DAM-DRUG Supplementary Figure S2 — CellOracle TF Perturbation Summary
=======================================================================
Panels:
  A–C. Bar charts: mean UMAP displacement magnitude per microglial state
       for IKZF1 KO, IRF8 KO, and SPI1 KO.
  D.   Side-by-side comparison: LateAD-DAM magnitude across three KOs.

Note: Quiver/vector-field PDFs (KO_*_quiver.pdf, KO_*_grid.pdf) are
pre-generated on TRUBA and available at:
  results/phase5/celloracle/KO_IKZF1_quiver.pdf
  results/phase5/celloracle/KO_IRF8_quiver.pdf
  results/phase5/celloracle/KO_SPI1_quiver.pdf

Runs locally:
  python code/phase6_figures/supp_fig_s2_celloracle.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

PROJECT = Path("/Volumes/PortableSSD/untitled folder/DAM-DRUG")
CO_DIR  = PROJECT / "results/phase5/celloracle"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})

STATE_ORDER  = ["Homeostatic", "DAM", "DAM_IRM", "IRM", "LAM", "LateAD_DAM"]
STATE_LABELS = {"Homeostatic": "Homeostatic", "DAM": "DAM", "DAM_IRM": "DAM-IRM",
                "IRM": "IRM", "LAM": "LAM", "LateAD_DAM": "LateAD-DAM"}
STATE_COLORS = {"Homeostatic": "#0072B2", "DAM": "#E69F00", "DAM_IRM": "#F0E442",
                "IRM": "#56B4E9", "LAM": "#CC79A7", "LateAD_DAM": "#D55E00"}

KOS = [
    ("IKZF1", "#D55E00"),   # DAM-DRUG primary target
    ("IRF8",  "#E69F00"),
    ("SPI1",  "#56B4E9"),
]

# ── Load delta summaries ───────────────────────────────────────────────────
def load_ko(tf):
    df = pd.read_csv(CO_DIR / f"KO_{tf}_delta_summary.csv")
    # ensure all states present (fill 0 if missing)
    df = df.set_index("state").reindex(STATE_ORDER, fill_value=0).reset_index()
    return df

data = {tf: load_ko(tf) for tf, _ in KOS}

# ── Figure layout ─────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 8))
gs  = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.35)

axes_ko   = [fig.add_subplot(gs[0, i]) for i in range(3)]
ax_comp   = fig.add_subplot(gs[1, :])

# ── Panels A–C: per-state bar chart per KO ────────────────────────────────
for ax, (tf, color) in zip(axes_ko, KOS):
    df     = data[tf]
    states = [STATE_LABELS[s] for s in STATE_ORDER]
    mags   = [df.loc[df["state"] == s, "magnitude"].values[0]
              if s in df["state"].values else 0 for s in STATE_ORDER]
    bar_colors = [STATE_COLORS[s] for s in STATE_ORDER]

    bars = ax.bar(states, mags, color=bar_colors, edgecolor="white", linewidth=0.4)

    # Highlight LateAD_DAM bar with red edge
    for bar, state in zip(bars, STATE_ORDER):
        if state == "LateAD_DAM":
            bar.set_edgecolor("#D55E00")
            bar.set_linewidth(1.5)

    ax.set_title(f"{tf} KO", fontsize=10, fontweight="bold", pad=4)
    ax.set_ylabel("Mean UMAP displacement\n(magnitude)", fontsize=7.5)
    ax.set_xticks(range(len(states)))
    ax.set_xticklabels(states, rotation=40, ha="right", fontsize=7.5)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
    ax.spines[["top", "right"]].set_visible(False)

    # Panel label
    panel_labels = ["A", "B", "C"]
    ax.text(-0.15, 1.08, panel_labels[KOS.index((tf, color))],
            transform=ax.transAxes, fontsize=12, fontweight="bold")

# ── Panel D: LateAD-DAM comparison across three KOs ──────────────────────
ko_names  = [tf for tf, _ in KOS]
ko_colors = [c for _, c in KOS]
lateAD_mags = []
for tf, _ in KOS:
    df = data[tf]
    val = df.loc[df["state"] == "LateAD_DAM", "magnitude"]
    lateAD_mags.append(val.values[0] if len(val) > 0 else 0)

bars = ax_comp.bar(ko_names, lateAD_mags, color=ko_colors,
                   edgecolor="white", linewidth=0.4, width=0.5)
ax_comp.set_title("Panel D — LateAD-DAM displacement magnitude by KO",
                   fontsize=9, pad=4)
ax_comp.set_ylabel("Mean UMAP displacement\n(magnitude, LateAD-DAM)", fontsize=8)
ax_comp.set_xticks(range(len(ko_names)))
ax_comp.set_xticklabels([f"{n} KO" for n in ko_names], fontsize=9)
ax_comp.set_ylim(0, max(lateAD_mags) * 1.25)
ax_comp.spines[["top", "right"]].set_visible(False)
ax_comp.text(-0.06, 1.05, "D", transform=ax_comp.transAxes,
             fontsize=12, fontweight="bold")

# Annotate values
for bar, val in zip(bars, lateAD_mags):
    ax_comp.text(bar.get_x() + bar.get_width() / 2,
                 val + max(lateAD_mags) * 0.02,
                 f"{val:.3f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

fig.suptitle(
    "Supplementary Figure S2 — CellOracle TF KO Perturbation Magnitudes",
    fontsize=10, y=1.02
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"supp_fig_S2_celloracle.{ext}",
                bbox_inches="tight", dpi=300)
    print(f"Saved → {OUT}/supp_fig_S2_celloracle.{ext}")

plt.close(fig)
print("\nNote: Full quiver/vector-field PDFs available at:")
print("  results/phase5/celloracle/KO_IKZF1_quiver.pdf")
print("  results/phase5/celloracle/KO_IRF8_quiver.pdf")
print("  results/phase5/celloracle/KO_SPI1_quiver.pdf")
