"""
DAM-DRUG Figure — CellChat L-R Interaction Results
====================================================
Panel A: Chord diagram (imported from chord_count.pdf)
Panel B: Top L-R pairs → Microglia dot/bar plot (from lr_summary.csv)
Panel C: Pathway-level summary bar chart (sum_prob, Microglia as receiver)

Output: results/figures/fig_cellchat.pdf  +  .png

Runs locally from synced CSVs and PDFs.
"""

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT  = Path(os.environ.get("DAM_DRUG_DIR",
                str(Path(__file__).resolve().parents[2])))
LR_DIR   = PROJECT / "results/phase2/LR/cellchat"
OUT_DIR  = PROJECT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CHORD_PDF = LR_DIR / "chord_count.pdf"
LR_CSV    = LR_DIR / "lr_summary.csv"

# ── CVD-safe palette (Okabe-Ito) for cell types ───────────────────────────────
CELL_COLORS = {
    "Microglia":    "#0072B2",
    "Astrocyte":    "#E69F00",
    "ExcNeuron":    "#56B4E9",
    "Pvalb":        "#009E73",
    "Sst":          "#F0E442",
    "Chandelier":   "#D55E00",
    "Lamp5":        "#CC79A7",
    "Lamp5 Lhx6":   "#999999",
    "OPC":          "#44AA99",
    "Oligodendrocyte": "#882255",
    "Vascular":     "#DDCC77",
    "Vip":          "#117733",
    "Sncg":         "#AA4499",
    "Pax6":         "#332288",
    "Sst Chodl":    "#6699CC",
}
DEFAULT_COLOR = "#888888"

# Marker shapes per sender — paired with colours for dual encoding
CELL_SHAPES = {
    "Microglia":       "o",
    "Astrocyte":       "s",
    "ExcNeuron":       "^",
    "Pvalb":           "D",
    "Sst":             "v",
    "Chandelier":      "P",
    "Lamp5":           "X",
    "Lamp5 Lhx6":      "h",
    "OPC":             "*",
    "Oligodendrocyte": "p",
    "Vascular":        "H",
    "Vip":             "<",
    "Sncg":            ">",
    "Pax6":            "8",
    "Sst Chodl":       "+",
}
DEFAULT_SHAPE = "o"

# Pathway colors (Okabe-Ito extended)
PATHWAY_COLORS = {
    "SLIT":       "#D55E00",
    "GAS":        "#0072B2",
    "TGFb":       "#E69F00",
    "SPP1":       "#009E73",
    "COMPLEMENT": "#CC79A7",
    "NRG":        "#56B4E9",
    "BMP":        "#F0E442",
    "TULP":       "#999999",
    "SLITRK":     "#44AA99",
    "CX3C":       "#882255",
}
DEFAULT_PATH_COLOR = "#BBBBBB"

def main():
    # ── Load data ─────────────────────────────────────────────────────────────
    print(f"Loading {LR_CSV.name} ...")
    df = pd.read_csv(LR_CSV)
    print(f"  {len(df):,} L-R pairs, {df['source'].nunique()} cell types")

    # Microglia as receiver
    recv = df[df["target"] == "Microglia"].copy()
    recv = recv.sort_values("prob", ascending=False)
    print(f"  {len(recv):,} pairs with Microglia as receiver")

    # Top 20 individual pairs
    top_pairs = recv.head(20).copy()
    top_pairs["label"] = (top_pairs["source"] + "\n→ " +
                          top_pairs["ligand"] + " : " + top_pairs["receptor"])

    # Pathway-level summary (sum prob, Microglia receiver)
    pathway_sum = (recv.groupby("pathway_name")["prob"]
                   .sum()
                   .sort_values(ascending=False)
                   .head(10)
                   .reset_index())
    pathway_sum.columns = ["pathway", "sum_prob"]

    # ── Convert chord PDF to image ────────────────────────────────────────────
    chord_img = None
    if CHORD_PDF.exists():
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                out_stem = str(Path(tmpdir) / "chord")
                subprocess.run(
                    ["pdftoppm", "-r", "200", "-png", "-singlefile",
                     str(CHORD_PDF), out_stem],
                    check=True, capture_output=True
                )
                png_path = Path(out_stem + ".png")
                if png_path.exists():
                    chord_img = mpimg.imread(str(png_path))
                    print(f"  Chord diagram loaded: {chord_img.shape}")
        except Exception as e:
            print(f"  WARNING: chord PDF import failed ({e}); Panel A will be blank")
    else:
        print(f"  WARNING: {CHORD_PDF} not found; Panel A will be blank")

    # ── Layout ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 10))
    gs  = GridSpec(2, 2, figure=fig,
                   left=0.07, right=0.97, top=0.93, bottom=0.08,
                   hspace=0.40, wspace=0.38,
                   height_ratios=[1.1, 1],
                   width_ratios=[1, 1.3])

    ax_chord   = fig.add_subplot(gs[0, 0])   # Panel A — chord diagram
    ax_pathway = fig.add_subplot(gs[0, 1])   # Panel C — pathway bar
    ax_pairs   = fig.add_subplot(gs[1, :])   # Panel B — top L-R pairs (full width)

    # ── Panel A: Chord diagram ────────────────────────────────────────────────
    ax_chord.set_facecolor("white")
    if chord_img is not None:
        ax_chord.imshow(chord_img)
        ax_chord.axis("off")
    else:
        ax_chord.text(0.5, 0.5, "chord_count.pdf\n(not found)",
                      ha="center", va="center", transform=ax_chord.transAxes,
                      fontsize=10, color="gray")
        ax_chord.set_xticks([]); ax_chord.set_yticks([])
    ax_chord.set_title("A  Cell–cell interaction network\n(number of L-R pairs)",
                       fontsize=10, fontweight="bold", loc="left", pad=4)

    # ── Panel C: Pathway bar chart ────────────────────────────────────────────
    colors_p = [PATHWAY_COLORS.get(p, DEFAULT_PATH_COLOR) for p in pathway_sum["pathway"]]
    bars = ax_pathway.barh(pathway_sum["pathway"][::-1],
                           pathway_sum["sum_prob"][::-1],
                           color=colors_p[::-1], edgecolor="white", linewidth=0.5)
    ax_pathway.set_xlabel("Cumulative interaction probability", fontsize=9)
    ax_pathway.set_title("C  Top signalling pathways\n→ Microglia (receiver)",
                         fontsize=10, fontweight="bold", loc="left", pad=4)
    ax_pathway.spines[["top", "right"]].set_visible(False)
    ax_pathway.tick_params(axis="y", labelsize=8)
    ax_pathway.tick_params(axis="x", labelsize=8)
    # Value labels
    for bar in bars:
        w = bar.get_width()
        ax_pathway.text(w + 0.003, bar.get_y() + bar.get_height() / 2,
                        f"{w:.3f}", va="center", fontsize=7)

    # ── Panel B: Top 20 L-R pairs — horizontal dot/bar plot ──────────────────
    y_pos    = np.arange(len(top_pairs))
    colors_b = [CELL_COLORS.get(s, DEFAULT_COLOR) for s in top_pairs["source"]]
    sizes    = (top_pairs["prob"] / top_pairs["prob"].max()) * 200

    # Horizontal bars (prob) — light fill behind dots
    ax_pairs.barh(y_pos, top_pairs["prob"].values,
                  color=colors_b, alpha=0.25, height=0.55)

    # Dots: plot per sender so each gets its own shape + colour (dual encoding)
    seen_sources = list(top_pairs["source"].unique())
    scatter_handles = []
    for src in seen_sources:
        mask = top_pairs["source"] == src
        idx  = top_pairs.index[mask]
        rows = top_pairs.loc[idx]
        yy   = np.where(top_pairs["source"] == src)[0]
        clr  = CELL_COLORS.get(src, DEFAULT_COLOR)
        mrk  = CELL_SHAPES.get(src, DEFAULT_SHAPE)
        sz   = (rows["prob"] / top_pairs["prob"].max()) * 200
        ax_pairs.scatter(rows["prob"].values, yy,
                         s=sz.values, c=clr, marker=mrk,
                         zorder=3, edgecolors="white", linewidths=0.4)
        scatter_handles.append(
            plt.scatter([], [], s=60, c=clr, marker=mrk, label=src,
                        edgecolors="white", linewidths=0.4)
        )

    # Y-tick labels
    pair_labels = (top_pairs["source"] + "  →  " +
                   top_pairs["ligand"] + " : " + top_pairs["receptor"] +
                   "  [" + top_pairs["pathway_name"] + "]")
    ax_pairs.set_yticks(y_pos)
    ax_pairs.set_yticklabels(pair_labels.values, fontsize=7.5)
    ax_pairs.invert_yaxis()
    ax_pairs.set_xlabel("Interaction probability", fontsize=9)
    ax_pairs.set_title("B  Top L-R pairs — Microglia as receiver",
                       fontsize=10, fontweight="bold", loc="left", pad=4)
    ax_pairs.spines[["top", "right"]].set_visible(False)
    ax_pairs.tick_params(axis="x", labelsize=8)

    # Legend: colour + shape per sender (dual encoding)
    ax_pairs.legend(handles=scatter_handles, title="Sender cell type",
                    fontsize=7, title_fontsize=7.5,
                    loc="lower right", framealpha=0.85,
                    ncol=2)

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ("pdf", "png"):
        out_path = OUT_DIR / f"fig_cellchat.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved {out_path}")

    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
