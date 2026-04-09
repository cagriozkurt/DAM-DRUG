"""
DAM-DRUG Figure 2 — Target Prioritization
==========================================
Panels:
  A. Pseudobulk DGE bubble: 9 TFs × 4 microglial states
  B. AUCell regulon heatmap: regulons × 6 states
  C. Pseudotime–regulon correlation waterfall
  D. GSE95587 bulk replication forest plot

Run locally:
  python code/phase6_figures/fig2_targets.py

Requires: pandas, numpy, matplotlib, seaborn, scipy
Aggregated GRN files (B, C): rsync from TRUBA first:
  rsync truba:/arf/scratch/mozkurt/DAM-DRUG/results/phase2/GRN/regulon_auc_by_state_aggregated.csv results/phase2/GRN/
  rsync truba:/arf/scratch/mozkurt/DAM-DRUG/results/phase2/GRN/regulon_pseudotime_corr.csv results/phase2/GRN/
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
GRN     = PROJECT / "results/phase2/GRN"
DGE_PB  = PROJECT / "results/phase1/DGE_pseudobulk"
GSE     = PROJECT / "results/phase5/gse95587"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42,   # editable text in Illustrator
})
STATE_ORDER  = ["Homeostatic", "DAM", "DAM-IRM", "IRM", "LateAD-DAM", "LAM"]
# Okabe-Ito CVD-safe palette (matches fig1)
STATE_COLORS = {"Homeostatic": "#0072B2", "DAM": "#E69F00", "DAM-IRM": "#F0E442",
                "IRM": "#56B4E9", "LateAD-DAM": "#D55E00", "LAM": "#CC79A7"}
TF_ORDER = ["BHLHE41", "IRF8", "IKZF1", "RUNX1", "CEBPB", "RELB",
            "PPARG", "SPI1", "BHLHE40"]

# ── Panel A: Pseudobulk DGE bubble ────────────────────────────────────────────
def panel_dge_bubble(ax):
    df = pd.read_csv(DGE_PB / "TF_pseudobulk_summary.csv")
    # Keep only main state contrasts
    contrasts_order = ["DAM_vs_Homeostatic", "IRM_vs_Homeostatic", "LAM_vs_Homeostatic"]
    contrast_labels = {"DAM_vs_Homeostatic": "DAM", "IRM_vs_Homeostatic": "IRM",
                       "LAM_vs_Homeostatic": "LAM"}
    df = df[df["contrast"].isin(contrasts_order)].copy()
    df["contrast_label"] = df["contrast"].map(contrast_labels)
    df["sig"] = df["padj"] < 0.05
    df["logFC_clipped"] = df["log2FoldChange"].clip(-4, 4)
    df["size"] = (df["log2FoldChange"].abs().clip(0, 4) * 30 + 20)
    df["alpha"] = np.where(df["sig"], 0.9, 0.25)

    x_map = {c: i for i, c in enumerate(["DAM", "IRM", "LAM"])}
    y_map = {tf: i for i, tf in enumerate(reversed(TF_ORDER))}

    cmap = plt.cm.RdBu_r
    norm = matplotlib.colors.TwoSlopeNorm(vcenter=0, vmin=-4, vmax=4)

    for _, row in df.iterrows():
        x = x_map.get(row["contrast_label"], None)
        y = y_map.get(row["names"], None)
        if x is None or y is None:
            continue
        color = cmap(norm(row["logFC_clipped"]))
        ax.scatter(x, y, s=row["size"], color=color,
                   alpha=row["alpha"], linewidths=0.5,
                   edgecolors="k" if row["sig"] else "none", zorder=3)

    ax.set_xticks(range(3))
    ax.set_xticklabels(["DAM", "IRM", "LAM"], fontsize=9)
    ax.set_yticks(range(len(TF_ORDER)))
    ax.set_yticklabels(list(reversed(TF_ORDER)), fontsize=9)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(-0.6, len(TF_ORDER) - 0.4)
    ax.grid(axis="both", color="lightgrey", linewidth=0.5, zorder=0)
    ax.set_title("A  Pseudobulk DGE", fontweight="bold", loc="left", fontsize=10)
    ax.set_xlabel("Microglial state vs Homeostatic")

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.02, aspect=15)
    cb.set_label("log₂FC", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # Size legend — placed below the axes to avoid overlapping the colorbar
    for lfc, label in [(1, "1"), (2, "2"), (3, "3")]:
        ax.scatter([], [], s=lfc * 30 + 20, color="grey", alpha=0.6, label=f"|log₂FC|={label}")
    ax.legend(fontsize=7, frameon=False, title="|log₂FC|", title_fontsize=7,
              bbox_to_anchor=(0.0, -0.14), loc="upper left", borderaxespad=0, ncol=3)


# ── Panel B: AUCell regulon heatmap ───────────────────────────────────────────
def panel_auc_heatmap(ax):
    # Prefer aggregated file, fall back to seed=42
    auc_file = GRN / "regulon_auc_by_state_aggregated.csv"
    if not auc_file.exists():
        auc_file = GRN / "regulon_auc_by_state.csv"
        print(f"  [Panel B] Using seed=42 AUC data. Rsync aggregated file for final figure.")

    auc = pd.read_csv(auc_file, index_col="state")

    # Reorder states
    states_present = [s for s in STATE_ORDER if s in auc.index]
    auc = auc.loc[states_present]

    # Select top 20 regulons by range across states (most informative)
    reg_range = auc.max() - auc.min()
    top20 = reg_range.nlargest(20).index.tolist()
    auc = auc[top20]

    # Sort by DAM AUC descending
    if "DAM" in auc.index:
        reg_order = auc.loc["DAM"].sort_values(ascending=False).index.tolist()
    else:
        reg_order = auc.columns.tolist()
    auc = auc[reg_order]

    # Normalize each regulon 0–1 for display
    auc_norm = (auc - auc.min()) / (auc.max() - auc.min() + 1e-9)

    sns.heatmap(
        auc_norm.T, ax=ax,
        cmap="YlOrRd", vmin=0, vmax=1,
        linewidths=0, linecolor="white",  # no cell borders — eliminates moiré
        cbar_kws={"label": "Scaled AUC", "shrink": 0.6, "pad": 0.01},
        yticklabels=True, xticklabels=True,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for coll in ax.collections:
        coll.set_rasterized(True)  # rasterize mesh; keep axis/text vector

    # Highlight target TFs in red
    target_tfs = {"SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
                  "RELB", "BHLHE40", "BHLHE41"}
    for label in ax.get_yticklabels():
        tf = label.get_text().split("(")[0]
        if tf in target_tfs:
            label.set_color("#D55E00")
            label.set_fontweight("bold")

    ax.set_title("B  Regulon AUCell activity", fontweight="bold", loc="left", fontsize=10)


# ── Panel C: Pseudotime correlation waterfall ─────────────────────────────────
def panel_pseudotime_corr(ax):
    corr_file = GRN / "regulon_pseudotime_corr.csv"
    if not corr_file.exists():
        ax.text(0.5, 0.5, "regulon_pseudotime_corr.csv\nnot found\n(rsync from TRUBA)",
                ha="center", va="center", transform=ax.transAxes, fontsize=9, color="grey")
        ax.set_title("C  Pseudotime–regulon correlation", fontweight="bold", loc="left", fontsize=10)
        return

    df = pd.read_csv(corr_file)
    # Keep top 20 by |rho| — matches the heatmap panel selection
    df = df.reindex(df["rho"].abs().nlargest(20).index).sort_values("rho")
    colors = ["#0072B2" if r < 0 else "#D55E00" for r in df["rho"]]
    edge   = ["#000000" if abs(r) >= 0.3 else "none" for r in df["rho"]]

    ax.barh(range(len(df)), df["rho"], color=colors, edgecolor=edge,
            linewidth=0.8, height=0.7)
    ax.axvline(0, color="k", linewidth=0.8)
    ax.axvline(0.3,  color="#D55E00", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.axvline(-0.3, color="#0072B2", linewidth=0.8, linestyle="--", alpha=0.6)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["regulon"], fontsize=6.5)
    ax.set_xlabel("Spearman ρ (AUC vs pseudotime)", fontsize=8)

    # Legend for bar colors
    pos_patch = mpatches.Patch(color="#D55E00", label="ρ > 0 (disease-up)")
    neg_patch = mpatches.Patch(color="#0072B2", label="ρ < 0 (disease-down)")
    ax.legend(handles=[pos_patch, neg_patch], fontsize=7, frameon=False,
              loc="lower right")

    # Bold target TF labels
    target_tfs = {"SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
                  "RELB", "BHLHE40", "BHLHE41"}
    for label in ax.get_yticklabels():
        tf = label.get_text().split("(")[0]
        if tf in target_tfs:
            label.set_color("#D55E00")
            label.set_fontweight("bold")

    ax.set_title("C  Pseudotime–regulon correlation", fontweight="bold", loc="left", fontsize=10)


# ── Panel D: GSE95587 forest plot ────────────────────────────────────────────
def panel_gse_forest(ax):
    df = pd.read_csv(GSE / "deseq2_target_tfs.csv")
    df = df.sort_values("log2FC", ascending=False).reset_index(drop=True)

    colors = ["#D55E00" if p < 0.05 else "grey" for p in df["padj"]]
    y = range(len(df))

    ax.barh(y, df["log2FC"], xerr=df["lfcSE"],
            color=colors, error_kw={"linewidth": 1, "capsize": 3},
            height=0.6, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="k", linewidth=0.8)
    ax.axvline(0.5,  color="#D55E00", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.axvline(-0.5, color="#0072B2", linewidth=0.8, linestyle="--", alpha=0.5)

    ax.set_yticks(list(y))
    ax.set_yticklabels(df["gene"], fontsize=9)
    ax.set_xlabel("log₂FC (AD vs Control)", fontsize=8)

    # padj annotations
    for i, (_, row) in enumerate(df.iterrows()):
        if row["padj"] < 0.001:
            stars = "***"
        elif row["padj"] < 0.01:
            stars = "**"
        elif row["padj"] < 0.05:
            stars = "*"
        else:
            stars = "ns"
        # Anchor to right edge of CI whisker; use transform offset for consistent gap
        x_pos = row["log2FC"] + row["lfcSE"]
        ax.annotate(stars, xy=(x_pos, i), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=8,
                    color="#D55E00" if stars != "ns" else "grey")

    sig_patch = mpatches.Patch(color="#D55E00", label="padj < 0.05")
    ns_patch  = mpatches.Patch(color="grey",    label="ns")
    ax.legend(handles=[sig_patch, ns_patch], fontsize=7, frameon=False, loc="upper left")
    ax.set_title("D  GSE95587 bulk replication\n(fusiform gyrus, n=117)",
                 fontweight="bold", loc="left", fontsize=10)


# ── Compose figure ────────────────────────────────────────────────────────────
def main():
    # 3-panel layout: A (top-left), C (bottom-left), B (right, full height)
    # Panel D (GSE95587 forest plot) moved to Figure 4A to avoid duplication.
    fig = plt.figure(figsize=(18, 10))
    gs  = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.35,
                            width_ratios=[1, 1.8])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[:, 1])   # spans both rows
    ax_c = fig.add_subplot(gs[1, 0])

    print("Panel A: DGE bubble...")
    panel_dge_bubble(ax_a)
    print("Panel B: AUCell heatmap...")
    panel_auc_heatmap(ax_b)
    print("Panel C: Pseudotime correlation...")
    panel_pseudotime_corr(ax_c)

    out_pdf = OUT / "fig2_target_prioritization.pdf"
    out_png = OUT / "fig2_target_prioritization.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    print(f"Saved:\n  {out_pdf}\n  {out_png}")


if __name__ == "__main__":
    main()
