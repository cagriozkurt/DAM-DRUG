"""
DAM-DRUG Figure 4 — Multi-Modal IKZF1 Validation
==================================================
Panels:
  A. GSE95587 forest plot  — bulk AD replication (fusiform gyrus, n=117)
  B. CellOracle KO delta   — perturbation magnitude per state (IKZF1/IRF8/SPI1)
  C. ATAC motif enrichment — fold-enrichment in AD-upregulated chromatin peaks
  D. IKZF1 convergent evidence summary — 6 modalities, lollipop

Run locally (all required CSVs present):
  python code/phase6_figures/fig4_validation.py

Panel C: loads results/phase2/atac_da/fimo_enrichment_summary.csv if present;
          falls back to hardcoded values from the TRUBA run.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
GSE     = PROJECT / "results/phase5/gse95587"
CO      = PROJECT / "results/phase5/celloracle"
ATAC    = PROJECT / "results/phase2/atac_da"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9,
    "axes.linewidth": 0.8, "pdf.fonttype": 42,
})

# Okabe-Ito — matches fig1/fig2
STATE_ORDER  = ["Homeostatic", "DAM", "DAM_IRM", "IRM", "LateAD_DAM", "LAM"]
STATE_LABELS = {"Homeostatic": "Homeostatic", "DAM": "DAM", "DAM_IRM": "DAM-IRM",
                "IRM": "IRM", "LateAD_DAM": "LateAD-DAM", "LAM": "LAM"}
STATE_COLORS = {"Homeostatic": "#0072B2", "DAM": "#E69F00", "DAM_IRM": "#F0E442",
                "IRM": "#56B4E9", "LateAD_DAM": "#D55E00", "LAM": "#CC79A7"}
KO_COLORS    = {"IKZF1": "#D55E00", "IRF8": "#0072B2", "SPI1": "#009E73"}


# ── Panel A: GSE95587 Forest Plot ─────────────────────────────────────────────
def panel_gse_forest(ax):
    df = pd.read_csv(GSE / "deseq2_target_tfs.csv")
    df = df.sort_values("log2FC", ascending=False).reset_index(drop=True)

    colors = ["#D55E00" if p < 0.05 else "#AAAAAA" for p in df["padj"]]
    y = np.arange(len(df))

    ax.barh(y, df["log2FC"], xerr=df["lfcSE"],
            color=colors, error_kw={"linewidth": 0.9, "capsize": 3},
            height=0.6, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="k", linewidth=0.8)
    ax.axvline(0.5,  color="#D55E00", linewidth=0.7, linestyle="--", alpha=0.5)
    ax.axvline(-0.5, color="#0072B2", linewidth=0.7, linestyle="--", alpha=0.5)

    ax.set_yticks(list(y))
    ax.set_yticklabels(df["gene"], fontsize=9)
    ax.set_xlabel("log₂FC (AD vs Control)", fontsize=8)

    # Significance annotations right of CI whisker
    for i, (_, row) in enumerate(df.iterrows()):
        if row["padj"] < 0.001:   stars = "***"
        elif row["padj"] < 0.01:  stars = "**"
        elif row["padj"] < 0.05:  stars = "*"
        else:                     stars = "ns"
        ax.annotate(stars, xy=(row["log2FC"] + row["lfcSE"], i),
                    xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=8,
                    color="#D55E00" if stars != "ns" else "#888888")

    sig_patch = mpatches.Patch(color="#D55E00", label="padj < 0.05")
    ns_patch  = mpatches.Patch(color="#AAAAAA", label="ns")
    ax.legend(handles=[sig_patch, ns_patch], fontsize=7, frameon=False,
              loc="upper left")
    ax.set_title("A  GSE95587 bulk replication\n(fusiform gyrus, n=117)",
                 fontweight="bold", loc="left", fontsize=10)


# ── Panel B: CellOracle KO delta magnitude ────────────────────────────────────
def panel_celloracle(ax):
    dfs = {}
    for tf in ["IKZF1", "IRF8", "SPI1"]:
        f = CO / f"KO_{tf}_delta_summary.csv"
        d = pd.read_csv(f).set_index("state")
        dfs[tf] = d["magnitude"]

    states = [s for s in STATE_ORDER if s in dfs["IKZF1"].index]
    x      = np.arange(len(states))
    width  = 0.24
    offsets = [-width, 0, width]

    for i, tf in enumerate(["IKZF1", "IRF8", "SPI1"]):
        vals = [dfs[tf].get(s, 0) for s in states]
        bars = ax.bar(x + offsets[i], vals, width=width * 0.92,
                      color=KO_COLORS[tf], label=tf,
                      edgecolor="white", linewidth=0.4)

    # Highlight LateAD_DAM column
    if "LateAD_DAM" in states:
        lx = states.index("LateAD_DAM")
        ax.axvspan(lx - width * 1.8, lx + width * 1.8,
                   color="#FFE0CC", alpha=0.35, zorder=0)

    ax.set_xticks(x)
    ax.set_xticklabels([STATE_LABELS[s] for s in states],
                       rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Mean |ΔUMAP| per cell", fontsize=8)
    ax.set_title("B  CellOracle in-silico KO perturbation\n(UMAP shift magnitude)",
                 fontweight="bold", loc="left", fontsize=10)
    ax.legend(fontsize=7, frameon=False, title="KO", title_fontsize=7,
              bbox_to_anchor=(1.01, 1), loc="upper left")

    # Annotate peak value (LateAD_DAM IKZF1); ylim headroom prevents spine collision
    if "LateAD_DAM" in states:
        lx = states.index("LateAD_DAM")
        peak = dfs["IKZF1"].get("LateAD_DAM", 0)
        ax.text(lx + offsets[0], peak + 0.001, f"{peak:.3f}",
                ha="center", va="bottom", fontsize=6.5,
                color=KO_COLORS["IKZF1"], fontweight="bold")
    ax.set_ylim(bottom=0, top=0.055)


# ── Panel C: ATAC motif enrichment ────────────────────────────────────────────
def panel_atac(ax):
    # Try to load from file (rsynced from TRUBA); fall back to known results
    fimo_file = ATAC / "fimo_enrichment_summary.csv"
    if fimo_file.exists():
        df = pd.read_csv(fimo_file)
    else:
        print("  [Panel C] fimo_enrichment_summary.csv not found — using hardcoded results")
        df = pd.DataFrame({
            "tf":          ["PPARG", "IKZF1", "RELB", "SPI1", "IRF8", "RUNX1", "CEBPB"],
            "fold_enrich": [4.5,      3.9,     3.9,    3.8,    3.5,    3.5,     3.3],
            "n_hits":      [8,        7,       6,      7,      6,      6,       5],
            "direction":   ["AD-up"] * 7,
        })
    df = df.sort_values("fold_enrich", ascending=True).reset_index(drop=True)

    y = np.arange(len(df))
    # Color bars by whether TF is a druggable target
    druggable = {"PPARG", "IKZF1", "IRF8", "SPI1", "RUNX1", "CEBPB", "RELB"}
    colors = ["#D55E00" if t in druggable else "#AAAAAA" for t in df["tf"]]

    ax.barh(y, df["fold_enrich"], color=colors, height=0.65,
            edgecolor="white", linewidth=0.4)
    ax.axvline(1.0, color="k", linewidth=0.8)
    ax.axvline(2.0, color="grey", linewidth=0.6, linestyle="--", alpha=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(df["tf"], fontsize=9)
    ax.set_xlabel("Fold-enrichment in AD-upregulated peaks\n(observed / expected FIMO hits)", fontsize=8)

    # Annotate fold values; xlim extended to give labels room
    ax.set_xlim(0, 5.5)
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["fold_enrich"] + 0.12, i, f"{row['fold_enrich']:.1f}×",
                va="center", fontsize=7.5, color="#333333")

    ax.set_title("C  ATAC chromatin accessibility\n(motif enrichment in 56 AD-up peaks, n=1,191 DAPs)",
                 fontweight="bold", loc="left", fontsize=10)


# ── Panel D: IKZF1 convergent evidence summary ────────────────────────────────
def panel_summary(ax):
    """Formatted evidence table — avoids plotting non-commensurate metrics on shared axis."""
    rows = [
        # [Modality, Method, Value, Significance]
        ["Differential expression",  "DESeq2 pseudobulk",     "log₂FC = +1.4",   "padj < 0.001"],
        ["Regulon activity (AUCell)", "pySCENIC (5-seed)",    "AUC = 0.153",     "peak: LateAD-DAM"],
        ["Pseudotime correlation",    "Spearman",              "ρ = +0.309",      "p < 0.001"],
        ["Bulk replication",          "DESeq2 (GSE95587)",     "log₂FC = +0.638", "padj = 0.004"],
        ["Chromatin accessibility",   "FIMO motif enrichment", "3.9× in AD-up",   "56/1,191 peaks"],
        ["In-silico KO",              "CellOracle GRN",        "Δmag = 0.047",    "largest perturbation"],
    ]
    col_labels = ["Modality", "Method", "Effect Size", "Significance"]

    ax.axis("off")
    tbl = ax.table(
        cellText=rows,
        colLabels=col_labels,
        loc="center",
        cellLoc="left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 2.0)

    # Header style
    for j in range(len(col_labels)):
        tbl[0, j].set_facecolor("#2B4A6F")
        tbl[0, j].set_text_props(color="white", fontweight="bold")

    # Alternate row shading; highlight IKZF1-specific rows with faint orange
    ikzf1_rows = {0, 2, 3, 4, 5}   # rows where IKZF1 is the primary finding
    for i in range(len(rows)):
        bg = "#FFF3EC" if i in ikzf1_rows else "#F5F5F5"
        if i % 2 == 1 and i not in ikzf1_rows:
            bg = "#FFFFFF"
        for j in range(len(col_labels)):
            tbl[i + 1, j].set_facecolor(bg)

    # Column widths
    col_widths = [0.28, 0.26, 0.22, 0.24]
    for (ri, ci), cell in tbl.get_celld().items():
        cell.set_width(col_widths[ci])
        cell.set_edgecolor("#DDDDDD")

    ax.set_title("D  IKZF1 convergent validation\n(6 independent lines of evidence)",
                 fontweight="bold", loc="left", fontsize=10)
    ax.text(0.0, -0.03,
            "* log₂FC (DGE) from internal pseudobulk DESeq2, DAM vs Homeostatic.",
            transform=ax.transAxes, fontsize=7, color="#666666", va="top")


# ── Compose ───────────────────────────────────────────────────────────────────
def main():
    fig = plt.figure(figsize=(18, 14))
    gs  = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.36,
                           width_ratios=[1, 1.1])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    print("Panel A: GSE95587 forest plot...")
    panel_gse_forest(ax_a)
    print("Panel B: CellOracle KO delta...")
    panel_celloracle(ax_b)
    print("Panel C: ATAC motif enrichment...")
    panel_atac(ax_c)
    print("Panel D: IKZF1 convergent evidence...")
    panel_summary(ax_d)

    out_pdf = OUT / "fig4_validation.pdf"
    out_png = OUT / "fig4_validation.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    print(f"Saved:\n  {out_pdf}\n  {out_png}")


if __name__ == "__main__":
    main()
