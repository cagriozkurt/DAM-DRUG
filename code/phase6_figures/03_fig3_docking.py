"""
DAM-DRUG Figure 3 — Structure-Based Drug Screening Results
===========================================================
Panels:
  A. Tier-1 MM-GBSA bar chart (26 Tier-1 hits, colored by target)
  B. Tier-2 top-6 MM-GBSA waterfall (validated hits)
  C. Selectivity plot: ΔΔG vs off-targets for each hit
  D. Summary scorecard table (top hits: target / drug / GBSA / SI)

Run locally:
  python code/phase6_figures/03_fig3_docking.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
PH4     = PROJECT / "results/phase4"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9,
    "axes.linewidth": 0.8, "pdf.fonttype": 42,
})

TARGET_COLORS = {
    "IKZF1_8RQC":       "#AA3377",
    "IKZF1_8RQC_CRBN":  "#AA3377",
    "IRF8_AF2_DBD":     "#EE6677",
    "IRF8_AF2_DBD_prep":"#EE6677",
    "BHLHE41_AF2_bHLH": "#CCBB44",
    "BHLHE41_AF2_bHLH_prep": "#CCBB44",
    "PPARG_1FM9_LBD":   "#228833",
    "PPARG_1FM9_LBD_prep": "#228833",
    "RUNX1_1LJM_Runt":  "#4477AA",
    "RUNX1_1LJM_Runt_prep": "#4477AA",
    "MAF_4EOT_bZIP":    "#66CCEE",
    "MAF_4EOT_bZIP_prep":"#66CCEE",
}
TARGET_LABELS = {
    "IKZF1_8RQC":       "IKZF1",
    "IKZF1_8RQC_CRBN":  "IKZF1",
    "IRF8_AF2_DBD_prep":"IRF8",
    "BHLHE41_AF2_bHLH_prep": "BHLHE41",
    "PPARG_1FM9_LBD_prep": "PPARG",
    "RUNX1_1LJM_Runt_prep": "RUNX1",
    "MAF_4EOT_bZIP_prep":"MAF",
}


# ── Panel A: Tier-1 MM-GBSA bar ───────────────────────────────────────────────
def panel_tier1(ax):
    df = pd.read_csv(PH4 / "mmpbsa/mmpbsa_summary.csv")
    df = df.sort_values("dg_gbsa").head(15).reset_index(drop=True)  # top 15; full list → supplementary
    df["color"] = df["target"].map(TARGET_COLORS).fillna("grey")
    df["label"] = df["target"].map(TARGET_LABELS).fillna(df["target"])

    # Clean up chembl IDs for display (use short ID)
    df["disp"] = df["chembl_id"].str.replace("CHEMBL", "")

    bars = ax.barh(range(len(df)), df["dg_gbsa"],
                   xerr=df["sem"], color=df["color"],
                   error_kw={"linewidth": 0.8, "capsize": 2},
                   height=0.7, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="k", linewidth=0.8)
    ax.axvline(-5, color="grey", linewidth=0.6, linestyle="--", alpha=0.5)

    ax.set_yticks(range(len(df)))
    # Strip redundant target name from y-labels — color coding handles it
    ax.set_yticklabels(
        [f"CHEMBL{r['disp']}" for _, r in df.iterrows()],
        fontsize=7
    )
    ax.set_xlabel("ΔG_bind MM-GBSA (kcal/mol)", fontsize=8)
    ax.set_title("A  Tier-1 MM-GBSA top 15 hits (full n=26 → Supplementary)",
                 fontweight="bold", loc="left", fontsize=10)

    # Legend for targets — outside plot to avoid bar occlusion
    seen = {}
    handles = []
    for t, c in TARGET_COLORS.items():
        lbl = TARGET_LABELS.get(t, t)
        if lbl not in seen and lbl in df["label"].values:
            seen[lbl] = True
            handles.append(mpatches.Patch(color=c, label=lbl))
    ax.legend(handles=handles, fontsize=7, frameon=False,
              bbox_to_anchor=(1.01, 1), loc="upper left",
              title="Target", title_fontsize=7)


# ── Panel B: Tier-2 top-6 waterfall ───────────────────────────────────────────
def panel_tier2(ax):
    df = pd.read_csv(PH4 / "vina_selectivity/selectivity_table.csv")

    # Sort by MM-GBSA
    df = df.sort_values("mmpbsa_dg").reset_index(drop=True)

    colors = []
    for _, row in df.iterrows():
        c = TARGET_COLORS.get(row["on_target"], "#888888")
        colors.append(c)

    bars = ax.bar(range(len(df)), df["mmpbsa_dg"],
                  color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.axhline(-10, color="grey", linewidth=0.6, linestyle="--", alpha=0.6)

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["pref_name"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("ΔG_bind MM-GBSA (kcal/mol)", fontsize=8)
    ax.set_title("B  Tier-2 validated hits (FDA-approved, n=6)",
                 fontweight="bold", loc="left", fontsize=10)

    # Annotate values
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(i, row["mmpbsa_dg"] - 0.8, f"{row['mmpbsa_dg']:.1f}",
                ha="center", va="top", fontsize=7, color="white", fontweight="bold")

    # Legend for target colors (replaces in-plot text labels)
    seen = {}
    handles = []
    for t, c in TARGET_COLORS.items():
        lbl = TARGET_LABELS.get(t, t)
        if lbl not in seen and lbl in [TARGET_LABELS.get(r["on_target"], r["on_target"]) for _, r in df.iterrows()]:
            seen[lbl] = True
            handles.append(mpatches.Patch(color=c, label=lbl))
    if handles:
        ax.legend(handles=handles, fontsize=7, frameon=False,
                  bbox_to_anchor=(1.01, 1), loc="upper left",
                  title="Target", title_fontsize=7)


# ── Panel C: Selectivity ΔΔG ─────────────────────────────────────────────────
def panel_selectivity(ax):
    df = pd.read_csv(PH4 / "vina_selectivity/selectivity_table.csv")

    off_targets = ["DRD2", "HTR2A", "hERG"]
    vina_cols   = ["vina_drd2", "vina_htr2a", "vina_herg"]
    vina_on_col = "vina_ontarget"
    # CVD-safe colors per off-target — direction (pos/neg) read from bar height
    OT_COLORS  = {"DRD2": "#E69F00", "HTR2A": "#56B4E9", "hERG": "#CC79A7"}
    # Hatching per off-target so bars are distinguishable without colour
    OT_HATCHES = {"DRD2": "",      "HTR2A": "///",   "hERG": "..."}

    x = np.arange(len(df))
    width = 0.25
    offsets = [-width, 0, width]

    for i, (ot, col) in enumerate(zip(off_targets, vina_cols)):
        ddg = df[vina_on_col] - df[col]
        ax.bar(x + offsets[i], ddg, width=width * 0.9,
               color=OT_COLORS[ot], alpha=0.85, label=ot,
               hatch=OT_HATCHES[ot], edgecolor="#333333", linewidth=0.8)

    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(df["pref_name"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("ΔΔG_Vina (on-target − off-target, kcal/mol)", fontsize=8)
    ax.set_title("C  Selectivity: Vina ΔΔG vs off-targets\n(positive = on-target preferred)",
                 fontweight="bold", loc="left", fontsize=10)

    # Legend: colour + hatch pattern per off-target
    ot_handles = [mpatches.Patch(facecolor=OT_COLORS[ot], hatch=OT_HATCHES[ot],
                                 edgecolor="#333333", label=ot) for ot in off_targets]
    ax.legend(handles=ot_handles, fontsize=7, frameon=False,
              title="Off-target", title_fontsize=7,
              bbox_to_anchor=(1.01, 1), loc="upper left")


# ── Panel D: Summary scorecard ────────────────────────────────────────────────
def panel_scorecard(ax):
    df = pd.read_csv(PH4 / "vina_selectivity/selectivity_table.csv")
    df = df.sort_values("mmpbsa_dg").reset_index(drop=True)

    tgt_labels = [TARGET_LABELS.get(t, t) for t in df["on_target"]]

    table_data = []
    for i, row in df.iterrows():
        tgt  = TARGET_LABELS.get(row["on_target"], row["on_target"])
        flag = "⚠" if row["pref_name"] == "CAPMATINIB" else ""
        table_data.append([
            row["pref_name"],
            tgt,
            f"{row['mmpbsa_dg']:.1f}",
            f"{row['min_si']:.2f}",
            flag,
        ])

    col_labels = ["Drug", "Target", "ΔG_bind\n(kcal/mol)", "Min SI", "Note"]
    ax.axis("off")
    tbl = ax.table(
        cellText=table_data,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 1.8)
    # Manual column widths: Drug wide, Target wide, rest normal
    col_widths = [0.22, 0.20, 0.15, 0.18, 0.08]
    for (row_idx, col_idx), cell in tbl.get_celld().items():
        cell.set_width(col_widths[col_idx])

    # Style header
    for j in range(len(col_labels)):
        tbl[0, j].set_facecolor("#2B4A6F")
        tbl[0, j].set_text_props(color="white", fontweight="bold")

    # Alternate row shading (neutral) — no target-color fill to avoid pink
    for i in range(len(df)):
        bg = "#F5F5F5" if i % 2 == 0 else "#FFFFFF"
        for j in range(len(col_labels)):
            tbl[i + 1, j].set_facecolor(bg)

    ax.set_title("D  Top hit scorecard", fontweight="bold", loc="left", fontsize=10)
    ax.text(0.0, -0.04,
            "Min SI = min(SI_DRD2, SI_5HT2A, SI_hERG); SI = ΔG_off-target / ΔG_on-target (Vina scores).  "
            "⚠ = CNS penetration concern (P-gp substrate, low logBB).",
            transform=ax.transAxes, fontsize=7, color="#555555", va="top", wrap=True)


# ── Compose ───────────────────────────────────────────────────────────────────
def main():
    fig = plt.figure(figsize=(18, 14))
    gs  = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32,
                            height_ratios=[1.4, 1])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    print("Panel A: Tier-1 MM-GBSA...")
    panel_tier1(ax_a)
    print("Panel B: Tier-2 waterfall...")
    panel_tier2(ax_b)
    print("Panel C: Selectivity ΔΔG...")
    panel_selectivity(ax_c)
    print("Panel D: Scorecard...")
    panel_scorecard(ax_d)

    out_pdf = OUT / "fig3_docking.pdf"
    out_png = OUT / "fig3_docking.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    print(f"Saved:\n  {out_pdf}\n  {out_png}")


if __name__ == "__main__":
    main()
