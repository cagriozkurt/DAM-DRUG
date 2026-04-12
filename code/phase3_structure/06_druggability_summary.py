"""
DAM-DRUG Phase 3 — Druggability ranking and summary plots
==========================================================
Reads pocket_summary.csv, selects best pocket per TF,
ranks targets, and generates summary figures.
"""

import os
import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
RES3     = PROJECT / "results/phase3"
RES3.mkdir(parents=True, exist_ok=True)
IN_CSV   = RES3 / "pocket_summary.csv"
RANK_CSV = RES3 / "druggability_ranking.csv"
FIG_OUT  = RES3 / "pocket_summary_plots.png"

df = pd.read_csv(IN_CSV)
log.info(f"Loaded {len(df)} structures from {IN_CSV.name}")

# ── Per-TF best structure selection ───────────────────────────────────────────
# Prefer PDB over AF2 (more reliable geometry).
# For PPARG, prefer 2PRG (apo, unbiased pocket geometry) for ranking;
# note 1FM9 volume separately as the ligand-expanded cavity.

def best_row(tf_df: pd.DataFrame) -> pd.Series:
    """Select best structure for a TF: prefer PDB, then highest drug_score."""
    pdb_rows = tf_df[tf_df["source"] == "pdb"]
    if not pdb_rows.empty:
        return pdb_rows.sort_values("drug_score", ascending=False).iloc[0]
    return tf_df.sort_values("drug_score", ascending=False).iloc[0]

ranking_rows = []
TF_ORDER = ["PPARG", "IKZF1", "IRF8", "BHLHE41", "RUNX1", "RUNX2", "MAF",
            "SPI1", "ACSL1", "PIK3CA"]

for tf in TF_ORDER:
    sub = df[df["tf"] == tf]
    if sub.empty:
        log.warning(f"  {tf}: no data in pocket_summary.csv")
        continue

    best = best_row(sub)
    row = best.to_dict()

    # For PPARG: also record 1FM9 volume (ligand-bound cavity size)
    if tf == "PPARG":
        fm9 = sub[sub["pdb_id"] == "1FM9"]
        row["note_1fm9_volume"] = fm9["volume_A3"].values[0] if not fm9.empty else np.nan

    # Flag low-confidence AF2 pockets
    if best["source"] == "af2" and not np.isnan(best.get("mean_plddt_pocket", np.nan)):
        row["af2_pocket_plddt_flag"] = "LOW" if best["mean_plddt_pocket"] < 70 else "OK"
    else:
        row["af2_pocket_plddt_flag"] = "N/A"

    ranking_rows.append(row)
    log.info(f"  {tf}: best={best['structure_id']}  "
             f"drug_score={best['drug_score']:.3f}  vol={best['volume_A3']:.0f}Å³")

rank_df = pd.DataFrame(ranking_rows)
rank_df = rank_df.sort_values("drug_score", ascending=False).reset_index(drop=True)
rank_df.index += 1  # rank from 1
rank_df.to_csv(RANK_CSV, index_label="rank")
log.info(f"\nSaved druggability ranking → {RANK_CSV}")

# ── Plots ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle("DAM-DRUG Phase 3 — TF Druggability Assessment", fontweight="bold")

tf_colors = {"PPARG": "#2196F3", "IKZF1": "#FF5722", "IRF8": "#4CAF50",
             "BHLHE41": "#9C27B0", "RUNX1": "#FF9800",
             "RUNX2": "#795548", "MAF": "#009688",
             "SPI1": "#E91E63", "ACSL1": "#607D8B", "PIK3CA": "#CDDC39"}
source_marker = {"pdb": "o", "af2": "^"}

# Plot 1: Drug score bar chart (all structures, grouped by TF)
ax1 = axes[0]
x_labels, x_pos, colors = [], [], []
bar_vals = []
for i, row in df.iterrows():
    x_labels.append(row["structure_id"].replace("_trimmed", "").replace("_LBD", "\nLBD"))
    bar_vals.append(row["drug_score"] if pd.notna(row["drug_score"]) else 0)
    colors.append(tf_colors.get(row["tf"], "#999999"))

x = range(len(bar_vals))
bars = ax1.bar(x, bar_vals, color=colors, edgecolor="white", linewidth=0.5)
ax1.set_xticks(list(x))
ax1.set_xticklabels(x_labels, rotation=75, ha="right", fontsize=7)
ax1.set_ylabel("fpocket Drug Score")
ax1.set_title("Drug Score by Structure")
ax1.axhline(0.5, color="red", linestyle="--", linewidth=1, label="score=0.5 threshold")
ax1.legend(fontsize=8)

patches = [mpatches.Patch(color=c, label=tf) for tf, c in tf_colors.items()]
ax1.legend(handles=patches, fontsize=8, loc="upper right")

# Plot 2: Pocket volume vs hydrophobicity scatter
ax2 = axes[1]
for _, row in df.iterrows():
    if pd.isna(row.get("volume_A3")) or pd.isna(row.get("hydrophobicity")):
        continue
    marker = source_marker.get(row["source"], "s")
    ax2.scatter(row["volume_A3"], row["hydrophobicity"],
                c=tf_colors.get(row["tf"], "#999"), marker=marker,
                s=80, alpha=0.8, edgecolors="k", linewidths=0.5)
    ax2.annotate(row["structure_id"].split("_")[1][:4],
                 (row["volume_A3"], row["hydrophobicity"]),
                 fontsize=6, xytext=(3, 3), textcoords="offset points")

ax2.set_xlabel("Pocket Volume (Å³)")
ax2.set_ylabel("Hydrophobicity Score")
ax2.set_title("Pocket Volume vs. Hydrophobicity")

patches = [mpatches.Patch(color=c, label=tf) for tf, c in tf_colors.items()]
legend1 = ax2.legend(handles=patches, fontsize=8, loc="upper right")
circle = plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
                    markersize=8, label="PDB")
tri   = plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="gray",
                    markersize=8, label="AF2")
ax2.legend(handles=[circle, tri], fontsize=8, loc="lower right")
ax2.add_artist(legend1)

plt.tight_layout()
plt.savefig(FIG_OUT, dpi=150, bbox_inches="tight")
log.info(f"Saved figure → {FIG_OUT}")

# ── Print ranking table ────────────────────────────────────────────────────────
log.info("\n" + "=" * 70)
log.info("DRUGGABILITY RANKING")
log.info("=" * 70)
cols = ["tf", "structure_id", "drug_score", "pocket_score", "volume_A3",
        "n_pockets_total", "af2_pocket_plddt_flag"]
log.info(rank_df[[c for c in cols if c in rank_df.columns]].to_string())
log.info(f"\nNext step: run 01_docking_prep.py")
