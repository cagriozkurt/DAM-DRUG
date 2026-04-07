"""
DAM-DRUG Supplementary Figure S6 — AlphaFold2 Model Quality (pLDDT)
=====================================================================
For each target TF with an AF2 model used in docking, plots the
per-residue pLDDT scores (stored as B-factor in AF2 PDB files).
Highlights the domain used for docking (if known) and marks the
mean pocket pLDDT from fpocket analysis.

Runs locally:
  python code/phase6_figures/supp_fig_s6_af2_quality.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

PROJECT  = Path("/Volumes/PortableSSD/untitled folder/DAM-DRUG")
AF2_DIR  = PROJECT / "data/structures/af2"
OUT      = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})

# ── pLDDT colour scale (AlphaFold standard) ───────────────────────────────
# >90: very high (dark blue), 70-90: high (light blue),
# 50-70: low (yellow), <50: very low (orange)
def plddt_color(score):
    if score >= 90:   return "#1F73B7"
    if score >= 70:   return "#5DAFDF"
    if score >= 50:   return "#F0D25E"
    return "#E06A2B"

# ── TF definitions: which AF2 PDB to use, domain residue range used ───────
# domain_range: (start, end) residue numbers of the domain fragment docked
# pocket_plddt: mean pLDDT of fpocket-identified pocket residues (from step 3.2)
TF_DEFS = [
    {"tf": "IKZF1",  "pdb": "IKZF1_AF2.pdb",  "label": "IKZF1\n(ZF2 domain used)",
     "domain": (258, 338), "pocket_plddt": None, "note": "ZF2: drug_score=0.001 (undruggable)"},
    {"tf": "IRF8",   "pdb": "IRF8_AF2.pdb",   "label": "IRF8\n(DBD domain used)",
     "domain": (1, 130),   "pocket_plddt": 83,  "note": "DBD: drug_score=0.659 (druggable)"},
    {"tf": "BHLHE41","pdb": "BHLHE41_AF2.pdb","label": "BHLHE41\n(bHLH domain used)",
     "domain": (1, 100),   "pocket_plddt": None, "note": "bHLH: drug_score=0.003"},
    {"tf": "SPI1",   "pdb": "SPI1_AF2.pdb",   "label": "SPI1\n(AF2 full; ETS docked via 8EE9)",
     "domain": None,       "pocket_plddt": None, "note": "ETS domain: 8EE9 used for docking"},
    {"tf": "ACSL1",  "pdb": "ACSL1_AF2.pdb",  "label": "ACSL1\n(AF2 full-length)",
     "domain": None,       "pocket_plddt": None, "note": "Full AF2; pocket too large (~6567 Å³)"},
    {"tf": "PIK3CA", "pdb": "PIK3CA_AF2.pdb", "label": "PIK3CA\n(AF2 KD fragment)",
     "domain": None,       "pocket_plddt": None, "note": "KD docked via 4OVU"},
    {"tf": "MAF",    "pdb": "MAF_AF2.pdb",    "label": "MAF\n(bZIP domain used)",
     "domain": None,       "pocket_plddt": None, "note": "bZIP: drug_score=0.321"},
    {"tf": "RUNX1",  "pdb": "RUNX1_AF2.pdb",  "label": "RUNX1\n(Runt domain; docked via 1LJM)",
     "domain": None,       "pocket_plddt": None, "note": "Runt: 1LJM PDB used for docking"},
]

# ── PDB parser: extract CA pLDDT per residue ──────────────────────────────
def parse_plddt(pdb_path):
    """Return arrays of (residue_number, pLDDT) from CA atoms in AF2 PDB."""
    resnums, plddts = [], []
    seen = set()
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            resnum = int(line[22:26].strip())
            bfactor = float(line[60:66].strip())
            if resnum not in seen:
                resnums.append(resnum)
                plddts.append(bfactor)
                seen.add(resnum)
    return np.array(resnums), np.array(plddts)

# ── Plot ───────────────────────────────────────────────────────────────────
n_tfs  = len(TF_DEFS)
ncols  = 2
nrows  = (n_tfs + 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(12, nrows * 2.8),
                          constrained_layout=True)
axes_flat = axes.flatten()

legend_patches = [
    mpatches.Patch(color="#1F73B7", label="pLDDT ≥90 (very high)"),
    mpatches.Patch(color="#5DAFDF", label="70–90 (confident)"),
    mpatches.Patch(color="#F0D25E", label="50–70 (low)"),
    mpatches.Patch(color="#E06A2B", label="<50 (very low)"),
]

for ax, tf_def in zip(axes_flat, TF_DEFS):
    pdb_path = AF2_DIR / tf_def["pdb"]
    resnums, plddts = parse_plddt(pdb_path)

    # Colour each position by pLDDT band
    colors = [plddt_color(p) for p in plddts]
    ax.bar(resnums, plddts, color=colors, width=1.0, linewidth=0)

    # Threshold lines
    ax.axhline(90, color="#1F73B7", linewidth=0.7, linestyle="--", alpha=0.6)
    ax.axhline(70, color="#F0D25E", linewidth=0.7, linestyle="--", alpha=0.6)
    ax.axhline(50, color="#E06A2B", linewidth=0.7, linestyle="--", alpha=0.6)

    # Shade docking domain
    if tf_def["domain"] is not None:
        d0, d1 = tf_def["domain"]
        ax.axvspan(d0, d1, alpha=0.12, color="#009E73", zorder=0)
        ax.text((d0 + d1) / 2, 92, "docked\ndomain", ha="center",
                fontsize=6, color="#006340", va="bottom")

    # Mean pLDDT line
    mean_plddt = plddts.mean()
    ax.axhline(mean_plddt, color="black", linewidth=0.8, linestyle=":",
               alpha=0.8, label=f"mean={mean_plddt:.1f}")
    ax.text(resnums[-1] * 0.98, mean_plddt + 1, f"μ={mean_plddt:.1f}",
            ha="right", fontsize=6.5, color="black")

    # Pocket pLDDT annotation
    if tf_def["pocket_plddt"] is not None:
        ax.axhline(tf_def["pocket_plddt"], color="#D55E00", linewidth=1.0,
                   linestyle="-.", alpha=0.8)
        ax.text(resnums[0] + 2, tf_def["pocket_plddt"] + 1.5,
                f"pocket μ={tf_def['pocket_plddt']}",
                fontsize=6, color="#D55E00", va="bottom")

    ax.set_title(tf_def["label"], fontsize=8, fontweight="bold")
    ax.set_xlabel("Residue", fontsize=7)
    ax.set_ylabel("pLDDT", fontsize=7)
    ax.set_ylim(0, 100)
    ax.set_xlim(resnums[0], resnums[-1])
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=6.5)

    # Note
    ax.text(0.99, 0.04, tf_def["note"], transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.5, color="#555555",
            style="italic")

# Hide unused axes
for ax in axes_flat[len(TF_DEFS):]:
    ax.set_visible(False)

# Global legend in last used axis or as figure legend
fig.legend(handles=legend_patches, loc="lower right", fontsize=7.5,
           framealpha=0.9, title="pLDDT confidence", title_fontsize=7.5,
           bbox_to_anchor=(0.98, 0.02))

fig.suptitle(
    "Supplementary Figure S6 — AlphaFold2 Per-Residue pLDDT Quality for DAM-DRUG TF Models",
    fontsize=10
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"supp_fig_S6_af2_quality.{ext}",
                bbox_inches="tight", dpi=300)
    print(f"Saved → {OUT}/supp_fig_S6_af2_quality.{ext}")

plt.close(fig)
