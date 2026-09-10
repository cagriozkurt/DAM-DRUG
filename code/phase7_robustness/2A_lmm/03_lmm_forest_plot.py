"""
Section 2A - Forest plot of LMM state coefficients
=================================================
Input:  results/phase7/lmm/<gene>_lmm.csv  (from 02_fit_lmm.R)
Output: results/phase7/lmm/ikzf1_state_forest.{pdf,png}

IKZF1 highlighted; IRF8/SPI1/BHLHE41/RUNX1/PPARG/CEBPB as small multiples.
Shows the 'min10' full-model state coefficients (ref = Homeostatic) with 95% CI.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
LMM = PROJECT / "results/phase7/lmm"

GENES = ["IKZF1", "IRF8", "SPI1", "BHLHE41", "RUNX1", "PPARG", "CEBPB"]
STATE_ORDER = ["stateDAM", "stateIRM", "stateLAM", "stateLateAD-DAM"]
STATE_LABEL = {"stateDAM": "DAM", "stateIRM": "IRM", "stateLAM": "LAM",
               "stateLateAD-DAM": "LateAD-DAM"}
VARIANT = "min10"


def load(gene):
    f = LMM / f"{gene}_lmm.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d.columns = [c.strip('"') for c in d.columns]
    for c in ("gene", "variant", "term"):
        d[c] = d[c].astype(str).str.strip('"')
    return d[d["variant"] == VARIANT]


def main():
    fig, axes = plt.subplots(1, len(GENES), figsize=(2.4 * len(GENES), 3.4),
                             sharey=True)
    for ax, gene in zip(axes, GENES):
        d = load(gene)
        ax.axvline(0, color="grey", lw=0.8, zorder=0)
        if d is None or d.empty:
            ax.set_title(f"{gene}\n(no fit)", fontsize=10)
            continue
        ys = np.arange(len(STATE_ORDER))[::-1]
        for y, term in zip(ys, STATE_ORDER):
            row = d[d["term"] == term]
            if row.empty:
                continue
            est = float(row["estimate"].iloc[0])
            lo, hi = float(row["ci_lo"].iloc[0]), float(row["ci_hi"].iloc[0])
            sig = lo > 0 or hi < 0
            color = "#C44E52" if gene == "IKZF1" else "#4C72B0"
            color = color if sig else "#AAAAAA"
            ax.plot([lo, hi], [y, y], color=color, lw=2)
            ax.plot(est, y, "o", color=color, ms=6 if gene == "IKZF1" else 5)
        ax.set_yticks(ys)
        ax.set_yticklabels([STATE_LABEL[t] for t in STATE_ORDER], fontsize=9)
        ttl_w = "bold" if gene == "IKZF1" else "normal"
        ax.set_title(gene, fontsize=11, fontweight=ttl_w)
        ax.set_xlabel("log1p(CPM) vs\nHomeostatic", fontsize=8)
        ax.tick_params(labelsize=8)

    fig.suptitle("LMM state effects with (1|donor) + (1|region)   [>=10 cells/group]",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    for ext in ("pdf", "png"):
        fig.savefig(LMM / f"ikzf1_state_forest.{ext}", dpi=200)
    print(f"Wrote {LMM}/ikzf1_state_forest.pdf / .png")


if __name__ == "__main__":
    main()
