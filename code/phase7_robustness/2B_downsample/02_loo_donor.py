"""
Section 2B.2 — Leave-one-donor-out robustness of the LateAD-DAM regulon peak
===========================================================================
Question: does the LateAD-DAM IKZF1(+) AUCell peak (~0.153) survive removal of
any single LateAD-DAM donor?

Input:  results/phase7/downsample/ikzf1_auc_cells.parquet
Output: results/phase7/downsample/loo_donor_lateaddam.csv
        results/phase7/downsample/loo_donor_lateaddam.png

Baseline metrics (all cells):
  - per-state mean IKZF1(+) AUCell
  - gap = mean(LateAD-DAM) - mean(Homeostatic)
  - rank of LateAD-DAM among the 6 states (1 = highest)
For each LateAD-DAM donor d: drop all of d's cells (every state), recompute.
Flag donors whose removal (a) drops LateAD-DAM out of rank 1, or (b) halves the gap.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
IN = PROJECT / "results/phase7/downsample/ikzf1_auc_cells.parquet"
OUT_CSV = PROJECT / "results/phase7/downsample/loo_donor_lateaddam.csv"
OUT_PNG = PROJECT / "results/phase7/downsample/loo_donor_lateaddam.png"

TARGET = "LateAD-DAM"
REF = "Homeostatic"


def state_means(df):
    return df.groupby("state")["ikzf1_auc"].mean()


def metrics(df):
    m = state_means(df)
    gap = m.get(TARGET, np.nan) - m.get(REF, np.nan)
    rank = int((m.sort_values(ascending=False).index.get_loc(TARGET)) + 1) if TARGET in m else np.nan
    return m.get(TARGET, np.nan), gap, rank


def main():
    df = pd.read_parquet(IN)
    base_mean, base_gap, base_rank = metrics(df)
    print(f"Baseline: LateAD-DAM mean AUC = {base_mean:.5f}  gap vs {REF} = {base_gap:.5f}  rank = {base_rank}/6")

    late_donors = sorted(df.loc[df["state"] == TARGET, "donor"].unique())
    print(f"{len(late_donors)} LateAD-DAM donors")

    rows = []
    for d in late_donors:
        sub = df[df["donor"] != d]
        n_late_cells = int((df["donor"] == d).sum())
        n_late_state = int(((df["donor"] == d) & (df["state"] == TARGET)).sum())
        mean_d, gap_d, rank_d = metrics(sub)
        rows.append({
            "donor_dropped": d,
            "n_cells_dropped": n_late_cells,
            "n_lateaddam_cells_dropped": n_late_state,
            "lateaddam_mean_auc": mean_d,
            "gap_vs_homeostatic": gap_d,
            "lateaddam_rank": rank_d,
            "rank_changed": rank_d != base_rank,
            "gap_halved": gap_d < 0.5 * base_gap,
        })
    res = pd.DataFrame(rows).sort_values("gap_vs_homeostatic")
    res.to_csv(OUT_CSV, index=False)

    n_rank_lost = int((res["lateaddam_rank"] > 1).sum()) if not np.isnan(base_rank) else np.nan
    n_gap_halved = int(res["gap_halved"].sum())
    print(f"\nLOO refits where LateAD-DAM lost rank 1: {n_rank_lost}/{len(res)}")
    print(f"LOO refits where gap more than halved: {n_gap_halved}/{len(res)}")
    print(f"gap range: [{res['gap_vs_homeostatic'].min():.5f}, {res['gap_vs_homeostatic'].max():.5f}]")
    print(f"\nMost influential donors (smallest resulting gap):\n{res.head(5).to_string(index=False)}")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axhline(base_gap, color="k", ls="--", lw=1, label=f"baseline gap = {base_gap:.4f}")
    ax.axhline(0, color="grey", lw=0.8)
    ax.bar(range(len(res)), res["gap_vs_homeostatic"], color="#4C72B0")
    ax.set_xlabel("LateAD-DAM donor removed (sorted by resulting gap)")
    ax.set_ylabel(f"mean AUCell gap\nLateAD-DAM - {REF}")
    ax.set_title("Leave-one-donor-out: LateAD-DAM IKZF1(+) regulon peak")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"\nWrote {OUT_CSV}\nWrote {OUT_PNG}")


if __name__ == "__main__":
    main()
