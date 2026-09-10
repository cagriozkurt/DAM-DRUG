"""
Section 2B.3 — Iterative per-donor downsampling of LateAD-DAM
============================================================
Question: is the LateAD-DAM IKZF1(+) AUCell peak an artefact of a few donors
contributing many cells? Cap every LateAD-DAM donor at k cells and recheck.

Input:  results/phase7/downsample/ikzf1_auc_cells.parquet
Output: results/phase7/downsample/downsample_distribution.csv
        results/phase7/downsample/downsample_distribution.png

For k in {1,2,3,4}: 1000 draws. Each draw samples <= k cells per LateAD-DAM
donor (other states untouched), recomputes mean LateAD-DAM AUCell and the gap
vs Homeostatic. Reports fraction of draws with gap > 0 and with LateAD-DAM
ranked highest.
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
OUT_CSV = PROJECT / "results/phase7/downsample/downsample_distribution.csv"
OUT_PNG = PROJECT / "results/phase7/downsample/downsample_distribution.png"

TARGET = "LateAD-DAM"
REF = "Homeostatic"
KS = [1, 2, 3, 4]
N_DRAWS = 1000
SEED = 42


def main():
    df = pd.read_parquet(IN)
    rng = np.random.default_rng(SEED)

    other_means = df[df["state"] != TARGET].groupby("state")["ikzf1_auc"].mean()
    ref_mean = other_means[REF]
    max_other = other_means.max()
    base_gap = df.loc[df["state"] == TARGET, "ikzf1_auc"].mean() - ref_mean
    print(f"Baseline LateAD-DAM gap vs {REF}: {base_gap:.5f}")
    print(f"Highest non-LateAD-DAM state mean: {max_other:.5f} ({other_means.idxmax()})")

    late = df[df["state"] == TARGET][["donor", "ikzf1_auc"]].copy()
    by_donor = {d: g["ikzf1_auc"].to_numpy() for d, g in late.groupby("donor")}

    rows = []
    draw_records = []
    for k in KS:
        gaps = np.empty(N_DRAWS)
        top = np.empty(N_DRAWS, dtype=bool)
        for i in range(N_DRAWS):
            vals = []
            for d, arr in by_donor.items():
                take = min(k, arr.size)
                idx = rng.choice(arr.size, size=take, replace=False)
                vals.append(arr[idx])
            m = np.concatenate(vals).mean()
            gaps[i] = m - ref_mean
            top[i] = m > max_other
        rows.append({
            "k_cells_per_donor": k,
            "mean_gap": gaps.mean(),
            "gap_ci_lo": np.percentile(gaps, 2.5),
            "gap_ci_hi": np.percentile(gaps, 97.5),
            "frac_gap_positive": (gaps > 0).mean(),
            "frac_lateaddam_top_state": top.mean(),
        })
        draw_records.append(pd.DataFrame({"k": k, "gap": gaps}))
        print(f"k={k}: mean gap {gaps.mean():.5f}  95%CI [{np.percentile(gaps,2.5):.5f}, "
              f"{np.percentile(gaps,97.5):.5f}]  frac gap>0 = {(gaps>0).mean():.3f}  "
              f"frac top state = {top.mean():.3f}")

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False)

    alld = pd.concat(draw_records, ignore_index=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    parts = ax.violinplot([alld.loc[alld["k"] == k, "gap"] for k in KS],
                          positions=KS, showmedians=True)
    ax.axhline(0, color="grey", lw=0.8)
    ax.axhline(base_gap, color="k", ls="--", lw=1, label=f"baseline gap = {base_gap:.4f}")
    ax.set_xlabel("cells per LateAD-DAM donor (cap)")
    ax.set_ylabel(f"mean AUCell gap  LateAD-DAM - {REF}")
    ax.set_title(f"Downsampling robustness ({N_DRAWS} draws / k)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"\nWrote {OUT_CSV}\nWrote {OUT_PNG}")


if __name__ == "__main__":
    main()
