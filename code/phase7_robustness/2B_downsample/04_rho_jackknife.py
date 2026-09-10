"""
Section 2B.4 — Donor jackknife / bootstrap of the IKZF1(+) AUCell-vs-pseudotime rho
==================================================================================
Baseline Spearman rho(IKZF1(+) AUCell, dpt_pseudotime) = +0.309 (manuscript;
recomputed here from the same 100K-cell AUCell subsample).

Input:  results/phase7/downsample/ikzf1_auc_cells.parquet
Output: results/phase7/downsample/rho_jackknife.csv   (one row per dropped donor)
        results/phase7/downsample/rho_bootstrap.csv   (summary)
        results/phase7/downsample/rho_jackknife.png

(1) leave-one-donor-out over ALL donors: recompute rho.
(2) 1000 donor-level bootstrap resamples (resample donors with replacement,
    pool their cells): rho distribution + 95% percentile CI.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
IN = PROJECT / "results/phase7/downsample/ikzf1_auc_cells.parquet"
OUT_JK = PROJECT / "results/phase7/downsample/rho_jackknife.csv"
OUT_BS = PROJECT / "results/phase7/downsample/rho_bootstrap.csv"
OUT_PNG = PROJECT / "results/phase7/downsample/rho_jackknife.png"

N_BOOT = 1000
SEED = 42


def spearman(a, b):
    return stats.spearmanr(a, b).statistic


def main():
    df = pd.read_parquet(IN).dropna(subset=["dpt", "ikzf1_auc"])
    auc = df["ikzf1_auc"].to_numpy()
    dpt = df["dpt"].to_numpy()
    donor = df["donor"].to_numpy()
    base = spearman(auc, dpt)
    print(f"Baseline rho (n={len(df):,}): {base:+.5f}")

    donors = np.array(sorted(df["donor"].unique()))

    # (1) jackknife
    jk = []
    for d in donors:
        mask = donor != d
        jk.append({"donor_dropped": d, "n_cells_dropped": int((~mask).sum()),
                   "rho": spearman(auc[mask], dpt[mask])})
    jk_df = pd.DataFrame(jk).sort_values("rho")
    jk_df["delta_vs_base"] = jk_df["rho"] - base
    jk_df.to_csv(OUT_JK, index=False)
    print(f"jackknife rho range: [{jk_df['rho'].min():+.5f}, {jk_df['rho'].max():+.5f}]  "
          f"sign stable: {(np.sign(jk_df['rho']) == np.sign(base)).all()}")

    # (2) donor bootstrap
    rng = np.random.default_rng(SEED)
    cells_by_donor = {d: np.where(donor == d)[0] for d in donors}
    boot = np.empty(N_BOOT)
    for i in range(N_BOOT):
        pick = rng.choice(donors, size=len(donors), replace=True)
        idx = np.concatenate([cells_by_donor[d] for d in pick])
        boot[i] = spearman(auc[idx], dpt[idx])
    ci = np.percentile(boot, [2.5, 97.5])
    bs_df = pd.DataFrame([{
        "baseline_rho": base,
        "boot_mean": boot.mean(),
        "boot_ci_lo": ci[0],
        "boot_ci_hi": ci[1],
        "frac_positive": (boot > 0).mean(),
        "n_boot": N_BOOT,
    }])
    bs_df.to_csv(OUT_BS, index=False)
    print(f"bootstrap rho: mean {boot.mean():+.5f}  95%CI [{ci[0]:+.5f}, {ci[1]:+.5f}]  "
          f"frac > 0 = {(boot > 0).mean():.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].bar(range(len(jk_df)), jk_df["rho"], color="#4C72B0")
    ax[0].axhline(base, color="k", ls="--", lw=1, label=f"baseline {base:+.3f}")
    ax[0].set_title("Leave-one-donor-out rho")
    ax[0].set_xlabel("donor removed (sorted)")
    ax[0].set_ylabel("Spearman rho")
    ax[0].legend(fontsize=8)
    ax[1].hist(boot, bins=40, color="#4C72B0")
    ax[1].axvline(base, color="k", ls="--", lw=1)
    ax[1].axvline(0, color="red", lw=1)
    ax[1].set_title(f"Donor bootstrap rho ({N_BOOT})")
    ax[1].set_xlabel("Spearman rho")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"\nWrote {OUT_JK}\nWrote {OUT_BS}\nWrote {OUT_PNG}")


if __name__ == "__main__":
    main()
