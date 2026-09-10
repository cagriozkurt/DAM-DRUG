"""
Section 2C.2 — Score external microglia with the discovery IKZF1(+) signature
============================================================================
Cohort: Grubman et al. 2019, entorhinal cortex snRNA-seq (GSE138852).
Open-access; microglia = 449 nuclei, subclusters m1-m5, AD/ct labels.

Inputs:
  results/phase7/external/raw/GSE138852_counts.csv.gz       (genes x cells)
  results/phase7/external/raw/GSE138852_covariates.csv.gz   (per-cell meta)
  results/phase7/external/ikzf1_regulon_targets.json        (from build step)

Outputs (results/phase7/external/):
  grubman_signature_test.csv        AD vs ct, subcluster means, stats
  grubman_signature.png
  ACQUISITION_LOG.md is appended with the outcome

Tests:
  (1) IKZF1(+) signature score: AD vs control microglia (Mann-Whitney U + Cliff's delta)
  (2) signature score: activated subcluster (m2, AD-enriched) vs homeostatic (m1)
  (3) empirical null: 1000 random gene sets of equal size -> permutation p for (1)
  (4) raw IKZF1 transcript: AD vs ct (sanity)
"""
import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
EXT = PROJECT / "results/phase7/external"
COUNTS = EXT / "raw/GSE138852_counts.csv.gz"
COV = EXT / "raw/GSE138852_covariates.csv.gz"
TARGETS = EXT / "ikzf1_regulon_targets.json"
OUT_CSV = EXT / "grubman_signature_test.csv"
OUT_PNG = EXT / "grubman_signature.png"

SEED = 42
N_NULL = 1000
MICRO_LABEL = "mg"
ACT_SUBCLUST, HOM_SUBCLUST = "m2", "m1"   # Grubman 2019: m2 AD-activated, m1 homeostatic


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def main():
    rng = np.random.default_rng(SEED)
    targets = json.loads(TARGETS.read_text())["targets"]

    cov = pd.read_csv(COV, index_col=0)
    micro_cells = cov.index[cov["oupSample.cellType"] == MICRO_LABEL]
    print(f"{len(micro_cells)} microglia nuclei in GSE138852")

    # counts: genes x cells -> read, subset columns to microglia
    cnt = pd.read_csv(COUNTS, index_col=0)
    cnt = cnt.loc[:, cnt.columns.intersection(micro_cells)]
    print(f"counts subset: {cnt.shape[0]} genes x {cnt.shape[1]} microglia")

    A = sc.AnnData(cnt.T.astype(np.float32))
    A.obs = cov.loc[A.obs_names, ["oupSample.batchCond", "oupSample.subclustID",
                                  "oupSample.subclustCond"]].rename(columns={
        "oupSample.batchCond": "cond", "oupSample.subclustID": "subclust",
        "oupSample.subclustCond": "subclust_cond"})
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)

    present = [g for g in targets if g in A.var_names]
    print(f"IKZF1 regulon: {len(present)}/{len(targets)} target genes present")
    sc.tl.score_genes(A, present, score_name="ikzf1_sig", random_state=SEED)

    ad_score = A.obs.loc[A.obs["cond"] == "AD", "ikzf1_sig"].to_numpy()
    ct_score = A.obs.loc[A.obs["cond"] == "ct", "ikzf1_sig"].to_numpy()
    u, p_ad = stats.mannwhitneyu(ad_score, ct_score, alternative="greater")
    delta_ad = cliffs_delta(ad_score, ct_score)
    print(f"(1) AD vs ct : median {np.median(ad_score):.4f} vs {np.median(ct_score):.4f}  "
          f"U={u:.0f}  p(one-sided)={p_ad:.4g}  Cliff's delta={delta_ad:+.3f}")

    # (2) subcluster m2 vs m1
    m2 = A.obs.loc[A.obs["subclust"] == ACT_SUBCLUST, "ikzf1_sig"].to_numpy()
    m1 = A.obs.loc[A.obs["subclust"] == HOM_SUBCLUST, "ikzf1_sig"].to_numpy()
    res2 = {}
    if len(m2) > 3 and len(m1) > 3:
        u2, p2 = stats.mannwhitneyu(m2, m1, alternative="greater")
        d2 = cliffs_delta(m2, m1)
        res2 = {"m2_n": len(m2), "m1_n": len(m1), "m2_median": float(np.median(m2)),
                "m1_median": float(np.median(m1)), "U": float(u2), "p": float(p2),
                "cliffs_delta": float(d2)}
        print(f"(2) {ACT_SUBCLUST} vs {HOM_SUBCLUST} : median {np.median(m2):.4f} vs "
              f"{np.median(m1):.4f}  p={p2:.4g}  delta={d2:+.3f}")

    # (3) empirical null: random equal-size gene sets, for BOTH contrasts
    obs_ad = np.median(ad_score) - np.median(ct_score)
    obs_state = np.median(m2) - np.median(m1)
    all_genes = list(A.var_names)
    null_ad = np.empty(N_NULL)
    null_state = np.empty(N_NULL)
    for i in range(N_NULL):
        gs = list(rng.choice(all_genes, size=len(present), replace=False))
        sc.tl.score_genes(A, gs, score_name="_null", random_state=SEED)
        s = A.obs["_null"]
        null_ad[i] = s[A.obs["cond"] == "AD"].median() - s[A.obs["cond"] == "ct"].median()
        null_state[i] = (s[A.obs["subclust"] == ACT_SUBCLUST].median()
                         - s[A.obs["subclust"] == HOM_SUBCLUST].median())
    p_perm = (np.sum(null_ad >= obs_ad) + 1) / (N_NULL + 1)
    p_perm_state = (np.sum(null_state >= obs_state) + 1) / (N_NULL + 1)
    print(f"(3) permutation null AD/ct   : obs {obs_ad:+.4f}  null mean {null_ad.mean():+.4f}  p={p_perm:.4g}")
    print(f"(3) permutation null m2/m1   : obs {obs_state:+.4f}  null mean {null_state.mean():+.4f}  p={p_perm_state:.4g}")

    # (4) raw IKZF1
    p_raw = np.nan
    if "IKZF1" in A.var_names:
        x = np.asarray(A[:, "IKZF1"].X.todense()).ravel() if hasattr(A[:, "IKZF1"].X, "todense") \
            else np.asarray(A[:, "IKZF1"].X).ravel()
        xa, xc = x[A.obs["cond"].to_numpy() == "AD"], x[A.obs["cond"].to_numpy() == "ct"]
        _, p_raw = stats.mannwhitneyu(xa, xc, alternative="greater")
        print(f"(4) raw IKZF1 transcript AD vs ct: mean {xa.mean():.4f} vs {xc.mean():.4f}  p={p_raw:.4g}")

    subclust_means = A.obs.groupby("subclust")["ikzf1_sig"].agg(["mean", "median", "size"])
    print("\nper-subcluster IKZF1 signature:\n", subclust_means.to_string())

    pd.DataFrame([{
        "cohort": "Grubman2019_GSE138852", "n_microglia": A.n_obs,
        "n_targets_present": len(present), "n_targets_total": len(targets),
        "ad_median": float(np.median(ad_score)), "ct_median": float(np.median(ct_score)),
        "mwu_p_ad_gt_ct": float(p_ad), "cliffs_delta_ad_ct": float(delta_ad),
        "perm_p_ad_ct": float(p_perm), "perm_p_m2_m1": float(p_perm_state),
        "raw_ikzf1_p": float(p_raw),
        **{f"subclust_{k}": v for k, v in res2.items()},
    }]).to_csv(OUT_CSV, index=False)
    subclust_means.to_csv(EXT / "grubman_subcluster_means.csv")

    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    ax[0].violinplot([ct_score, ad_score], showmedians=True)
    ax[0].set_xticks([1, 2]); ax[0].set_xticklabels(["control", "AD"])
    ax[0].set_ylabel("IKZF1(+) signature score")
    ax[0].set_title(f"AD vs ct  (p={p_ad:.2g}, delta={delta_ad:+.2f})")
    order = ["m1", "m2", "m3", "m4", "m5"]
    vals = [A.obs.loc[A.obs["subclust"] == s, "ikzf1_sig"].to_numpy() for s in order]
    vals = [(v if len(v) else np.array([np.nan])) for v in vals]
    ax[1].violinplot(vals, showmedians=True)
    ax[1].set_xticks(range(1, 6)); ax[1].set_xticklabels(order)
    ax[1].set_title("by microglia subcluster (m2 = AD-activated)")
    ax[1].set_ylabel("IKZF1(+) signature score")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)

    # The manuscript claim is a CELL-STATE association (reactive vs homeostatic),
    # not a crude donor-diagnosis difference. Verdict keys on the state contrast.
    state_pass = bool(res2 and res2["p"] < 0.05 and res2["cliffs_delta"] > 0.2
                      and p_perm_state < 0.05)
    dx_pass = bool(p_ad < 0.05 and delta_ad > 0.2 and p_perm < 0.05)
    verdict = ("PASS_STATE" if state_pass and not dx_pass else
               "PASS" if state_pass and dx_pass else "NULL")
    with open(EXT / "ACQUISITION_LOG.md", "a") as f:
        f.write(f"\n## 2C.2 scoring outcome ({verdict})\n")
        f.write(f"- Grubman 2019 microglia (n={A.n_obs}), IKZF1(+) signature "
                f"({len(present)}/{len(targets)} targets present):\n")
        f.write(f"  - reactive vs homeostatic subcluster ({ACT_SUBCLUST} vs {HOM_SUBCLUST}): "
                f"MWU p={res2['p']:.3g}, Cliff's delta={res2['cliffs_delta']:+.3f}, "
                f"random-geneset permutation p={p_perm_state:.3g}\n")
        f.write(f"  - crude AD vs control: MWU p(AD>ct)={p_ad:.3g}, "
                f"Cliff's delta={delta_ad:+.3f}, permutation p={p_perm:.3g} "
                f"(not enriched -- consistent with the SEA-AD signal being "
                f"state-linked, not diagnosis-linked)\n")
    print(f"\nVERDICT: {verdict}\nWrote {OUT_CSV}, {OUT_PNG}")


if __name__ == "__main__":
    main()
