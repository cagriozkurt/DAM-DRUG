"""
Section 2A.1 — Multi-region pseudobulk for regional-confounding LMMs
==================================================================
Builds donor x substate x brain-region pseudobulk from the trajectory object
(raw integer counts in X, 36,601 genes, all 10 SEA-AD regions, 236,002 microglia).

Input:  results/phase1/trajectory/microglia_trajectory.h5ad
Outputs (results/phase7/lmm/):
  pseudobulk_multiregion.csv      long: group_id, donor, state, region, gene, count, cpm, log1p_cpm
  pseudobulk_group_meta.csv       group_id, donor, state, region, n_cells, total_umi,
                                  age, age_z, sex, braak, adnc, cognitive_status
  state_x_region_group_counts.csv surviving groups per state under each filter

Two variants are flagged per group via `passes_min10` (>= 10 cells).
DAM-IRM is dropped (analytically ambiguous hybrid, per manuscript primary DGE).
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
H5AD = PROJECT / "results/phase1/trajectory/microglia_trajectory.h5ad"
OUT = PROJECT / "results/phase7/lmm"
OUT.mkdir(parents=True, exist_ok=True)

TARGET_GENES = ["IKZF1", "IKZF2", "IKZF3", "IRF8", "SPI1", "BHLHE41",
                "RUNX1", "CEBPB", "PPARG", "ACTB", "GAPDH", "B2M"]
DROP_STATES = {"DAM-IRM"}
MIN_CELLS = 10
CHUNK = 20000


def main():
    A = ad.read_h5ad(H5AD, backed="r")
    var_names = list(A.var_names)
    n_obs = A.n_obs
    print(f"{n_obs:,} cells x {A.n_vars:,} genes")

    genes = [g for g in TARGET_GENES if g in var_names]
    missing = set(TARGET_GENES) - set(genes)
    if missing:
        print(f"  WARNING: genes not found, skipped: {sorted(missing)}")
    gidx = [var_names.index(g) for g in genes]

    # --- per-cell total UMI (all genes) + target-gene counts, chunked ---
    total_umi = np.zeros(n_obs, dtype=np.float64)
    tgt = np.zeros((n_obs, len(genes)), dtype=np.float64)
    for start in range(0, n_obs, CHUNK):
        stop = min(start + CHUNK, n_obs)
        X = A.X[start:stop]
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        if not np.allclose(X.data, np.round(X.data)):
            raise SystemExit("X is not raw integer counts — aborting")
        total_umi[start:stop] = np.asarray(X.sum(axis=1)).ravel()
        tgt[start:stop, :] = X[:, gidx].toarray()
        print(f"  rows {start:,}-{stop:,}", end="\r")
    print()

    obs = A.obs.copy()
    A.file.close()

    meta = pd.DataFrame({
        "donor": obs["Donor ID"].astype(str).values,
        "state": obs["state"].astype(str).values,
        "region": obs["Brain Region"].astype(str).values,
        "age": pd.to_numeric(obs["Age at Death"], errors="coerce").values,
        "sex": obs["Sex"].astype(str).values,
        "braak": obs["Braak"].astype(str).values,
        "adnc": obs["Overall AD neuropathological Change"].astype(str).values,
        "cognitive_status": obs["Cognitive Status"].astype(str).values,
        "total_umi": total_umi,
    })
    for j, g in enumerate(genes):
        meta[f"__cnt_{g}"] = tgt[:, j]

    keep = ~meta["state"].isin(DROP_STATES) & meta["state"].notna() & (meta["state"] != "nan")
    meta = meta[keep].reset_index(drop=True)
    print(f"after dropping {DROP_STATES}: {len(meta):,} cells, states = {sorted(meta['state'].unique())}")

    grp_keys = ["donor", "state", "region"]
    agg = {"total_umi": "sum", "age": "first", "sex": "first",
           "braak": "first", "adnc": "first", "cognitive_status": "first"}
    for g in genes:
        agg[f"__cnt_{g}"] = "sum"
    gb = meta.groupby(grp_keys, observed=True)
    grp = gb.agg(agg)
    grp["n_cells"] = gb.size()
    grp = grp.reset_index()
    grp["group_id"] = (grp["donor"] + "|" + grp["state"] + "|" + grp["region"])
    grp["age_z"] = (grp["age"] - grp["age"].mean()) / grp["age"].std(ddof=0)
    grp["passes_min10"] = grp["n_cells"] >= MIN_CELLS
    print(f"{len(grp):,} donor x state x region groups; {grp['passes_min10'].sum():,} pass >= {MIN_CELLS} cells")

    # --- long table ---
    long_rows = []
    for _, r in grp.iterrows():
        for g in genes:
            cnt = r[f"__cnt_{g}"]
            cpm = cnt / r["total_umi"] * 1e6 if r["total_umi"] > 0 else np.nan
            long_rows.append({
                "group_id": r["group_id"], "donor": r["donor"], "state": r["state"],
                "region": r["region"], "gene": g, "count": cnt, "cpm": cpm,
                "log1p_cpm": np.log1p(cpm), "n_cells": r["n_cells"],
                "passes_min10": r["passes_min10"],
            })
    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(OUT / "pseudobulk_multiregion.csv", index=False)

    meta_cols = ["group_id", "donor", "state", "region", "n_cells", "total_umi",
                 "age", "age_z", "sex", "braak", "adnc", "cognitive_status", "passes_min10"]
    grp[meta_cols].to_csv(OUT / "pseudobulk_group_meta.csv", index=False)

    # --- coverage summary ---
    cov = []
    for st, s in grp.groupby("state", observed=True):
        cov.append({
            "state": st,
            "n_groups_unconstrained": len(s),
            "n_regions_unconstrained": s["region"].nunique(),
            "n_donors_unconstrained": s["donor"].nunique(),
            "n_groups_min10": int(s["passes_min10"].sum()),
            "n_regions_min10": s.loc[s["passes_min10"], "region"].nunique(),
            "n_donors_min10": s.loc[s["passes_min10"], "donor"].nunique(),
        })
    cov_df = pd.DataFrame(cov)
    cov_df.to_csv(OUT / "state_x_region_group_counts.csv", index=False)
    print("\nCoverage by state:\n", cov_df.to_string(index=False))
    print(f"\nWrote 3 files to {OUT}")


if __name__ == "__main__":
    main()
