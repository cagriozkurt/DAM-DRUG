"""
WP1.4 — Benchmark the JASPAR-2026 rescued regulons and run negative null models
=============================================================================
Runs inside scenic.sif (loompy, scipy, pandas, numpy, statsmodels).

Inputs (results/phase7/grn_jaspar2026/ unless noted):
  scenic_auc_jaspar2026.loom            per-cell AUCell for the rescued regulons
  regulons_jaspar2026.csv               pyscenic ctx output (for target-gene lists)
  ../../phase2/GRN/microglia_raw.loom   per-cell expression (IKZF1/2/3 traces, target-set mean)
  ../../phase1/trajectory/pseudotime.csv
  data/resources/jaspar2026/tf_tg_curated.tsv   FIXME: JASPAR 2026 curated TF-TG (cols: TF, target)

Outputs (results/phase7/grn_jaspar2026/):
  benchmark_regulon_pseudotime.csv     per regulon: rho, padj, peak_state, is_target_TF
  benchmark_tf_ranking.csv             IKZF1 vs BHLHE41/IRF8/SPI1/RUNX1/CEBPB/PPARG
  null_permutation_fdr.csv             1000-perm empirical FDR for |rho|>0.30
  paralogue_disentangling.csv          IKZF1-regulon-target-mean vs IKZF1/IKZF2/IKZF3
  ikzf1_target_hypergeometric.csv      enrichment vs JASPAR 2026 curated IKZF1 targets
  CONCLUSION.md
"""
import os
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
import loompy
from scipy import stats
from statsmodels.stats.multitest import multipletests

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
OUT = PROJECT / "results/phase7/grn_jaspar2026"
GRN = PROJECT / "results/phase2/GRN"
TRAJ = PROJECT / "results/phase1/trajectory"
AUC_LOOM = OUT / "scenic_auc_jaspar2026.loom"
REG_CSV = OUT / "regulons_jaspar2026.csv"
EXPR_LOOM = GRN / "microglia_raw.loom"
PT_CSV = TRAJ / "pseudotime.csv"
# Curated IKZF1 target ground truth (one HGNC symbol per line) built by
# code/slurm/paper2/s1_00_fetch_curated_targets.sh
# (ChEA3 Literature+ENCODE+ReMap ChIP-seq union + DoRothEA A/B/C via OmniPath).
CURATED_TARGETS = PROJECT / "data/references/curated_ikzf1_targets.txt"

RHO_THRESH = 0.30
N_PERM = 1000
SEED = 42
FOCUS_TFS = ["IKZF1", "BHLHE41", "BHLHE40", "IRF8", "SPI1", "RUNX1", "CEBPB", "PPARG"]


def load_auc():
    with loompy.connect(str(AUC_LOOM), mode="r", validate=False) as ds:
        a = ds.ca["RegulonsAUC"]
        regs = list(a.dtype.names)
        df = pd.DataFrame({r: a[r] for r in regs}, index=list(ds.ca["CellID"]))
        df["state"] = list(ds.ca["state"])
        df["donor"] = list(ds.ca["Donor_ID"]) if "Donor_ID" in ds.ca.keys() else np.nan
    return df, regs


def load_expr(genes):
    with loompy.connect(str(EXPR_LOOM), mode="r", validate=False) as ds:
        gsym = np.array(ds.ra["Gene"] if "Gene" in ds.ra.keys() else ds.ra[list(ds.ra.keys())[0]])
        cids = list(ds.ca["CellID"])
        idx = {g: int(np.where(gsym == g)[0][0]) for g in genes if g in gsym}
        mat = {g: ds[idx[g], :] for g in idx}
    return pd.DataFrame(mat, index=cids)


def regulon_targets(tf):
    df = pd.read_csv(REG_CSV, header=None, skiprows=3,
                     names=["TF", "MotifID", "AUC", "NES", "MotifSimQ", "OrthID",
                            "Annotation", "Context", "TargetGenes", "RankAtMax"])
    genes = set()
    for tg in df.loc[df["TF"] == tf, "TargetGenes"]:
        try:
            for item in ast.literal_eval(tg):
                genes.add(item[0])
        except Exception:
            pass
    return sorted(genes)


def main():
    rng = np.random.default_rng(SEED)
    auc, regs = load_auc()
    pt = pd.read_csv(PT_CSV, index_col="cell_id")["dpt_pseudotime"]
    m = auc.join(pt, how="inner")
    ptv = m["dpt_pseudotime"].to_numpy()
    donor = m["donor"].to_numpy()
    n = len(m)
    print(f"{n} cells, {len(regs)} rescued regulons")

    # ── (a,b,c) per-regulon pseudotime correlation + TF ranking ──────────────
    rows = []
    for r in regs:
        rho, p = stats.spearmanr(m[r].to_numpy(), ptv)
        rows.append({"regulon": r, "TF": r.split("(")[0], "rho": rho, "pval": p})
    bdf = pd.DataFrame(rows)
    bdf["padj"] = multipletests(bdf["pval"], method="fdr_bh")[1]
    state_mean = m.groupby("state")[regs].mean()
    bdf["peak_state"] = bdf["regulon"].map(lambda r: state_mean[r].idxmax())
    bdf = bdf.sort_values("rho", ascending=False)
    bdf.to_csv(OUT / "benchmark_regulon_pseudotime.csv", index=False)

    rank = bdf.reset_index(drop=True)
    rank["rho_rank"] = rank.index + 1
    tf_rank = rank[rank["TF"].isin(FOCUS_TFS)][["regulon", "TF", "rho", "padj",
                                               "peak_state", "rho_rank"]]
    tf_rank.to_csv(OUT / "benchmark_tf_ranking.csv", index=False)
    ik = bdf[bdf["TF"] == "IKZF1"]
    ikzf1_rho = float(ik["rho"].iloc[0]) if len(ik) else np.nan
    print(f"IKZF1 regulon rho = {ikzf1_rho:+.3f} (rank {tf_rank.loc[tf_rank.TF=='IKZF1','rho_rank'].tolist()})")

    # ── (d) permutation null: shuffle pseudotime, count regulons |rho|>0.30 ──
    obs_count = int((bdf["rho"].abs() >= RHO_THRESH).sum())
    obs_max = float(bdf["rho"].abs().max())
    auc_mat = m[regs].to_numpy()
    donors = pd.unique(donor)
    d_idx = {d: np.where(donor == d)[0] for d in donors}

    def perm_stats(pt_perm):
        rr = np.array([stats.spearmanr(auc_mat[:, j], pt_perm).statistic
                       for j in range(len(regs))])
        return int((np.abs(rr) >= RHO_THRESH).sum()), float(np.abs(rr).max())

    null_unrestr = np.empty((N_PERM, 2))
    null_block = np.empty((N_PERM, 2))
    for i in range(N_PERM):
        null_unrestr[i] = perm_stats(rng.permutation(ptv))
        # donor-block: permute whole-donor pseudotime blocks among donors
        perm = ptv.copy()
        shuf = rng.permutation(donors)
        for src, dst in zip(donors, shuf):
            take = d_idx[src]
            fill = d_idx[dst]
            k = min(len(take), len(fill))
            perm[fill[:k]] = ptv[take[:k]]
        null_block[i] = perm_stats(perm)

    fdr = pd.DataFrame([{
        "scheme": "unrestricted",
        "obs_n_regulons_ge_thresh": obs_count, "obs_max_abs_rho": obs_max,
        "null_mean_n": null_unrestr[:, 0].mean(),
        "null_p95_n": np.percentile(null_unrestr[:, 0], 95),
        "emp_p_maxrho_ge_ikzf1": (np.sum(null_unrestr[:, 1] >= abs(ikzf1_rho)) + 1) / (N_PERM + 1),
        "emp_fdr_at_thresh": null_unrestr[:, 0].mean() / max(obs_count, 1),
    }, {
        "scheme": "donor_block",
        "obs_n_regulons_ge_thresh": obs_count, "obs_max_abs_rho": obs_max,
        "null_mean_n": null_block[:, 0].mean(),
        "null_p95_n": np.percentile(null_block[:, 0], 95),
        "emp_p_maxrho_ge_ikzf1": (np.sum(null_block[:, 1] >= abs(ikzf1_rho)) + 1) / (N_PERM + 1),
        "emp_fdr_at_thresh": null_block[:, 0].mean() / max(obs_count, 1),
    }])
    fdr.to_csv(OUT / "null_permutation_fdr.csv", index=False)
    print(fdr.to_string(index=False))

    # ── (e) IKZF1 vs IKZF2 vs IKZF3 paralogue disentangling ─────────────────
    targets = regulon_targets("IKZF1")
    expr = load_expr(sorted(set(targets) | {"IKZF1", "IKZF2", "IKZF3"}))
    expr = expr.loc[expr.index.intersection(m.index)]
    tgt_present = [g for g in targets if g in expr.columns]
    tgt_mean = expr[tgt_present].mean(axis=1)
    par_rows = []
    for para in ["IKZF1", "IKZF2", "IKZF3"]:
        if para in expr.columns:
            rho, p = stats.spearmanr(tgt_mean, expr[para])
            par_rows.append({"paralogue": para, "rho_targetmean_vs_expr": rho,
                             "pval": p, "n_targets": len(tgt_present)})
    pd.DataFrame(par_rows).to_csv(OUT / "paralogue_disentangling.csv", index=False)
    print(pd.DataFrame(par_rows).to_string(index=False))

    # ── (f) hypergeometric: IKZF1 regulon targets vs JASPAR curated targets ─
    if CURATED_TARGETS.exists():
        curated = set(g.strip().upper() for g in CURATED_TARGETS.read_text().splitlines()
                      if g.strip())
        # universe = all genes GRNBoost2/pyscenic could have called (loom gene set)
        with loompy.connect(str(EXPR_LOOM), mode="r", validate=False) as ds:
            universe = set(str(x).upper() for x in ds.ra[list(ds.ra.keys())[0]])
        pred = set(t.upper() for t in targets)
        M = len(universe); n_c = len(curated & universe)
        N_p = len(pred & universe); k = len(pred & curated)
        p_hyper = stats.hypergeom.sf(k - 1, M, n_c, N_p)
        hg = {"source": "ChEA3 Lit/ENCODE/ReMap ChIP-seq + DoRothEA A/B/C (OmniPath)",
              "n_curated_total": len(curated), "universe": M,
              "n_curated_in_universe": n_c, "n_predicted_in_universe": N_p,
              "overlap": k,
              "fold_enrichment": (k / N_p) / (n_c / M) if N_p and n_c else np.nan,
              "hypergeom_p": p_hyper}
    else:
        hg = {"note": "run code/slurm/paper2/s1_00_fetch_curated_targets.sh first "
                      "(data/references/curated_ikzf1_targets.txt missing)"}
    pd.DataFrame([hg]).to_csv(OUT / "ikzf1_target_hypergeometric.csv", index=False)
    print(hg)

    # ── CONCLUSION ─────────────────────────────────────────────────────────
    (OUT / "CONCLUSION.md").write_text(f"""# WP1 — JASPAR 2026 GRN benchmark: CONCLUSION

- Rescued regulons: {len(regs)} (was 46 with cisTarget v10).
- IKZF1 regulon pseudotime rho = {ikzf1_rho:+.3f}; rank among all regulons =
  {tf_rank.loc[tf_rank.TF=='IKZF1','rho_rank'].tolist()}.
- Focus-TF ranking: see benchmark_tf_ranking.csv. Does BHLHE41/IRF8 now match
  or beat IKZF1 for trajectory coupling? -> {'BHLHE41 present' if (bdf.TF=='BHLHE41').any() else 'BHLHE41 still not retained'}.
- Permutation null ({N_PERM} perms): unrestricted empirical p(max|rho| >= IKZF1) =
  {fdr.loc[0,'emp_p_maxrho_ge_ikzf1']:.4g}; donor-block =
  {fdr.loc[1,'emp_p_maxrho_ge_ikzf1']:.4g}. Empirical FDR at |rho|>{RHO_THRESH}:
  unrestricted {fdr.loc[0,'emp_fdr_at_thresh']:.3f}, donor-block
  {fdr.loc[1,'emp_fdr_at_thresh']:.3f}.
- Paralogue test: IKZF1-regulon-target-mean correlates with IKZF1
  ({par_rows[0]['rho_targetmean_vs_expr']:+.3f}) vs IKZF2 / IKZF3
  ({[round(r['rho_targetmean_vs_expr'],3) for r in par_rows[1:]]}).
- Curated-target hypergeometric: {hg}

Manuscript: fold into Paper #2 Results (new section) + Table 1 evidence block.
""")
    print(f"\nWrote outputs + CONCLUSION.md to {OUT}")


if __name__ == "__main__":
    main()
