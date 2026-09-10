"""
WP2.4 — Collate multi-region CellChat + expression-matched permutation null
=========================================================================
Determines whether SLIT2 -> ROBO2 (interneuron -> microglia) is a pan-cortical
axis or MTG-specific, and whether the observed MTG interaction probability
(sum_prob = 0.318 in the accepted paper) exceeds an expression-matched
background of 1,000 random ligand-receptor pairs.

Inputs:
  results/phase7/cellchat_multiregion/<R>/lr_interactions_<R>.csv   (all 9 regions)
  results/phase2/LR/  (MTG reference tables from the accepted paper, if present)
  per-region prep/counts + cell_meta for the expression-matched null

Outputs -> results/phase7/cellchat_multiregion/
  slit2_robo2_by_region.csv
  permutation_null_by_region.csv
  CONCLUSION.md
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
CC = PROJECT / "results/phase7/cellchat_multiregion"
REGIONS = ["AnG", "DFC", "FI", "HIP", "ITG", "LEC", "MEC", "STG", "V1C"]
LIGAND, RECEPTOR = "SLIT2", "ROBO2"
SENDER, RECEIVER = "InhNeuron", "Microglia"
N_PERM = 1000
SEED = 42
MTG_SUM_PROB = 0.318   # accepted-paper value for reference


def load_counts(region):
    p = CC / region / "prep" / "counts_raw.h5"
    with h5py.File(p, "r") as f:
        shape = tuple(f.attrs["shape"])          # (ncells, ngenes)
        from scipy.sparse import csc_matrix
        M = csc_matrix((f["data"][:], f["indices"][:], f["indptr"][:]), shape=shape)
        genes = [g.decode() for g in f["gene_names"][:]]
        barcodes = [b.decode() for b in f["barcodes"][:]]
    meta = pd.read_csv(CC / region / "prep" / "cell_meta.csv", index_col=0)
    meta = meta.loc[barcodes]
    return M.tocsc(), np.array(genes), meta


def mean_expr_by_group(M, genes, meta, group_col="cell_type_broad"):
    # M: cells x genes CSC
    out = {}
    for grp, idx in meta.groupby(group_col).groups.items():
        rows = [meta.index.get_loc(b) for b in idx]
        out[grp] = np.asarray(M[rows, :].mean(axis=0)).ravel()
    return pd.DataFrame(out, index=genes)


def main():
    rng = np.random.default_rng(SEED)

    # ── 1. SLIT2->ROBO2 rank + prob per region ─────────────────────────────
    rows = []
    for R in REGIONS:
        f = CC / R / f"lr_interactions_{R}.csv"
        if not f.exists():
            rows.append({"region": R, "status": "missing"})
            continue
        lr = pd.read_csv(f)
        lr = lr.sort_values("prob", ascending=False).reset_index(drop=True)
        lr["overall_rank"] = lr.index + 1
        hit = lr[(lr["ligand"] == LIGAND) & (lr["receptor"] == RECEPTOR)]
        specific = hit[(hit["source"] == SENDER) & (hit["target"] == RECEIVER)]
        rows.append({
            "region": R,
            "slit2_robo2_detected": len(hit) > 0,
            "sum_prob_all_pairs": float(hit["prob"].sum()) if len(hit) else 0.0,
            "best_source": hit["source"].iloc[0] if len(hit) else None,
            "best_target": hit["target"].iloc[0] if len(hit) else None,
            "best_prob": float(hit["prob"].iloc[0]) if len(hit) else 0.0,
            "best_overall_rank": int(hit["overall_rank"].iloc[0]) if len(hit) else None,
            "inh_to_micro_prob": float(specific["prob"].iloc[0]) if len(specific) else 0.0,
            "n_total_interactions": len(lr),
        })
    by_region = pd.DataFrame(rows)
    by_region.to_csv(CC / "slit2_robo2_by_region.csv", index=False)
    print(by_region.to_string(index=False))

    # ── 2. expression-matched permutation null per region ──────────────────
    perm_rows = []
    for R in REGIONS:
        f = CC / R / f"lr_interactions_{R}.csv"
        cph = CC / R / "prep" / "counts_raw.h5"
        if not f.exists() or not cph.exists():
            continue
        lr = pd.read_csv(f)
        obs = lr[(lr["ligand"] == LIGAND) & (lr["receptor"] == RECEPTOR)]["prob"].sum()

        M, genes, meta = load_counts(R)
        mg = mean_expr_by_group(M, genes, meta)
        if SENDER not in mg.columns or RECEIVER not in mg.columns:
            perm_rows.append({"region": R, "status": f"missing {SENDER}/{RECEIVER}"})
            continue
        gidx = {g: i for i, g in enumerate(genes)}
        if LIGAND not in gidx or RECEPTOR not in gidx:
            perm_rows.append({"region": R, "status": "SLIT2/ROBO2 not in matrix"})
            continue

        lig_lvl = mg.loc[LIGAND, SENDER]
        rec_lvl = mg.loc[RECEPTOR, RECEIVER]
        # candidate ligands/receptors: genes expressed in sender/receiver within
        # +/- 0.5 log-units of the observed levels (expression-matched)
        sender_expr = mg[SENDER]
        recv_expr = mg[RECEIVER]
        cand_lig = sender_expr.index[(sender_expr > 0) &
                                     (np.abs(np.log1p(sender_expr) - np.log1p(lig_lvl)) < 0.5)]
        cand_rec = recv_expr.index[(recv_expr > 0) &
                                   (np.abs(np.log1p(recv_expr) - np.log1p(rec_lvl)) < 0.5)]
        cand_lig = np.array(cand_lig)
        cand_rec = np.array(cand_rec)
        if len(cand_lig) < 20 or len(cand_rec) < 20:
            perm_rows.append({"region": R, "status": "too few expression-matched genes"})
            continue

        # null: random L-R pairs scored by the CellChat probability proxy
        # prob ~ (L in sender) * (R in receiver) / (Kh + product)  [triMean-agnostic proxy]
        def proxy(l, r):
            a = mg.loc[l, SENDER]; b = mg.loc[r, RECEIVER]
            return (a * b) / (0.5 + a * b)
        obs_proxy = proxy(LIGAND, RECEPTOR)
        null = np.array([proxy(rng.choice(cand_lig), rng.choice(cand_rec))
                         for _ in range(N_PERM)])
        p_emp = (np.sum(null >= obs_proxy) + 1) / (N_PERM + 1)
        perm_rows.append({
            "region": R, "obs_sum_prob": float(obs),
            "obs_proxy": float(obs_proxy), "null_mean_proxy": float(null.mean()),
            "null_p95_proxy": float(np.percentile(null, 95)),
            "perm_p": float(p_emp),
            "n_cand_ligands": int(len(cand_lig)), "n_cand_receptors": int(len(cand_rec)),
        })
    perm = pd.DataFrame(perm_rows)
    perm.to_csv(CC / "permutation_null_by_region.csv", index=False)
    print(perm.to_string(index=False))

    # ── 3. verdict ────────────────────────────────────────────────────────
    detected = by_region[by_region.get("slit2_robo2_detected", False) == True]
    n_regions_detected = len(detected)
    n_inh_micro = int((by_region.get("inh_to_micro_prob", pd.Series(dtype=float)) > 0).sum())
    sig = perm[perm.get("perm_p", 1) < 0.05] if "perm_p" in perm.columns else pd.DataFrame()
    verdict = ("PAN-CORTICAL" if n_inh_micro >= 5 and len(sig) >= 5
               else "PARTIAL" if n_inh_micro >= 2
               else "MTG-RESTRICTED")
    (CC / "CONCLUSION.md").write_text(f"""# WP2 — SLIT2->ROBO2 multi-region triage: CONCLUSION

**Verdict: {verdict}**

- SLIT2->ROBO2 detected in {n_regions_detected}/9 non-MTG regions;
  InhNeuron->Microglia specifically in {n_inh_micro}/9.
- Reference (accepted paper, MTG): sum_prob = {MTG_SUM_PROB}.
- Expression-matched permutation null: significant (p<0.05) in {len(sig)}/9 regions.

See slit2_robo2_by_region.csv and permutation_null_by_region.csv.

## Manuscript action
- If PAN-CORTICAL: SLIT2->ROBO2 becomes a validated pan-cortical axis in Paper #2.
- If MTG-RESTRICTED / PARTIAL: downgrade to an exploratory, region-restricted
  MTG hypothesis in the abstract and discussion (TODO.md Section 3 item 4).
""")
    print(f"\nVerdict: {verdict}. Wrote CONCLUSION.md")


if __name__ == "__main__":
    main()
