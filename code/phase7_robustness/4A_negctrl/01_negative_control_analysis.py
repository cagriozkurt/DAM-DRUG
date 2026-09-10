"""
Section 4A — Tafamidis & diflunisal as MM-GBSA negative controls
==============================================================
Assembles one scorecard from EXISTING phase-3/4 outputs (no new computation)
showing that implicit-solvent MM-GBSA assigned favourable dG / top rank to two
ligands that then leave the pocket in explicit-solvent 100 ns MD — an
implicit-solvent + shallow-pocket failure mode.

Inputs (all present):
  results/phase4/mmpbsa/mmpbsa_summary.csv
  results/phase4/vina_scores_*.csv, results/phase4/gnina_scores_*.csv
  results/phase4/consensus_hits.csv
  results/phase3/druggability_ranking.csv
MD core-RMSD values are from the manuscript / Supp Fig S8 (archived TRUBA
trajectories; gen_seed = -1, not code-reproducible).

Outputs:
  results/phase7/negctrl/negative_control_scorecard.csv
  results/phase7/negctrl/mmgbsa_vs_md_retention.png
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
P4 = PROJECT / "results/phase4"
P3 = PROJECT / "results/phase3"
OUT = PROJECT / "results/phase7/negctrl"
OUT.mkdir(parents=True, exist_ok=True)

# CHEMBL id -> (name, on-target label used in phase 4)
CONTROLS = {
    "CHEMBL2103837": ("Tafamidis", "IRF8_AF2_DBD_prep"),
    "CHEMBL898":     ("Diflunisal", "PPARG_1FM9_LBD_prep"),
}
# MD core-RMSD (last 20 ns), from manuscript lines 481-483 / Supp Fig S8
MD_RMSD = {
    "CHEMBL2103837": {"target": "IRF8_AF2_DBD_prep", "core_rmsd_A": 25.17, "sd": 3.56,
                      "interpretation": "ligand egress"},
    "CHEMBL898":     {"target": "PPARG_1FM9_LBD_prep", "core_rmsd_A": 83.42, "sd": 2.89,
                      "interpretation": "complete egress"},
}
# other compounds with 100 ns MD trajectories (md100ns_manifest.txt) — context
MD_CONTEXT = ["CHEMBL42", "CHEMBL490"]   # IRF8; retained (no egress reported)


def load_gnina(target):
    f = P4 / f"gnina_scores_{target}.csv"
    df = pd.read_csv(f, names=["chembl_id", "vina", "cnn", "cnn_aff"], header=0)
    return df.set_index("chembl_id")


def main():
    mm = pd.read_csv(P4 / "mmpbsa" / "mmpbsa_summary.csv")
    drug = pd.read_csv(P3 / "druggability_ranking.csv").set_index("structure_id")

    rows = []
    for cid, (name, tgt) in CONTROLS.items():
        sub = mm[mm["target"] == tgt].sort_values("dg_gbsa").reset_index(drop=True)
        r = sub[sub["chembl_id"] == cid]
        dg = float(r["dg_gbsa"].iloc[0]) if len(r) else np.nan
        sem = float(r["sem"].iloc[0]) if len(r) else np.nan
        rank_gbsa = int(r["rank_gbsa"].iloc[0]) if len(r) else np.nan
        n_in_set = len(sub)
        gn = load_gnina(tgt)
        vina = float(gn.loc[cid, "vina"]) if cid in gn.index else np.nan
        cnn = float(gn.loc[cid, "cnn"]) if cid in gn.index else np.nan
        pk = drug.loc[tgt] if tgt in drug.index else None
        rows.append({
            "chembl_id": cid, "compound": name, "on_target": tgt,
            "vina_kcal_mol": vina, "gnina_cnn": cnn,
            "mmgbsa_dg_kcal_mol": dg, "mmgbsa_sem": sem,
            "mmgbsa_rank_in_target_set": f"{rank_gbsa}/{n_in_set}",
            "md_core_rmsd_A": MD_RMSD[cid]["core_rmsd_A"],
            "md_core_rmsd_sd": MD_RMSD[cid]["sd"],
            "md_outcome": MD_RMSD[cid]["interpretation"],
            "pocket_volume_A3": float(pk["volume_A3"]) if pk is not None else np.nan,
            "pocket_hydrophobicity": float(pk["hydrophobicity"]) if pk is not None else np.nan,
            "pocket_drug_score": float(pk["drug_score"]) if pk is not None else np.nan,
            "pocket_plddt_flag": pk["af2_pocket_plddt_flag"] if pk is not None else "",
        })

    # context compounds (IRF8 MD, retained)
    gn_irf8 = load_gnina("IRF8_AF2_DBD_prep")
    mm_irf8 = mm[mm["target"] == "IRF8_AF2_DBD_prep"].sort_values("dg_gbsa").reset_index(drop=True)
    for cid in MD_CONTEXT:
        r = mm_irf8[mm_irf8["chembl_id"] == cid]
        rows.append({
            "chembl_id": cid, "compound": f"(MD context: {cid})",
            "on_target": "IRF8_AF2_DBD_prep",
            "vina_kcal_mol": float(gn_irf8.loc[cid, "vina"]) if cid in gn_irf8.index else np.nan,
            "gnina_cnn": float(gn_irf8.loc[cid, "cnn"]) if cid in gn_irf8.index else np.nan,
            "mmgbsa_dg_kcal_mol": float(r["dg_gbsa"].iloc[0]) if len(r) else np.nan,
            "mmgbsa_sem": float(r["sem"].iloc[0]) if len(r) else np.nan,
            "mmgbsa_rank_in_target_set": f"{int(r['rank_gbsa'].iloc[0])}/{len(mm_irf8)}" if len(r) else "NA",
            "md_core_rmsd_A": np.nan, "md_core_rmsd_sd": np.nan,
            "md_outcome": "retained (no egress reported)",
            "pocket_volume_A3": np.nan, "pocket_hydrophobicity": np.nan,
            "pocket_drug_score": np.nan, "pocket_plddt_flag": "",
        })

    sc = pd.DataFrame(rows)
    sc.to_csv(OUT / "negative_control_scorecard.csv", index=False)
    print(sc.to_string(index=False))

    # plot: MM-GBSA dG (x) vs MD outcome (y) for the compounds with MD
    md_pts = sc[sc["md_core_rmsd_A"].notna() | sc["md_outcome"].str.contains("retained")]
    fig, ax = plt.subplots(figsize=(7, 4))
    for _, r in sc.iterrows():
        if pd.isna(r["mmgbsa_dg_kcal_mol"]):
            continue
        y = r["md_core_rmsd_A"] if pd.notna(r["md_core_rmsd_A"]) else 2.0
        egress = pd.notna(r["md_core_rmsd_A"]) and r["md_core_rmsd_A"] > 5
        ax.scatter(r["mmgbsa_dg_kcal_mol"], y,
                   color="#C44E52" if egress else "#4C72B0", s=70, zorder=3)
        ax.annotate(r["compound"].replace("(MD context: ", "").rstrip(")"),
                    (r["mmgbsa_dg_kcal_mol"], y), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(5, color="grey", ls="--", lw=1, label="egress threshold (core-RMSD 5 Å)")
    ax.set_xlabel("MM-GBSA ΔG (kcal/mol) — more negative = 'better'")
    ax.set_ylabel("100 ns MD core-RMSD (Å, last 20 ns)")
    ax.set_title("Implicit-solvent MM-GBSA ΔG does not predict MD pocket retention")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "mmgbsa_vs_md_retention.png", dpi=200)
    print(f"\nWrote {OUT}/negative_control_scorecard.csv and mmgbsa_vs_md_retention.png")


if __name__ == "__main__":
    main()
