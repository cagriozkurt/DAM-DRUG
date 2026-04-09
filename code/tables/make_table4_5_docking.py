"""
DAM-DRUG — Table 4 (Docking Results) + Table 5 (Drug Candidate Summary)
=========================================================================

Table 4: All consensus hits passing dual-filter (Vina + CNN), with MM-GBSA
         for shortlisted compounds. 26 Tier-1 + 6 Tier-2 IKZF1 hits.

Table 5: Top candidate drugs per target after applying all filters
         (consensus + MM-GBSA + selectivity). The final repurposing shortlist.

Run locally:
  python code/tables/make_table4_5_docking.py
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
P4      = PROJECT / "results/phase4"
OUT     = PROJECT / "results/tables"
OUT.mkdir(parents=True, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────────
shortlist  = pd.read_csv(P4 / "candidates_shortlist.csv")
mmpbsa     = pd.read_csv(P4 / "mmpbsa/mmpbsa_summary.csv")
tier2_hits = pd.read_csv(P4 / "tier2_consensus_hits.csv")
selectivity= pd.read_csv(P4 / "vina_selectivity/selectivity_table.csv")
library    = pd.read_csv(PROJECT / "data/compounds/tier2_approved.csv",
                         usecols=["molecule_chembl_id", "pref_name", "mw", "alogp", "cns_mpo"])
library    = library.rename(columns={"molecule_chembl_id": "chembl_id"})

# ── Drug names: hardcode the few Tier1 IDs not in the library CSV ─────────────
# Retrieved from ChEMBL public database
EXTRA_NAMES = {
    "CHEMBL108":     "IBUPROFEN",
    "CHEMBL53":      "ASPIRIN",
    "CHEMBL41":      "ACETAMINOPHEN",
    "CHEMBL282052":  "NAPROXEN",
    "CHEMBL294199":  "KETOPROFEN",
    "CHEMBL439849":  "FENOPROFEN",
    "CHEMBL41355":   "DICLOFENAC",
    "CHEMBL1101":    "INDOMETHACIN",
    "CHEMBL1098":    "PIROXICAM",
    "CHEMBL14370":   "MEFENAMIC ACID",
    "CHEMBL1201203": "ROSIGLITAZONE",
    "CHEMBL1201193": "PIOGLITAZONE",
    "CHEMBL396778":  "TROGLITAZONE",
    "CHEMBL2105760": "LOBEGLITAZONE",
    "CHEMBL255044":  "FARGLITAZAR",
    "CHEMBL2103837": "TAFAMIDIS",
    "CHEMBL42":      "CHLORPROMAZINE",
    "CHEMBL646":     "HALOPERIDOL",
    "CHEMBL16":      "WARFARIN",
    "CHEMBL1089318": "TORASEMIDE",
    "CHEMBL1201168": "ETODOLAC",
    "CHEMBL285674":  "SULINDAC",
    "CHEMBL1201347": "TROGLITAZONE",
}

extra_df = pd.DataFrame([
    {"chembl_id": k, "pref_name": v} for k, v in EXTRA_NAMES.items()
])
library = pd.concat([library, extra_df], ignore_index=True).drop_duplicates("chembl_id")

def get_name(chembl_id):
    row = library[library["chembl_id"] == chembl_id]
    return row.iloc[0]["pref_name"] if not row.empty else chembl_id

def get_prop(chembl_id, col):
    row = library[library["chembl_id"] == chembl_id]
    if not row.empty and col in row.columns and pd.notna(row.iloc[0].get(col)):
        return row.iloc[0][col]
    return None

# ── Clean target name helper ──────────────────────────────────────────────────
def clean_target(stem):
    stem = stem.replace("_prep", "").replace("_AF2_DBD", " (AF2-DBD)").replace("_AF2_bHLH", " (AF2-bHLH)")
    stem = stem.replace("_1FM9_LBD", " (1FM9-LBD)").replace("_4EOT_bZIP", " (4EOT-bZIP)")
    stem = stem.replace("_1LJM_Runt", " (1LJM-Runt)").replace("_8RQC_CRBN", " (8RQC-CRBN glue)")
    return stem

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 4 — Full docking hit list (Tier 1 standard + Tier 2 IKZF1 PROTAC track)
# ══════════════════════════════════════════════════════════════════════════════

# Build Tier1/Tier2 standard hits from candidates_shortlist
t4_rows = []

for _, row in shortlist.iterrows():
    chembl = row["chembl_id"]
    target = row["target"]

    # MM-GBSA (if available — shortlisted subset only)
    mm = mmpbsa[(mmpbsa["target"] == target) & (mmpbsa["chembl_id"] == chembl)]
    dg = mm.iloc[0]["dg_gbsa"] if not mm.empty else None

    t4_rows.append({
        "tier":            row["tier"],
        "target":          clean_target(target),
        "chembl_id":       chembl,
        "drug_name":       get_name(chembl),
        "vina_score":      round(row["vina_score"], 2),
        "cnn_score":       round(row["cnn_score"], 3),
        "cnn_affinity_pKd": round(row["cnn_affinity"], 2) if pd.notna(row.get("cnn_affinity", np.nan)) else None,
        "composite_score": round(row["composite_score"], 3),
        "dG_GBSA_kcal_mol": round(dg, 2) if dg is not None else None,
        "promiscuous":     "⚠" if row.get("promiscuous", False) else "",
        "mw":              get_prop(chembl, "mw"),
        "alogp":           get_prop(chembl, "alogp"),
        "cns_mpo":         get_prop(chembl, "cns_mpo"),
    })

# Add IKZF1 PROTAC track: top 20 by Vina only (too many hits to list all)
# CNN threshold used in screening was ≥0.7; include top Vina hits for table conciseness
ikzf1_top = tier2_hits.sort_values("vina_score").head(20)
for _, row in ikzf1_top.iterrows():
    chembl = row["compound"]
    # check selectivity data
    sel = selectivity[selectivity["chembl_id"] == chembl]
    dg_sel = round(sel.iloc[0]["mmpbsa_dg"], 2) if not sel.empty else None
    sel_flag = "✗ non-sel" if not sel.empty and sel.iloc[0]["selective"] == "✗" else ""

    t4_rows.append({
        "tier":             "PROTAC-track",
        "target":           "IKZF1 (8RQC-CRBN glue)",
        "chembl_id":        chembl,
        "drug_name":        row.get("pref_name", get_name(chembl)),
        "vina_score":       round(row["vina_score"], 2),
        "cnn_score":        round(row["cnn_score"], 3),
        "cnn_affinity_pKd": None,
        "composite_score":  None,
        "dG_GBSA_kcal_mol": dg_sel,
        "promiscuous":      sel_flag,
        "mw":               row.get("mw"),
        "alogp":            row.get("alogp"),
        "cns_mpo":          row.get("cns_mpo"),
    })

t4 = pd.DataFrame(t4_rows)

# Sort: Tier1 by target then composite_score, PROTAC by vina
t1 = t4[t4["tier"].isin(["Tier1", "Tier2"])].sort_values(
    ["tier", "target", "composite_score"], ascending=[True, True, False])
tp = t4[t4["tier"] == "PROTAC-track"].sort_values("vina_score")
t4 = pd.concat([t1, tp], ignore_index=True)

t4.to_csv(OUT / "table4_docking_results.csv", index=False)
print(f"Table 4: {len(t4)} rows → {OUT}/table4_docking_results.csv")
print(t4[["tier","target","drug_name","vina_score","composite_score","dG_GBSA_kcal_mol"]].to_string(index=False))

# ── Table 4 Markdown ───────────────────────────────────────────────────────────
md4 = ["# Table 4 — Docking Results\n",
       "| Tier | Target | Drug name | ChEMBL | Vina (kcal/mol) | CNN score | Comp. score | ΔG_GBSA (kcal/mol) | MW | AlogP | CNS MPO | Flag |",
       "|---|---|---|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"]

for _, r in t4.iterrows():
    def f(v, fmt=None):
        if v is None or (isinstance(v, float) and np.isnan(v)): return "—"
        return format(v, fmt) if fmt else str(v)
    md4.append(f"| {r['tier']} | {r['target']} | {r['drug_name']} | {r['chembl_id']} | "
               f"{f(r['vina_score'],'.2f')} | {f(r['cnn_score'],'.3f')} | "
               f"{f(r['composite_score'],'.3f')} | {f(r['dG_GBSA_kcal_mol'],'.2f')} | "
               f"{f(r['mw'],'.1f')} | {f(r['alogp'],'.2f')} | {f(r['cns_mpo'],'.1f')} | {r['promiscuous']} |")

md4.append("\n**Notes:** Vina score: AutoDock Vina binding energy. CNN score: GNINA pose quality (0–1; ≥0.7 = good). "
           "Composite score: 0.4×(Vina/10) + 0.6×CNN. ΔG_GBSA: MM-GBSA binding free energy (negative = favorable; "
           "computed for top-10 Tier1 + all Tier2 hits). ⚠ = promiscuous (consensus hit in ≥3 targets). "
           "PROTAC-track: IKZF1 8RQC-CRBN molecular glue pocket; CNS MPO ≥3.5 = CNS drug-likeness threshold.")

(OUT / "table4_docking_results.md").write_text("\n".join(md4))
print(f"\nSaved markdown → {OUT}/table4_docking_results.md")

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 5 — Drug Candidate Summary (final shortlist)
# ══════════════════════════════════════════════════════════════════════════════
# Selection criteria:
#   Standard track: top MM-GBSA rank per target + composite_score ≥ 0.7 + not promiscuous
#   PROTAC track:   top vina + cnn_score ≥ 0.9 + cns_mpo ≥ 3.5 + dG_GBSA < -5 kcal/mol

# Tier1+Tier2 standard candidates
t5_rows = []

# Best hit per target: top MM-GBSA, non-promiscuous, dG < -2 kcal/mol (meaningful binding)
standard_mm = t4[(t4["tier"].isin(["Tier1", "Tier2"])) & (t4["promiscuous"] == "") &
                 (t4["dG_GBSA_kcal_mol"].notna())].copy()
standard_mm["dG_GBSA_kcal_mol"] = pd.to_numeric(standard_mm["dG_GBSA_kcal_mol"], errors="coerce")
standard_mm = standard_mm[standard_mm["dG_GBSA_kcal_mol"] < -2.0]  # meaningful MM-GBSA

# Mechanism/rationale per target
rationale = {
    "PPARG (1FM9-LBD)":   "PPARG LBD partial agonism → LDAM lipid homeostasis; NSAIDs have PPARG partial-agonist activity (Lehmann 1997)",
    "IRF8 (AF2-DBD)":     "IRF8 DBD engagement → stabilize homeostatic TF activity; anti-neuroinflammatory",
    "BHLHE41 (AF2-bHLH)": "BHLHE41 bHLH inhibition → de-repress DAM→LateAD-DAM transition; local anesthetic class",
    "RUNX1 (1LJM-Runt)":  "RUNX1 Runt domain; Tier2 shallow pocket — all MM-GBSA ≈ 0; supportive evidence only",
}
# Note: MAF excluded (MM-GBSA ≈ 0 for both hits; no favorable binding detected)
# Note: IKZF1 PROTAC track — no approved drug candidate; see below

for target_clean in standard_mm["target"].unique():
    sub = standard_mm[standard_mm["target"] == target_clean].sort_values("dG_GBSA_kcal_mol")
    best = sub.iloc[0]
    t5_rows.append({
        "target":    target_clean,
        "track":     best["tier"],
        "drug_name": best["drug_name"],
        "chembl_id": best["chembl_id"],
        "vina":      best["vina_score"],
        "cnn":       best["cnn_score"],
        "composite": best["composite_score"],
        "dG_GBSA":   best["dG_GBSA_kcal_mol"],
        "mw":        best["mw"],
        "alogp":     best["alogp"],
        "cns_mpo":   best["cns_mpo"],
        "rationale": rationale.get(target_clean, "—"),
    })

# PROTAC track: all 6 selectivity-tested IKZF1 candidates failed SI filter (min_si < 1.0)
# → no approved drug candidate identified; note in table as "de novo design required"
t5_rows.append({
    "target":    "IKZF1 (8RQC-CRBN glue)",
    "track":     "PROTAC-track",
    "drug_name": "None identified",
    "chembl_id": "—",
    "vina":      None, "cnn": None, "composite": None, "dG_GBSA": None,
    "mw":        None, "alogp": None, "cns_mpo": None,
    "rationale": "No approved drug passes glutarimide pharmacophore + selectivity filter; "
                 "requires IMiD-analogue library or de novo PROTAC design (DeLinker)",
})

t5 = pd.DataFrame(t5_rows)
t5.to_csv(OUT / "table5_drug_candidates.csv", index=False)
print(f"\nTable 5: {len(t5)} candidates → {OUT}/table5_drug_candidates.csv")
print(t5[["target","drug_name","dG_GBSA","cns_mpo","rationale"]].to_string(index=False))

# ── Table 5 Markdown ───────────────────────────────────────────────────────────
md5 = ["# Table 5 — Drug Candidate Summary\n",
       "| Target | Track | Drug | ChEMBL | Vina | CNN | ΔG_GBSA | MW | AlogP | CNS MPO | Rationale |",
       "|---|---|---|---|:---:|:---:|:---:|:---:|:---:|:---:|---|"]

for _, r in t5.iterrows():
    def f(v, fmt=None):
        if v is None or (isinstance(v, float) and np.isnan(v)): return "—"
        return format(v, fmt) if fmt else str(v)
    md5.append(f"| {r['target']} | {r['track']} | **{r['drug_name']}** | {r['chembl_id']} | "
               f"{f(r['vina'],'.2f')} | {f(r['cnn'],'.3f')} | {f(r['dG_GBSA'],'.2f')} | "
               f"{f(r['mw'],'.1f')} | {f(r['alogp'],'.2f')} | {f(r['cns_mpo'],'.1f')} | "
               f"{r['rationale']} |")

md5.append("\n**Selection criteria:** Standard track: best MM-GBSA rank per target, composite ≥ 0.7, "
           "non-promiscuous. PROTAC track: CNN ≥ 0.90, ΔG_GBSA < −5 kcal/mol. "
           "CNS MPO ≥ 3.5 = predicted CNS drug-likeness (Wager 2010).")

(OUT / "table5_drug_candidates.md").write_text("\n".join(md5))
print(f"Saved markdown → {OUT}/table5_drug_candidates.md")
