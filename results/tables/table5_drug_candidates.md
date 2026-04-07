# Table 5 — Drug Candidate Summary

| Target | Track | Drug | ChEMBL | Vina | CNN | ΔG_GBSA | MW | AlogP | CNS MPO | Rationale |
|---|---|---|---|:---:|:---:|:---:|:---:|:---:|:---:|---|
| IRF8 (AF2-DBD) | Tier1 | **TAFAMIDIS** | CHEMBL2103837 | -6.71 | 0.741 | -9.50 | 308.1 | 4.50 | 4.0 | IRF8 DBD engagement → stabilize homeostatic TF activity; anti-neuroinflammatory |
| PPARG (1FM9-LBD) | Tier1 | **DIFLUNISAL** | CHEMBL898 | -8.77 | 0.955 | -2.80 | 250.2 | 3.04 | 4.0 | PPARG LBD partial agonism → LDAM lipid homeostasis; NSAIDs have PPARG partial-agonist activity (Lehmann 1997) |
| BHLHE41 (AF2-bHLH) | Tier2 | **DIBUCAINE** | CHEMBL1086 | -6.00 | 0.805 | -21.89 | 343.5 | 3.49 | 3.5 | BHLHE41 bHLH inhibition → de-repress DAM→LateAD-DAM transition; local anesthetic class |
| IKZF1 (8RQC-CRBN glue) | PROTAC-track | **None identified** | — | — | — | — | — | — | — | No approved drug passes glutarimide pharmacophore + selectivity filter; requires IMiD-analogue library or de novo PROTAC design (DeLinker) |

**Selection criteria:** Standard track: best MM-GBSA rank per target, composite ≥ 0.7, non-promiscuous. PROTAC track: CNN ≥ 0.90, ΔG_GBSA < −5 kcal/mol. CNS MPO ≥ 3.5 = predicted CNS drug-likeness (Wager 2010).