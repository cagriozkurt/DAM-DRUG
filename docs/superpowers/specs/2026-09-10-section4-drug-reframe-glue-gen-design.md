# Section 4 — Drug Pipeline Reframe + CRBN Molecular-Glue Generation: Design Spec

**Date:** 2026-09-10
**Scope:** `TODO.md` Section 4, the two LOCAL sub-tasks only —
(4A) reframe tafamidis/diflunisal as MM-GBSA negative controls, and
(4B) RDKit/BRICS molecular-glue generation toward the IKZF1 ZF2 degron.
**Explicitly out of scope:** explicit-solvent MD, T-REMD, hERG QSAR,
counter-docking, in-vitro gate definitions (TODO §4 items 4–6 — deferred /
partly TRUBA).
**Compute:** local. RDKit work runs in the existing `lipogate` env
(rdkit 2026.03.1); no new MD.
**Branch:** continue on `robustness/section2-decoupling` (or a fresh
`robustness/section4-drug` — decide at approval).

---

## Context (verified)

- Manuscript already softened these to "low-confidence computational
  hypotheses" (lines 464, 659–660) but still calls them "top-ranked" and lists
  them as candidates in Tables 3/4/5 and Fig 5A (lines 479–483, 498, 522, 1186).
- **Tafamidis = CHEMBL2103837 → IRF8**: MM-GBSA ΔG = −9.5 kcal/mol (rank #1 of
  the IRF8 MM-GBSA set), 100 ns MD core-RMSD = 25.17 ± 3.56 Å (last 20 ns) —
  ligand egress.
- **Diflunisal = CHEMBL898 → PPARG**: ΔG = −2.8 kcal/mol, MD core-RMSD =
  83.42 ± 2.89 Å — complete egress (Supp Fig S8).
- Druggability (`results/phase3/druggability_ranking.csv`): IRF8 AF2-DBD pocket
  drug_score 0.659 but volume only 500 Å³, hydrophobicity 11.8, `LOW` pLDDT
  flag; IKZF1 orthosteric ZF pocket drug_score **0.001** (undruggable).
- CRBN redirect target: **PDB 8RQC**. `data/structures/pdb/8RQC.pdb` has the
  CELMoD ligand **QFC** (~42 heavy atoms); `8RQC_CRBN_ZF2.pdb` = CRBN (chains
  A/D, res 64–425) + IKZF1 **ZF2 degron** (chains B/E, res 144–170) + 4 Zn,
  ligand stripped.
- Existing CNS-MPO method (`code/slurm/26_fetch_tier2.py:133`): a 5-check
  proxy (no pKa), each 0/0.5/1 — MW ≤360, cLogP ≤3, TPSA 40–90, HBD ==0,
  RTB ≤8; max 5; Tier-2 cut was ≥3.0.
- MD `.xvg` time series are **not** local (only `md100ns_manifest.txt`); exact
  RMSD values come from the manuscript / Supp table / archived TRUBA trajectories.

---

## Task 4A — Tafamidis & diflunisal as MM-GBSA negative controls

**Directory:** `code/phase7_robustness/4A_negctrl/`

### 4A.1 `01_negative_control_analysis.py`
Assemble one scorecard from existing outputs (no new computation):

| field | source |
|---|---|
| MM-GBSA ΔG, SEM, rank within target | `results/phase4/mmpbsa/mmpbsa_summary.csv` |
| Vina score, CNN, composite | `results/phase4/vina_scores_*`, `gnina_scores_*` |
| Selectivity index (SI) | `results/phase4/vina_selectivity/` |
| MD core-RMSD (egress) | manuscript / Supp Fig S8 values, cited as archived-trajectory-derived |
| target pocket volume, hydrophobicity, polarity, pLDDT flag | `results/phase3/druggability_ranking.csv` |

Output: `results/phase7/negctrl/negative_control_scorecard.csv`
(rows: tafamidis, diflunisal; + the 4 compounds with MD as context).

### 4A.2 Discordance framing
- Rank the MM-GBSA set by ΔG; overlay MD retention (the 4 compounds with
  trajectories: IRF8 CHEMBL2103837 / CHEMBL42 / CHEMBL490, PPARG CHEMBL898).
  Show MM-GBSA ΔG rank does **not** predict MD pocket retention (n=4,
  qualitative — stated as such).
- Root-cause narrative: implicit-solvent MM-GBSA overstabilises ligands on
  shallow, solvent-exposed TF surface pockets (IRF8 DBD 500 Å³, low
  hydrophobicity, AF2 `LOW` confidence; PPARG hit not in the canonical LBD
  sub-site). Both "hits" leave the pocket within 100 ns of explicit-solvent MD.

### 4A.3 Output — `results/phase7/negctrl/CONCLUSION.md`
Manuscript patch:
- Retract "top-ranked" / candidate framing for tafamidis & diflunisal in
  Abstract, Results, Tables 3–5, Fig 5A legend.
- Relabel them **"MM-GBSA sensitivity benchmarks"** / negative controls that
  expose the implicit-solvent + shallow-pocket failure mode.
- Keep the honest existing caveats; make the reframe explicit and consistent.
- The CRBN-track hits (celecoxib, teriflunomide, …) are a separate matter —
  not touched here beyond noting they remain exploratory.

---

## Task 4B — CRBN molecular-glue generation (RDKit / BRICS)

**Directory:** `code/phase7_robustness/4B_glue_gen/`  •  **env:** `lipogate`

### 4B.1 `01_extract_anchor.py`
- Extract ligand **QFC** from `data/structures/pdb/8RQC.pdb` → reference SMILES
  (via RDKit from the coordinates + bond perception, cross-checked against the
  PDB chemical component `QFC` if fetchable; fallback = manual SMILES).
- Define the canonical CRBN warhead: lenalidomide / pomalidomide
  **isoindolinone–glutarimide** (3-(1-oxoisoindolin-2-yl)piperidine-2,6-dione),
  SMILES `O=C1CCC(N2Cc3ccccc3C2=O)C(=O)N1`.
- Mark the **growth exit vector** at the isoindolinone 4- (or 5-) aryl position
  (where CELMoD linkers attach) with a dummy atom → anchor fragment
  `O=C1CCC(N2Cc3cccc([*])c3C2=O)C(=O)N1`.
- Save `results/phase7/glue_design/anchor.json` (warhead SMILES, exit vector,
  QFC reference, tri-Trp cage residues His378/Trp380/Trp386/Trp400 for docking).

### 4B.2 `02_build_fragment_pool.py`
- BRICS-decompose `data/compounds/tier1_cns_approved.csv` (285 CNS-approved) +
  `data/compounds/tier2_approved.csv` (1,677) → fragment multiset.
- Keep fragments: 1 attachment point OR 2 (linkers), MW 40–250, no PAINS,
  no reactive/unstable groups (acyl halides, epoxides, Michael acceptors via
  SMARTS), ≤2 rings. Dedupe by canonical SMILES; keep frequency.
- Output `results/phase7/glue_design/fragment_pool.csv` (smiles, n_attach,
  mw, source_count).

### 4B.3 `03_generate_glues.py`
- Growth = attach 1–2 pool fragments at the anchor exit vector via
  RDKit `BRICS.BRICSBuild` seeded with `[anchor fragment] + pool`
  (and a Reaction-SMARTS fallback for a controlled 1-fragment grow).
- Enumerate up to 20,000 products; keep only those whose Murcko/substructure
  match retains the **intact glutarimide + isoindolinone** anchor
  (`HasSubstructMatch`).
- Sanitize, dedupe (InChIKey), ETKDGv3 embed 1 conformer, MMFF optimise.
- Output `results/phase7/glue_design/generated_raw.csv` + `.sdf`.
- Seed `random`/`numpy` = 42; RDKit `BRICSBuild(seed=42)`.

### 4B.4 `04_filter_cns.py`
Per molecule compute (RDKit): MW, cLogP (Crippen `MolLogP`), TPSA, HBD, HBA,
RTB, aromatic-ring count, QED, fraction Csp3.
- **CNS-MPO (repo 5-check proxy)** — reproduce `26_fetch_tier2.py` exactly for
  continuity.
- **CNS-MPO (Wager 2016, 6-parameter)** — cLogP, cLogD (approximated as
  cLogP − 0.5·(fraction ionisable), documented), MW, TPSA, HBD, most-basic pKa
  (rule-based estimate; documented as an approximation).
- **Gates (from TODO §4):** MW < 450 AND cLogP 2–4 AND TPSA < 90 AND
  CNS-MPO(proxy) ≥ 4.0 AND PAINS-free AND anchor retained.
- Rank passers by CNS-MPO(proxy) then QED.

### 4B.5 `05_report.py`
- `results/phase7/glue_design/generated_library.csv` (every molecule, all
  props, per-gate pass flags).
- `results/phase7/glue_design/glue_candidates_top.csv` (passers, ranked).
- Top-24 2D grid PNG; top-N (≤50) SDF with 3D conformers.
- `results/phase7/glue_design/CONCLUSION.md` — methods paragraph + a Results
  sentence; explicitly frames this as a **hypothesis-generating scaffold
  library, not validated binders**, consistent with the manuscript's
  "methodological demonstration" framing of the screening section.

### 4B decision rule
Deliverable, not a test. Success = ≥ 10 anchor-preserving, CNS-druglike
candidates with documented provenance + a methods paragraph. If < 10 pass all
gates, relax the tightest single gate (likely cLogP 2–4 → 1–4), record which,
and report both counts.

### 4B optional (decide at approval)
Dock the top ≤ 25 candidates into the 8RQC ternary interface (Vina box on the
QFC centroid) as a positional-plausibility check — cheap, not MD. Include or skip?

---

## Deliverables

```
code/phase7_robustness/4A_negctrl/01_negative_control_analysis.py
code/phase7_robustness/4B_glue_gen/{01_extract_anchor.py,02_build_fragment_pool.py,
  03_generate_glues.py,04_filter_cns.py,05_report.py}
code/phase7_robustness/README.md          (append Section 4)
results/phase7/negctrl/{negative_control_scorecard.csv,CONCLUSION.md}
results/phase7/glue_design/{anchor.json,fragment_pool.csv,generated_raw.csv/.sdf,
  generated_library.csv,glue_candidates_top.csv,glue_grid.png,glue_top.sdf,CONCLUSION.md}
```

## Risks
- **QFC SMILES from coordinates** may be imperfect (bond orders) — mitigate with
  the PDB CCD lookup and a hand-checked fallback; QFC is only a reference, not
  load-bearing.
- **BRICSBuild explosion / long runtime** — cap enumeration, time-box, sample.
- **pKa approximation** in Wager MPO — clearly labelled; the repo 5-check proxy
  is the gate of record, Wager is secondary.
- **Over-claiming** — CONCLUSION.md and manuscript text must state these are
  unvalidated in-silico scaffolds; no binding claims.
