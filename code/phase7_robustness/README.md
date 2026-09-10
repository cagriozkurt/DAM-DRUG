# Phase 7 — Robustness campaign

Analyses that harden the DAM-DRUG findings against "artefact" / confounding
critiques for a bioRxiv update / stronger resubmission. Not part of the
original JAD pipeline (phases 1–6). `TODO.md` at repo root is the roadmap.

Specs: `docs/superpowers/specs/2026-09-10-section*-design.md`

- **Section 2** (regional LMM, donor downsampling, external replication,
  epistemic reclassification) — branch `robustness/section2-decoupling`.
- **Section 4** (this branch, `robustness/section4-drug`) — drug-pipeline
  reframe + CRBN molecular-glue generation.

All scripts use `DAM_DRUG_DIR` (fallback = cwd), seed 42.

## Section 4 — env

RDKit / Vina / meeko work runs in the **`lipogate`** conda env
(rdkit 2026.03.1, vina 1.2.7, meeko 0.7.1). No new MD.

```
export DAM_DRUG_DIR=/path/to/DAM-DRUG
P=/opt/homebrew/Caskroom/miniconda/base/envs/lipogate/bin/python
```

## Section 4 — run order

| step | script | output |
|---|---|---|
| 4A.1 | `4A_negctrl/01_negative_control_analysis.py` | `results/phase7/negctrl/negative_control_scorecard.csv`, `.png`, `CONCLUSION.md` |
| 4B.1 | `4B_glue_gen/01_extract_anchor.py` | `results/phase7/glue_design/anchor.json`, `anchor.png` |
| 4B.2 | `4B_glue_gen/02_build_fragment_pool.py` | `fragment_pool.csv` (561 BRICS fragments) |
| 4B.3 | `4B_glue_gen/03_generate_glues.py` | `generated_raw.csv/.sdf` (2,750 anchor-preserving products; ~2 min) |
| 4B.4 | `4B_glue_gen/04_filter_cns.py` | `generated_library.csv` (gate flags) |
| 4B.5 | `4B_glue_gen/05_report.py` | `glue_candidates_top.csv`, `glue_top.sdf`, `glue_grid.png` |
| 4B.6 | `4B_glue_gen/06_dock_glues.py` | `glue_docking.csv`, `docked/*.pdbqt` (Vina into 8RQC box; ~10 min) |

## Section 4 result summary

- **4A — reframe.** Tafamidis (CHEMBL2103837, MM-GBSA rank 1/5 for IRF8) and
  diflunisal (CHEMBL898, rank 1/10 for PPARG) both leave the pocket in 100 ns
  explicit-solvent MD (core-RMSD 25.17 Å / 83.42 Å) while lower-ranked IRF8
  compounds stay bound → implicit-solvent + shallow/mis-assigned-pocket
  artefact. Recast as MM-GBSA sensitivity benchmarks / negative controls;
  retract as candidates. See `results/phase7/negctrl/CONCLUSION.md`.
- **4B — generation.** Lenalidomide isoindolinone–glutarimide anchor grown at
  the 4-amino position with 561 BRICS fragments × 5 linkers → 2,750
  anchor-preserving products. **The strict TODO gates (TPSA < 90, CNS-MPO ≥ 4)
  are unsatisfiable**: the warhead alone has TPSA 92.5 Å² and cLogP ≈ 0, so
  no lenalidomide-based glue can reach them. 184 candidates pass a
  pre-registered fallback (TPSA < 120, CNS-MPO(proxy) ≥ 3.5, MW < 450,
  cLogP 2–4, PAINS-free). Framed as unvalidated scaffolds, not binders.
  See `results/phase7/glue_design/CONCLUSION.md`.
