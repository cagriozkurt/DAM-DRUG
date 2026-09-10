# Phase 7 — Robustness & de novo glue campaign (Paper #2)

Analyses that harden the DAM-DRUG findings against "artefact" / confounding
critiques and add a de novo cereblon-glue design arm. **Not part of the
accepted JAD paper** (phases 1–6, frozen at tag `v1.0-jad-accepted`). This work
lives on branch `paper2-robustness-glue` and feeds a second manuscript.
`TODO.md` at repo root is the roadmap; verdicts are in
`results/phase7/*/CONCLUSION.md`.

Specs: `docs/superpowers/specs/2026-09-10-section{2,4}-*-design.md`

## Environments

| used by | env | key packages |
|---|---|---|
| Section 2 (genomics) | `damdrug` (`envs/damdrug.yml`) | scanpy 1.11.5, anndata 0.12, loompy 3.0.8, statsmodels; R 4.6 + lme4/lmerTest/broom.mixed/glmmTMB |
| Section 4 (chem) | `lipogate` | rdkit 2026.03.1, vina 1.2.7, meeko 0.7.1 |

All scripts use `DAM_DRUG_DIR` (fallback = cwd), seed 42.

```
export DAM_DRUG_DIR=/path/to/DAM-DRUG
conda env create -f envs/damdrug.yml
```

## Section 2 — Decoupling evidence & cross-cohort replication

| step | script | runtime | output |
|---|---|---|---|
| 2A.1 | `2A_lmm/01_build_pseudobulk.py` | ~1 min | `results/phase7/lmm/pseudobulk_*.csv` |
| 2A.2 | `2A_lmm/02_fit_lmm.R` | ~30 s | `results/phase7/lmm/<gene>_lmm.csv`, `lmm_model_comparison.csv` |
| 2A.3 | `2A_lmm/04_lmm_sensitivity.R` | ~10 s | `lmm_leave_one_region_out.csv`, `lmm_region_fixed_effect.csv`, `lmm_nbglmm_ikzf1.csv` |
| 2A.fig | `2A_lmm/03_lmm_forest_plot.py` | ~5 s | `ikzf1_state_forest.{pdf,png}` |
| 2B.1 | `2B_downsample/01_assemble_auc_table.py` | ~20 s | `results/phase7/downsample/ikzf1_auc_cells.parquet` |
| 2B.2 | `2B_downsample/02_loo_donor.py` | ~5 s | `loo_donor_lateaddam.{csv,png}` |
| 2B.3 | `2B_downsample/03_downsample_cells.py` | ~15 s | `downsample_distribution.{csv,png}` |
| 2B.4 | `2B_downsample/04_rho_jackknife.py` | ~10 s | `rho_jackknife.csv`, `rho_bootstrap.csv`, `rho_jackknife.png` |
| 2C.1 | `2C_external/01_fetch_open_access.sh` | ~1 min | `results/phase7/external/raw/GSE138852_*` |
| 2C.2 | `2C_external/02_score_external.py` | ~1 min | `grubman_signature_test.csv`, `grubman_signature.png` |
| 2D.1 | `2D_epistemic/01_evidence_table.py` | ~5 s | `results/phase7/epistemic/evidence_summary.{csv,png}` |

### Section 2 result summary

- **2A regional confounding — PASS.** IKZF1 DAM (+0.57) and LateAD-DAM (+0.55)
  survive `(1|region)`, region-as-fixed, 11/11 leave-one-region-out, NB-GLMM.
- **2B outlier / survivorship — PASS.** LateAD-DAM peak and pseudotime rho
  stable to any single donor and to 1-cell/donor downsampling.
- **2C external replication — PASS (state-level).** IKZF1(+) signature enriched
  in reactive vs homeostatic microglia in Grubman 2019 (delta +0.72, perm
  p .001); not an AD-vs-control effect (consistent with a state-linked signal).
- **2D** — evidence table regrouped; internal SEA-AD analyses = one
  non-independent block.

### Section 2 deviations from the spec

- **2A pseudobulk source:** built from `results/phase1/trajectory/microglia_trajectory.h5ad`
  (raw counts, 36,601 genes, 236,002 microglia, all 10 regions, carries the
  `state` label) instead of the 3 GB `SEA-AD_Microglia_multi-regional_*.h5ad`
  (which lacks `state` and adds monocyte/lymphocyte contaminants). Same cells,
  lighter, no barcode join needed.
- **2C cohort:** Grubman 2019 (GSE138852) used as the open-access fallback;
  Olah 2020 full matrix (inside `GSE146639_RAW.tar`) deferred; Mathys/Sun not
  accessed (credentialed).

## Section 4 — Drug-pipeline reframe + CRBN molecular-glue generation

| step | script | output |
|---|---|---|
| 4A.1 | `4A_negctrl/01_negative_control_analysis.py` | `results/phase7/negctrl/negative_control_scorecard.csv`, `.png`, `CONCLUSION.md` |
| 4B.1 | `4B_glue_gen/01_extract_anchor.py` | `results/phase7/glue_design/anchor.json`, `anchor.png` |
| 4B.2 | `4B_glue_gen/02_build_fragment_pool.py` | `fragment_pool.csv` (561 BRICS fragments) |
| 4B.3 | `4B_glue_gen/03_generate_glues.py` | `generated_raw.csv/.sdf` (2,750 anchor-preserving products; ~2 min) |
| 4B.4 | `4B_glue_gen/04_filter_cns.py` | `generated_library.csv` (gate flags) |
| 4B.5 | `4B_glue_gen/05_report.py` | `glue_candidates_top.csv`, `glue_top.sdf`, `glue_grid.png` |
| 4B.6 | `4B_glue_gen/06_dock_glues.py` | `glue_docking.csv`, `docked/*.pdbqt` (Vina into 8RQC box; ~10 min) |

### Section 4 result summary

- **4A — reframe.** Tafamidis (CHEMBL2103837, MM-GBSA rank 1/5 for IRF8) and
  diflunisal (CHEMBL898, rank 1/10 for PPARG) both leave the pocket in 100 ns
  explicit-solvent MD (core-RMSD 25.17 Å / 83.42 Å) while lower-ranked IRF8
  compounds stay bound → implicit-solvent + shallow/mis-assigned-pocket
  artefact. Recast as MM-GBSA sensitivity benchmarks / negative controls.
  See `results/phase7/negctrl/CONCLUSION.md`.
- **4B — generation.** Lenalidomide isoindolinone–glutarimide anchor grown at
  the 4-amino position with 561 BRICS fragments × 5 linkers → 2,750
  anchor-preserving products. **The strict TODO gates (TPSA < 90, CNS-MPO ≥ 4)
  are unsatisfiable**: the warhead alone has TPSA 92.5 Å² and cLogP ≈ 0. 184
  candidates pass a pre-registered fallback (TPSA < 120, CNS-MPO(proxy) ≥ 3.5,
  MW < 450, cLogP 2–4, PAINS-free); top 25 dock into the 8RQC ternary interface
  at Vina −5.9…−7.0 kcal/mol (bare anchor −5.2). Unvalidated scaffolds, not
  binders. See `results/phase7/glue_design/CONCLUSION.md`.

## Not yet done (other TODO.md sections)

Section 1 (pySCENIC / JASPAR 2026 — TRUBA), Section 3 (SLIT2→ROBO2 multi-region
CellChat — TRUBA), Section 4 heavy tail (explicit-solvent MD / T-REMD, hERG
QSAR, counter-docking), Section 5 (Zenodo, container digests).
Paper #2 manuscript: `Manuscript_Paper2.md` (to be drafted).
