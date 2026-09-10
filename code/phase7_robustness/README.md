# Phase 7 — Robustness campaign

Analyses that harden the DAM-DRUG findings against "artefact" / confounding
critiques for a bioRxiv update / stronger resubmission. Not part of the
original JAD pipeline (phases 1-6).

Design spec: `docs/superpowers/specs/2026-09-10-section2-decoupling-design.md`

## Environment

```
conda env create -f envs/damdrug.yml      # scanpy 1.11.5, anndata 0.12, loompy 3.0.8, statsmodels
# R 4.6: install.packages(c("lme4","lmerTest","broom.mixed","glmmTMB"))
export DAM_DRUG_DIR=/path/to/DAM-DRUG
```

All scripts use `DAM_DRUG_DIR` (fallback = cwd), seed 42.

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

Verdicts are in `results/phase7/*/CONCLUSION.md`; manuscript edits in
`results/phase7/epistemic/manuscript_patches.md`.

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

## Deviations from the design spec

- **2A pseudobulk source:** built from `results/phase1/trajectory/microglia_trajectory.h5ad`
  (raw counts, 36,601 genes, 236,002 microglia, all 10 regions, carries the
  `state` label) instead of the 3 GB `SEA-AD_Microglia_multi-regional_*.h5ad`
  (which lacks `state` and adds monocyte/lymphocyte contaminants). Same cells,
  lighter, no barcode join needed.
- **2C cohort:** Grubman 2019 (GSE138852) used as the open-access fallback;
  Olah 2020 full matrix (inside `GSE146639_RAW.tar`) deferred; Mathys/Sun not
  accessed (credentialed).

## Not yet done (other TODO.md sections)

Section 1 (pySCENIC / JASPAR 2026), Section 3 (SLIT2->ROBO2 multi-region),
Section 4 (drug reframe + molecular-glue generation), Section 5 (Zenodo).
