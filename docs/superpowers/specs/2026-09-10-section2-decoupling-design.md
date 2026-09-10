# Section 2 — Decoupling Evidence & Cross-Cohort Replication: Design Spec

**Date:** 2026-09-10
**Scope:** `TODO.md` Section 2 only (one sub-project of the DAM-DRUG robustness campaign).
**Target:** bioRxiv preprint update / stronger resubmission. No JAD Oct-26 deadline dependency.
**Compute:** local (this Mac) for 2A/2B/2D; 2C best-effort open-access, execution deferred.

---

## Objective

Resolve two reviewer-anticipated attacks on the *IKZF1* late-microglial-state finding:

1. **Regional confounding** — the *IKZF1* trajectory association could reflect
   which brain regions contribute which substates, not a disease-state transition.
2. **Survivorship / outlier-donor bias** — the LateAD-DAM regulon peak
   (AUCell = 0.153) rests on a sparse state (median 4 cells/donor, 66/84 donors)
   and could be driven by one or two donors.

Plus: make the non-independence of the internal SEA-AD evidence lines explicit
in the manuscript (2D), and add whatever external single-cell replication is
reachable without credentialed data access (2C).

The manuscript already names 2A + 2C as planned future work
(`Manuscript.md` Discussion, lines 745–762) and already carries honest
non-independence caveats (lines 419–433, 716–726). This sub-project executes
those promises and upgrades the caveats to demonstrated robustness.

---

## Data inventory (verified present)

| Asset | Path | Use |
|---|---|---|
| Multi-region microglia RNA | `data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad` (3.0 GB) | 2A pseudobulk |
| Trajectory object | `results/phase1/trajectory/microglia_trajectory.h5ad` | 2A/2B labels: `state`, `dpt_pseudotime`, `Brain Region`, `Donor ID`, `Supertype`, `Age at Death`, `Sex`, `Braak`, `Overall AD neuropathological Change` |
| AUCell scores | `results/phase2/GRN/scenic_auc_aggregated.loom` (100K cells × 46 regulons; `ca["RegulonsAUC"]`, `ca["CellID"]`, `ca["state"]`) | 2B |
| Pseudotime table | `results/phase1/trajectory/pseudotime.csv` (`cell_id`, `dpt_pseudotime`) | 2B jackknife |
| Supertype→state map | hard-coded in `code/phase1_QC/05_donor_region_stats.py:57-64` | all |
| Existing pseudobulk logic | `code/slurm/04_pseudobulk_deseq2.slurm` (design `~group + sex_bin + age_z`, min 10 cells/donor-state) | 2A reference |
| Existing ρ correlation | `code/phase2_GRN/03_regulon_pseudotime_correlation.py` (ρ = +0.30892) | 2B jackknife reference |

`Brain Region` categories observed: AnG, DFC, FI, HIP, ITG, LEC, MEC, MTG (+ up to 10 total).
No AUCell columns in any h5ad `obs` — AUCell lives only in the looms.
`data/raw/Mathys2019/` is **empty**. Mathys 2023 / Sun 2023 microglia are Synapse-credentialed.

---

## Environment (prerequisite — Task 0)

New conda env `damdrug`:

```
conda create -n damdrug -c conda-forge -c bioconda python=3.11 \
  scanpy anndata h5py loompy pandas numpy scipy scikit-learn statsmodels matplotlib
```

R 4.6.1 already present with working `lme4`; add `lmerTest`, `broom.mixed` via
`install.packages`. Record exact resolved versions into
`envs/damdrug.yml` and `envs/damdrug_R_sessionInfo.txt`.

---

## Task 2A — Regional confounding via linear mixed-effects models

**Directory:** `code/phase7_robustness/2A_lmm/`

### 2A.1 Build multi-region pseudobulk — `01_build_pseudobulk.py`
- Load multi-region microglia h5ad (`backed="r"`; iterate to control memory).
- Map `Supertype` → 6 substates via `SUPERTYPE_TO_STATE`. Drop `DAM-IRM`
  (manuscript treats it as an analytically ambiguous hybrid — same handling
  as the primary DGE).
- Group by (`Donor ID`, `state`, `Brain Region`). Sum raw counts. Emit **two
  variants**: (a) *unconstrained* — all donor×state×region bins retained
  regardless of size; (b) *filtered* — bins with ≥ 10 cells only (matches
  manuscript pseudobulk minimum). Because LateAD-DAM has a cohort-wide median
  of ~4 cells/donor, very few LateAD-DAM donor×region bins survive (b); record
  the exact surviving-unit count per state in `state_x_region_group_counts.csv`.
- Primary regional LMM target contrast is **DAM vs Homeostatic** across all
  regions (well-powered). LateAD-DAM is reported from **both** the unconstrained
  and filtered pseudobulk, with the surviving donor–region unit count stated
  explicitly alongside every LateAD-DAM coefficient.
- Normalise: CPM then `log1p`. Also emit raw summed counts (for a NB-GLMM
  sensitivity check).
- Covariate table per group: `Age at Death` (→ z-score `age_z`), `Sex`,
  `Braak`, `Overall AD neuropathological Change`, `n_cells`.
- Output: `results/phase7/lmm/pseudobulk_multiregion.csv` (long: group_id,
  donor, state, region, gene set restricted to `IKZF1, IRF8, SPI1, BHLHE41,
  RUNX1, CEBPB, PPARG` + housekeeping controls `ACTB, GAPDH`),
  `pseudobulk_group_meta.csv`, and a coverage table
  `state_x_region_group_counts.csv`.

### 2A.2 Fit LMMs — `02_fit_lmm.R`
For each target gene, fit with `lmerTest::lmer` on `log1p(CPM)`:

- **Full model:** `expr ~ state + age_z + Sex + (1 | donor) + (1 | region)`
- **Reduced model:** `expr ~ state + age_z + Sex + (1 | donor)`
- **No-state null:** `expr ~ age_z + Sex + (1 | donor) + (1 | region)`

Report, to `results/phase7/lmm/<gene>_lmm.csv`:
- Fixed-effect coefficients for each `state` level (ref = Homeostatic),
  SE, 95 % CI (Wald + profile), Satterthwaite p-value.
- Random-effect variance components; ICC for donor and for region.
- `anova(reduced, full)` LRT (does adding `(1|region)` matter?) and
  `anova(no_state_null, full)` LRT (does state matter after region?).
- Marginal & conditional R².

Companion: `03_lmm_forest_plot.py` → `results/phase7/lmm/ikzf1_state_forest.pdf/.png`
(state coefficients ± CI, full vs reduced side by side, for IKZF1 with the
other 6 TFs as a small-multiple).

### 2A.3 Sensitivity — `04_lmm_sensitivity.R`
- NB-GLMM on raw counts with `offset(log(total_umi))` (`glmmTMB` if available,
  else `lme4::glmer.nb`) for IKZF1.
- Region as fixed effect instead of random (few-levels robustness).
- Leave-one-region-out refits (drop MTG; drop each region once) — does the
  DAM / LateAD-DAM IKZF1 coefficient stay positive and CI-exclusive of 0?

### Decision rule (record in `results/phase7/lmm/CONCLUSION.md`)
- **PASS** — DAM and/or LateAD-DAM IKZF1 coefficient remains > 0 with 95 % CI
  excluding 0 in the full model *and* in ≥ 8/ n leave-one-region-out refits
  → manuscript claim stands; add "robust to regional composition (LMM,
  `(1|region)`, ΔAIC / LRT reported)".
- **PARTIAL** — significant in full model but sensitive to dropping one region
  → keep claim, add explicit region-specific caveat naming that region.
- **FAIL** — coefficient loses significance once `(1|region)` added
  → downgrade to "region-associated (MTG-weighted)"; propagate to Abstract,
  Results, Discussion, Fig 4.

---

## Task 2B — LateAD-DAM donor downsampling & outlier robustness

**Directory:** `code/phase7_robustness/2B_downsample/`

### 2B.1 Assemble per-cell table — `01_assemble_auc_table.py`
- Read `scenic_auc_aggregated.loom`: `IKZF1(+)` AUCell + `CellID` + `state`.
- Join donor + `dpt_pseudotime` from trajectory h5ad `obs` (by CellID) and
  `pseudotime.csv`. **Barcode reconciliation:** pySCENIC loom export commonly
  drops/alters trailing suffixes (`-1`, batch prefixes). If the raw join
  overlap is < 90 %, retry after normalising both sides
  (`barcode.split('-')[0]`, strip known batch prefixes); abort loudly if still
  < 90 %.
- Output: `results/phase7/downsample/ikzf1_auc_cells.parquet`
  (cell_id, donor, state, ikzf1_auc, dpt).

### 2B.2 Leave-one-donor-out — `02_loo_donor.py`
- Baseline: per-state mean `IKZF1(+)` AUCell; confirm LateAD-DAM ≈ 0.153 peak.
- For each of the 66 LateAD-DAM donors: drop that donor's cells, recompute
  per-state means and the LateAD-DAM − Homeostatic gap.
- Output: `loo_donor_lateaddam.csv` (donor_dropped, lateaddam_mean, gap,
  rank_of_lateaddam_among_states); flag any donor whose removal drops
  LateAD-DAM out of the top state or halves the gap.

### 2B.3 Iterative downsampling — `03_downsample_cells.py`
- For `k in {1,2,3,4}` cells/donor: 1,000 draws, sample ≤ k cells per
  LateAD-DAM donor, recompute LateAD-DAM mean AUCell and gap vs Homeostatic.
- Output: `downsample_distribution.csv` + violin plot; report the fraction of
  draws where the gap stays > 0 and where LateAD-DAM remains the top state.

### 2B.4 Pseudotime ρ jackknife — `04_rho_jackknife.py`
- Recompute the IKZF1(+) AUCell vs `dpt_pseudotime` Spearman ρ
  (baseline +0.309), dropping one donor at a time (all donors, not only
  LateAD-DAM), and via 1,000 donor-bootstrap resamples.
- Output: `rho_jackknife.csv` (ρ distribution, min/max, 95 % percentile CI),
  `rho_jackknife.pdf`.

### Decision rule (`results/phase7/downsample/CONCLUSION.md`)
- **PASS** — LateAD-DAM stays the top-AUCell state in ≥ 63/66 LOO refits,
  downsampled gap > 0 in ≥ 95 % of draws, and ρ 95 % CI excludes 0
  → add "not driven by outlier donors (LOO + downsampling + bootstrap CI)".
- **FAIL** — a single donor's removal collapses the peak, or ρ CI crosses 0
  → reframe LateAD-DAM regulon peak as exploratory / underpowered; soften
  Results item and Fig 4.

---

## Task 2C — External single-cell replication (best-effort, open-access only)

**Directory:** `code/phase7_robustness/2C_external/`
**Execution deferred** — build scripts now, run when/if data is reachable.

### 2C.1 Acquisition probe — `01_fetch_open_access.sh`
Attempt, in order, without credentials:
- **Primary open-access fallbacks** (both have pre-processed expression
  matrices + author-defined reactive/homeostatic cluster labels on GEO FTP):
  - **Olah et al. 2020 — GSE146639** (live-sorted human microglia scRNA-seq).
  - **Gerrits et al. 2021 — GSE148822** (human microglia, AD parietal + occipital cortex).
- Mathys 2019 (GSE138852-adjacent) / Sun 2023 GEO or Synapse mirror if any
  matrix is reachable without login.
- UCSC Cell Browser processed `exprMatrix.tsv.gz` + `meta.tsv` for any AD
  microglia dataset with reactive/homeostatic annotation.
- Log exactly what resolved to `results/phase7/external/ACQUISITION_LOG.md`.
If nothing usable resolves: stop, document, leave 2C as an explicit open
limitation (unchanged from current manuscript wording).

### 2C.2 Signature scoring — `02_score_external.py` (runs only if 2C.1 succeeds)
- Build the `IKZF1(+)` target gene set from
  `results/phase2/GRN/regulons_aggregated.csv`.
- On each external microglia matrix: `scanpy` normalise, `sc.tl.score_genes`
  (and AUCell-style rank scoring as cross-check).
- Test: signature score in the cohort's annotated reactive / DAM / lipid
  cluster vs homeostatic — Mann–Whitney U, Cliff's delta, plot.
- If ≥ 500 microglia and a usable embedding: independent PAGA + DPT, correlate
  signature vs pseudotime.
- Output: `results/phase7/external/<cohort>_signature_test.csv`, figure.

### Decision rule
- **PASS** — signature enriched (p < 0.05, delta > 0.2) in the reactive
  cluster of ≥ 1 external cohort → add as a genuine second external line
  (alongside GSE95587); update 2D block.
- **NULL / no data** — state plainly: "external single-nucleus replication
  was attempted with open-access cohorts X, Y; [enrichment not observed /
  no annotated cohort was reachable]." No claim inflation.

---

## Task 2D — Epistemic classification update (manuscript + figures)

**Directory:** `code/phase7_robustness/2D_epistemic/`
**Do last**, after 2A–2C conclusions are known.

- Regenerate the Fig 4 / Fig 3D evidence-summary table so internal SEA-AD RNA
  modalities (pseudobulk DGE, RcisTarget regulon, DPT correlation, CellOracle,
  ATAC) render as **one non-independent block**, with GSE95587 and any 2C
  cohort as separate external rows. Script: `01_evidence_table.py` →
  `results/phase7/epistemic/evidence_summary.csv` + updated figure panel.
- Text patches (tracked in `results/phase7/epistemic/manuscript_patches.md`,
  applied to `Manuscript.md` as a separate reviewed commit):
  - Abstract: state "internal analyses are non-independent" explicitly.
  - Results "Multi-evidence prioritisation" section: the block relabelling.
  - Discussion: fold in 2A + 2B outcomes.
  - Key Points, Fig 2 / Fig 4 legends: consistent language.
- No numeric claim changes here unless 2A/2B decision rules trigger a
  downgrade — those are separate commits with their own rationale.

---

## Deliverables

```
envs/damdrug.yml, envs/damdrug_R_sessionInfo.txt
code/phase7_robustness/2A_lmm/{01_build_pseudobulk.py,02_fit_lmm.R,03_lmm_forest_plot.py,04_lmm_sensitivity.R}
code/phase7_robustness/2B_downsample/{01_assemble_auc_table.py,02_loo_donor.py,03_downsample_cells.py,04_rho_jackknife.py}
code/phase7_robustness/2C_external/{01_fetch_open_access.sh,02_score_external.py}
code/phase7_robustness/2D_epistemic/01_evidence_table.py
code/phase7_robustness/README.md   (run order, env, decision-rule summary)
results/phase7/lmm/**, results/phase7/downsample/**, results/phase7/external/**, results/phase7/epistemic/**
results/phase7/lmm/CONCLUSION.md, results/phase7/downsample/CONCLUSION.md
```

Each analysis script: `DAM_DRUG_DIR` env-var path convention (matches the
repo's post-audit standard), fixed `random_state=42`, logging, `--help`.

## Out of scope for this sub-project
- pySCENIC / JASPAR 2026 re-run (Section 1).
- SLIT2→ROBO2 multi-region CellChat (Section 3).
- Any molecular dynamics or drug-discovery work (Section 4).
- Credentialed Synapse downloads.
- Zenodo deposition (Section 5, done once campaign complete).

## Risks
- **Memory:** 3 GB h5ad on a Mac — mitigate with `backed="r"` + chunked
  pseudobulk accumulation; never densify the full matrix.
- **CellID join mismatch** between loom and trajectory h5ad (barcode suffix
  formatting) — 2B.1 must verify overlap ≥ 90 % and abort loudly otherwise
  (same guard as `03_regulon_pseudotime_correlation.py:75`).
- **2C dry well** — likely; treated as an acceptable documented null, not a
  blocker.
- **LMM singular fits** for sparse state×region cells (LateAD-DAM) —
  expected; report `isSingular`, fall back to region-as-fixed-effect, keep
  LateAD-DAM interpretation cautious.
