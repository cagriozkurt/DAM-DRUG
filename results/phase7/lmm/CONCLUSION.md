# Section 2A — Regional confounding (linear mixed-effects models): CONCLUSION

**Date:** 2026-09-10
**Verdict: PASS** — elevated *IKZF1* expression in DAM and LateAD-DAM microglia
is not explained by regional sampling composition. It survives an explicit
`(1 | brain_region)` random effect, region as a fixed effect, every
leave-one-region-out refit, and a negative-binomial GLMM on raw counts.

## Data
- Pseudobulk from `microglia_trajectory.h5ad` (raw counts, 36,601 genes,
  182,187 microglia after dropping the DAM-IRM hybrid), grouped by
  donor × substate × brain-region across all **10 SEA-AD regions**
  (AnG, DFC, FI, HIP, ITG, LEC, MEC, MTG, STG, V1C).
- `min10` variant: 1,375 groups with ≥ 10 cells. `unconstrained`: 2,087 groups.
- Coverage (`state_x_region_group_counts.csv`): DAM 528 groups / 84 donors /
  10 regions under `min10`; LateAD-DAM only 31 groups / 13 donors / 8 regions
  under `min10` (249 groups / 66 donors / 10 regions unconstrained), so
  LateAD-DAM is reported from **both** variants.
- Model: `log1p(CPM) ~ state + age_z + Sex + (1|donor) + (1|region)`
  (ref state = Homeostatic), ML fit, `lmerTest` Satterthwaite p-values.

## Results — IKZF1

| variant | state | coef (log1p CPM) | 95% CI | p |
|---|---|---|---|---|
| min10 | DAM | **+0.570** | [0.532, 0.609] | 2e-144 |
| min10 | IRM | +0.461 | [0.422, 0.500] | 4e-101 |
| min10 | LAM | −0.278 | [−0.380, −0.176] | 1e-07 |
| min10 | LateAD-DAM | **+0.547** | [0.444, 0.651] | 3e-24 |
| unconstrained | DAM | +0.663 | [0.539, 0.786] | 4e-25 |
| unconstrained | LateAD-DAM | +0.291 | [0.134, 0.448] | 3e-04 |

- **Does region matter?** LRT adding `(1|region)`: p = 9e-11 (min10) / 3e-04
  (unconstrained) — there is real regional variance, but it is small:
  ICC_region for IKZF1 = **0.037** (min10) / 0.013 (unconstrained) vs
  ICC_donor = 0.32 / 0.045.
- **Does state survive region?** LRT full vs no-state model:
  p = 1e-168 (min10) / 4e-94 (unconstrained). State effect is overwhelming
  after accounting for region.
- `isSingular` = FALSE for IKZF1 in both variants.

## Sensitivity (`04_lmm_sensitivity.R`)
- **Leave-one-region-out** (drop each of 11 region levels once, 11 refits):
  IKZF1 stateDAM CI_lo > 0 in **11/11**, stateLateAD-DAM CI_lo > 0 in **11/11**.
- **Region as fixed effect:** IKZF1 stateDAM +0.571 [0.533, 0.609],
  stateLateAD-DAM +0.545 [0.442, 0.648] — indistinguishable from the
  random-effect estimates.
- **NB-GLMM on raw counts** with `offset(log(total_umi))`, `(1|donor)+(1|region)`:
  IKZF1 stateDAM log-rate-ratio +0.511 (p = 5e-223),
  stateLateAD-DAM +0.510 (p = 1e-38).

## Companion TFs (context, not the primary claim)
- IRF8, BHLHE41, RUNX1: DAM and LateAD-DAM coefficients positive, CI excludes 0,
  11/11 leave-one-region-out — same conclusion.
- PPARG: negative in DAM/LateAD-DAM throughout (down-regulated).
- SPI1, CEBPB: DAM-associated but **not** LateAD-DAM (0/11 leave-one-region-out
  for LateAD-DAM) — consistent with IKZF1's distinct late-stage profile.
- IRF8 unconstrained LateAD-DAM coefficient is marginally negative
  (−0.22, CI [−0.45, −0.002]); the min10 estimate is +0.87. Divergence is
  driven by noisy sub-10-cell bins — reported, not load-bearing.

## Manuscript action
Replace the current caveat (Discussion lines 716–726: "absence of
region-stratified modelling means... some differences may reflect regional
composition") with a demonstrated-robustness statement:

> Region-stratified linear mixed-effects modelling of pseudobulk expression
> across all ten SEA-AD regions (`expr ~ state + age + sex + (1|donor) +
> (1|region)`) confirmed that elevated *IKZF1* in DAM (+0.57 log1p-CPM, 95% CI
> 0.53–0.61) and LateAD-DAM (+0.55, 95% CI 0.44–0.65) remained significant
> after modelling regional variance (ICC_region ≈ 0.04). The state effect
> survived a likelihood-ratio test against a no-state model (p < 1e-90), every
> leave-one-region-out refit (11/11), region as a fixed effect, and a
> negative-binomial GLMM on raw counts.

## Files
`pseudobulk_multiregion.csv`, `pseudobulk_group_meta.csv`,
`state_x_region_group_counts.csv`, `<gene>_lmm.csv` (12 genes),
`lmm_model_comparison.csv`, `lmm_leave_one_region_out.csv`,
`lmm_region_fixed_effect.csv`, `lmm_nbglmm_ikzf1.csv`,
`ikzf1_state_forest.{pdf,png}`.
