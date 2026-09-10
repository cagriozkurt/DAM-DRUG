# Section 2B — LateAD-DAM outlier / survivorship robustness: CONCLUSION

**Date:** 2026-09-10
**Verdict: PASS** — the LateAD-DAM IKZF1(+) regulon peak and the
AUCell-vs-pseudotime correlation are not driven by outlier donors or by
uneven per-donor cell contribution.

## Inputs
- `scenic_auc_aggregated.loom` — 100,000-cell pySCENIC subsample, `IKZF1(+)` AUCell.
- `microglia_trajectory.h5ad` obs — state, donor, `dpt_pseudotime` (100% join).
- LateAD-DAM in this subsample: 637 cells, 57 donors, median 4 cells/donor.
- Baseline LateAD-DAM mean AUCell = **0.15265** (manuscript: 0.153); gap vs
  Homeostatic = **0.03260**; LateAD-DAM is the top-AUCell state (rank 1/6).
- Baseline Spearman rho(IKZF1(+) AUCell, dpt) = **+0.30892** (manuscript: +0.309).

## Results

### 2B.2 Leave-one-donor-out (57 LateAD-DAM donors)
- LateAD-DAM lost rank 1: **0 / 57** refits.
- Gap vs Homeostatic more than halved: **0 / 57**.
- Gap range across all refits: **[0.03124, 0.03316]** (baseline 0.03260).
- Most influential single donor (H21.33.002, 95 LateAD-DAM cells): gap still 0.0312, rank still 1.

### 2B.3 Iterative per-donor downsampling (1,000 draws per cap)
| cap k | mean gap | 95% CI | frac gap > 0 | frac LateAD-DAM top state |
|---|---|---|---|---|
| 1 | 0.03049 | [0.02705, 0.03389] | 1.000 | 1.000 |
| 2 | 0.03048 | [0.02768, 0.03298] | 1.000 | 1.000 |
| 3 | 0.03020 | [0.02792, 0.03240] | 1.000 | 1.000 |
| 4 | 0.03032 | [0.02825, 0.03239] | 1.000 | 1.000 |

Even at 1 cell per LateAD-DAM donor, LateAD-DAM remains the highest-AUCell
state in every one of 1,000 draws.

### 2B.4 Pseudotime rho jackknife / bootstrap
- Leave-one-donor-out (all 84 donors): rho range **[+0.302, +0.327]**, sign stable in all 84.
- Donor-level bootstrap (1,000 resamples): mean **+0.309**, 95% CI **[+0.266, +0.347]**, 100% positive.

## Manuscript action
Add to Results / Methods (Section 2B robustness supplement):

> The LateAD-DAM IKZF1(+) regulon peak was robust to donor sparsity: across
> leave-one-donor-out refits (57 LateAD-DAM donors) LateAD-DAM remained the
> highest-AUCell substate in every case (gap vs homeostatic 0.0312-0.0332),
> and capping every LateAD-DAM donor at a single cell preserved the ranking
> in 1,000/1,000 resamples. The IKZF1(+) AUCell-pseudotime correlation was
> stable under a donor jackknife (rho 0.302-0.327) and a donor bootstrap
> (95% CI 0.27-0.35).

No numeric claims require downgrading.

## Files
`ikzf1_auc_cells.parquet`, `loo_donor_lateaddam.{csv,png}`,
`downsample_distribution.{csv,png}`, `rho_jackknife.csv`, `rho_bootstrap.csv`,
`rho_jackknife.png`.
