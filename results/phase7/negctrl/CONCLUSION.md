# Section 4A — Tafamidis & diflunisal as MM-GBSA negative controls: CONCLUSION

**Date:** 2026-09-10
**Action: reframe.** Retract tafamidis and diflunisal as repurposing
candidates. Recast them as MM-GBSA sensitivity benchmarks that expose the
implicit-solvent + shallow/mis-assigned-pocket failure mode of the screening
pipeline.

## The discordance (from existing phase-3/4 outputs — no new computation)

| compound | on-target | Vina | gnina CNN | MM-GBSA ΔG | MM-GBSA rank in target set | 100 ns MD core-RMSD | outcome |
|---|---|---|---|---|---|---|---|
| **Tafamidis** (CHEMBL2103837) | IRF8 AF2-DBD | −6.71 | 0.74 | **−9.5 ± 0.05** | **1 / 5** | 25.17 ± 3.56 Å | ligand egress |
| **Diflunisal** (CHEMBL898) | PPARG 1FM9-LBD | −8.77 | 0.96 | −2.8 ± 0.21 | **1 / 10** | 83.42 ± 2.89 Å | complete egress |
| CHEMBL42 (context) | IRF8 AF2-DBD | −6.70 | 0.74 | −6.9 ± 0.14 | 3 / 5 | — | retained |
| CHEMBL490 (context) | IRF8 AF2-DBD | −6.68 | 0.75 | −7.1 ± 0.32 | 2 / 5 | — | retained |

MD core-RMSD values are from the manuscript / Supp Fig S8 (archived TRUBA
trajectories; `gen_seed = -1`, not code-reproducible).

## Root cause

Both molecules were **rank 1** by MM-GBSA ΔG within their target's compound set,
yet both leave the binding site within 100 ns of explicit-solvent MD, while the
2nd/3rd-ranked IRF8 compounds stay bound.

- **Tafamidis / IRF8:** favourable implicit-solvent ΔG on a small (500 Å³),
  weakly hydrophobic (11.8) pocket on an AlphaFold2 model flagged `LOW`
  confidence at the pocket. Implicit solvent removes the desolvation penalty
  that, in explicit water, ejects the ligand from a shallow surface groove.
- **Diflunisal / PPARG:** the entire PPARG MM-GBSA set has ΔG ≈ 0 (−2.8 to
  +0.4); "rank 1" is meaningless. The composite score nonetheless promoted it
  (Vina −8.77, CNN 0.96), and it then shows the largest egress of any compound
  simulated.

Neither result is evidence about IRF8 or PPARG biology; both are artefacts of
scoring shallow / non-native TF surface sites with implicit-solvent MM-GBSA and
a docking composite.

## Manuscript action (patch, separate reviewed commit)

1. **Abstract** (lines 33, 39, 779-780): remove "Tafamidis (→IRF8) and
   diflunisal (→PPARG)" as nominated repurposing hypotheses. Replace with a
   sentence that the virtual-screening arm is a methodological demonstration
   and that its two apparent top hits fail explicit-solvent MD.
2. **Results** (lines 479-486, 659-660): relabel tafamidis/diflunisal as
   "MM-GBSA sensitivity benchmarks / negative controls", not "→ IRF8" /
   "→ PPARG" candidates. State the rank-1-but-egress discordance explicitly.
3. **Tables 3 / 4 / 5** and **Figure 5A legend** (lines 498, 510, 520-524,
   1186-1187): mark the tafamidis and diflunisal rows as
   "benchmark (fails MD) — not advanced"; drop "top-ranked".
4. **Discussion / Limitations** (lines 726-728): keep the existing sentence
   that implicit-solvent MM-GBSA gives relative rankings not affinities;
   add that the two nominal top hits were used to quantify this limitation.
5. Do **not** touch the CRBN-track hits here (celecoxib, teriflunomide, …) —
   they remain exploratory and are addressed by 4B / Section 4 CRBN pivot.

## Files
`negative_control_scorecard.csv`, `mmgbsa_vs_md_retention.png`.
