# Section 4 tail — In Silico Off-Target Selectivity (CONCLUSION)

- Candidates counter-docked: 25 (top-25 by rank, results/phase7/glue_design/glue_candidates_top.csv)
- SI > 5.0 (all 3 off-targets: DRD2/HTR2A/HERG) after ddG->SI correction: 0/25
- hERG QSAR clean (predicted liability prob < 0.5): 23/25
- Pass BOTH gates: 0/25

## Methods note: SI formula correction
Phase 1's selectivity_table.csv used SI = vina_ontarget / vina_offtarget (raw ratio of two Vina kcal/mol scores of similar magnitude), which is mathematically bounded near 1 and cannot reach a SI > 5.0 gate for any molecule -- all 6 of its own compounds scored 0.58-0.88. This script uses SI = exp(ddG / RT) (ddG = deltaG_offtarget - deltaG_ontarget, RT = 0.593 kcal/mol @ 298 K), the standard free-energy-to-affinity-ratio proxy, so the gate is achievable and physically interpretable.

## Top 10 by min_si
| glue_id   |   vina_ontarget_8rqc |   vina_drd2 |   vina_htr2a |   vina_herg |     min_si |   herg_liability_prob | passes_selectivity_gate   |
|:----------|---------------------:|------------:|-------------:|------------:|-----------:|----------------------:|:--------------------------|
| GLUE1189  |               -6.855 |      -9.015 |       -9.759 |      -8.785 | 0.00746796 |              0.410458 | False                     |
| GLUE0872  |               -6.453 |      -9.068 |       -9.375 |      -8.066 | 0.00724468 |              0.422549 | False                     |
| GLUE0141  |               -6.172 |      -9.626 |       -9.602 |      -8.821 | 0.00295392 |              0.420224 | False                     |
| GLUE0870  |               -6.082 |      -9.619 |       -9.896 |      -8.2   | 0.00160971 |              0.445603 | False                     |
| GLUE0335  |               -6.563 |     -10.411 |       -9.884 |      -8.428 | 0.00152001 |              0.515709 | False                     |
| GLUE2036  |               -6.67  |     -10.548 |       -9.427 |      -8.675 | 0.00144502 |              0.442755 | False                     |
| GLUE1185  |               -6.104 |      -8.904 |      -10.038 |      -8.461 | 0.00131481 |              0.443372 | False                     |
| GLUE1089  |               -5.922 |      -9.232 |       -9.88  |      -8.774 | 0.00126266 |              0.445833 | False                     |
| GLUE0059  |               -6.684 |     -10.652 |       -9.797 |      -8.676 | 0.00124154 |              0.625097 | False                     |
| GLUE0246  |               -7.006 |     -10.993 |      -10.179 |      -8.715 | 0.00120239 |              0.395386 | False                     |
