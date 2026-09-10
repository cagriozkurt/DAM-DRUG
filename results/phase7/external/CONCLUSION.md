# Section 2C — External single-cell replication (open-access): CONCLUSION

**Date:** 2026-09-10
**Verdict: PASS (state-level), NULL (diagnosis-level)** — in an independent
open-access cohort the discovery IKZF1(+) regulon signature is strongly
enriched in reactive vs homeostatic microglia, but is *not* elevated by crude
AD-vs-control status. This is the pattern the SEA-AD analysis predicts: the
signal tracks the microglial activation trajectory, not donor diagnosis.

## Data acquired (no credentials)
- **Grubman et al. 2019, GSE138852** — entorhinal cortex snRNA-seq.
  449 microglia nuclei, subclusters m1-m5, per-cell AD/control labels.
  `GSE138852_counts.csv.gz`, `GSE138852_covariates.csv.gz`.
- Not attempted (Synapse-credentialed): Mathys 2019 (syn18485175), Mathys 2023, Sun 2023.
- Olah 2020 (GSE146639): full matrix only inside `GSE146639_RAW.tar` (per-sample);
  the GEO-level `readinCBC.csv.gz` is a bare barcode list. Deferred.

## Signature
IKZF1(+) regulon from `regulons_aggregated.csv` — 1,271 target genes
(union over 4 motif rows); 1,072 present in the Grubman matrix.
Scored with `scanpy.tl.score_genes` (control gene bins, seed 42).

## Results
| contrast | median (test / ref) | MWU p | Cliff's delta | random-geneset perm p |
|---|---|---|---|---|
| reactive subcluster **m2 vs m1** (homeostatic) | 0.163 / 0.099 | **4.1e-17** | **+0.72** | **0.001** |
| crude **AD vs control** microglia | 0.140 / 0.153 | 0.999 (n.s.) | −0.18 | 1.0 |
| raw *IKZF1* transcript, AD vs ct | 0.60 / 0.63 | 0.59 | — | — |

Per-subcluster signature mean: m1 0.101 < m4 0.122 < m5 0.144 < m2 0.164 ≈ m3 0.168.
The activated subclusters (m2, m3) carry the highest IKZF1(+) activity;
homeostatic m1 the lowest.

## Interpretation
- The IKZF1 regulon signature **replicates as a microglial-state marker** in an
  independent cohort, brain region (entorhinal cortex vs SEA-AD MTG-anchored
  trajectory), and platform — reactive microglia score ~0.65 SD higher than
  homeostatic, beyond any equal-size random gene set (permutation p = 0.001).
- It does **not** separate AD from control donors in Grubman. This is
  concordant with the SEA-AD result, where the association is with diffusion
  pseudotime / substate (ρ = +0.31) rather than with donor diagnosis, and it
  argues against the SEA-AD signal being a disease-status batch artefact.

## Manuscript action
Add a sentence to Results (external replication) and update the 2D evidence table:

> In an independent open-access cohort (Grubman et al. 2019, entorhinal cortex
> snRNA-seq; GSE138852), the IKZF1(+) regulon signature was enriched in the
> reactive microglial subcluster relative to the homeostatic subcluster
> (Cliff's delta +0.72; Mann-Whitney p = 4e-17; random-gene-set permutation
> p = 0.001), but was not elevated by crude AD-versus-control status —
> mirroring the state-linked (rather than diagnosis-linked) association seen in
> SEA-AD. Credentialed single-nucleus cohorts (Mathys 2019/2023, Sun 2023)
> were not accessed; this remains a single open-access replication.

## Files
`ikzf1_regulon_targets.json`, `ACQUISITION_LOG.md`, `grubman_signature_test.csv`,
`grubman_subcluster_means.csv`, `grubman_signature.png`.
