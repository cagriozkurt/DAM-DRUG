# Section 1 (WP1) — JASPAR 2026 cisTarget re-run: CONCLUSION

**Date:** 2026-09-11/12
**Verdict:** Under a properly-diagnosed JASPAR2026 cisTarget re-run, motif
enrichment does **not** support a direct, motif-driven regulatory link
between IKZF1 (or any of the other 10 comparator TFs — SPI1, RUNX1, IRF8,
PPARG, CEBPB, IKZF2, IKZF3, RELB, BHLHE40, BHLHE41) and its GRNBoost2
co-expression-inferred targets. This is a genuine limitation of the
motif-based validation approach, not a pipeline artifact — see the
diagnostic below. It is a distinct, separate question from Section 2's
regional-LMM and donor-robustness results, which stand independently of
this motif analysis and are unaffected by this finding.

## What was built and run

- `s1_02`: two cisTarget rankings databases from JASPAR2026 CORE vertebrate
  non-redundant motifs (1,097 motifs / 887 unique TFs) scored against
  aertslab's two canonical hg38 region windows — wide (10kbp-up/10kbp-down
  full transcript) and narrow (500bp-up/100bp-down full transcript) —
  scored together as aertslab's own production pipelines do (fixes an
  initial run that used the wide window alone).
- `s1_03`: `pyscenic ctx` re-pruned the SAME 5-seed GRNBoost2 consensus
  adjacency table used in the accepted paper against both windows, then
  `aucell` on the retained regulons.
- `s1_04`: benchmark/null-permutation suite submitted (job 6349095) for
  completeness against whatever regulons this DB retains.

## Result

- Single wide window (job 6346282): 9 regulons retained total, **0/11**
  target TFs (E2F1/E2F2/E2F7/E2F8/ETV6/FOSL2/NFIB/PBX3/ZNF148 passed instead
  — all strong, universal, easily-detected motifs unrelated to the
  hypothesis).
- Dual wide+narrow window (job 6348013): 11 regulons retained — down from
  **46** retained under Paper 1's original cisTarget v10/HOCOMOCO v11
  combination — and **still 0/11** target TFs.
  `results/phase7/grn_jaspar2026/TF_regulon_summary_jaspar2026.csv` confirms
  all 11 target TFs (including IKZF1) as `NOT_RETAINED`.
- `s1_04` benchmark/null suite (job 6349095, COMPLETED) is fully concordant:
  IKZF1 regulon pseudotime rho = NaN (0 target genes, never retained as a
  regulon so there is nothing to correlate); paralogue test and curated
  hypergeometric enrichment likewise degenerate (NaN/undefined) for the same
  reason. No inconsistency with the NES diagnosis below.

## Diagnosis (ruling out DB-construction artifacts before accepting the result)

Two DB-construction hypotheses were tested and ruled out as sufficient
explanations:

1. **Window-width dilution** — a single wide window disproportionately
   dilutes enrichment for TFs with large GRNBoost2 target sets (IKZF1 has
   5,967 candidate targets, well above the 455-target median across all
   1,784 TFs in the adjacency table). Fixed by adding the matched narrow
   proximal-promoter window aertslab pairs with it in production — result
   barely changed (9→11 total regulons, target TFs still 0/11).
2. **Motif-universe size** — `s1_01` deliberately used JASPAR2026 CORE
   **non-redundant** vertebrates (1,097 motifs) rather than the full CORE
   (2,633 PFMs) + UNVALIDATED collections the original TODO spec named,
   giving `ctx`'s per-module null-AUC estimate a smaller motif sample to
   draw from. Plausible contributing factor, not ruled out by direct test,
   but superseded by finding 3 below.
3. **Direct NES measurement** (`s1_05_diagnose_nes.py`, reproducing
   `pyscenic ctx`'s own recovery-curve/NES calculation via `ctxcore`'s
   public API for the single most favourable case — IKZF1's own top-50
   targets by GRNBoost2 importance, i.e. the tightest, highest-confidence
   module available): IKZF1's own motif (MA1508.2) scores **NES = −0.02**
   (wide window) and **NES = +0.57** (narrow window) — indistinguishable
   from the null (NES≈0), nowhere near `ctx`'s default `nes_threshold=3.0`.
   This is not a near-miss recoverable by adjusting thresholds or expanding
   the motif set: at the best-case target module, IKZF1's motif shows
   essentially zero recovery signal. Notably, *other* (non-IKZF1-annotated)
   motifs DO reach NES>3 for this same target-gene set (e.g. MA2126.1,
   MA2628.1, MA2546.1 for the wide window) — the targets share some
   regulatory signature, just not one matching IKZF1's JASPAR2026 PWM.

## Interpretation

GRNBoost2 infers TF-target edges from expression co-variation, which need
not reflect direct TF-DNA binding at a target's promoter/enhancer region.
The regional-LMM (§2.1) and donor-robustness (§2.2) results, which operate
directly on IKZF1 expression and AUCell regulon activity (not on motif
enrichment), are unaffected by this finding — they answer a different
question (is the association robust to regional/donor confounding) than
this one (does IKZF1's motif explain its own inferred target set). This
result should be reported as a distinct, honest limitation: cisTarget-style
motif enrichment, even with a correctly-constructed dual-window JASPAR2026
database, does not corroborate a direct motif-driven mechanism for IKZF1's
(or the 10 comparator TFs') co-expression network in this dataset.

## Files

- `TF_regulon_summary_jaspar2026.csv`, `regulon_auc_by_state_jaspar2026.csv`
  (11 retained regulons, all non-target TFs)
- `code/phase7_robustness/s1_grn/s1_05_diagnose_nes.py` (NES diagnostic)
- Large intermediates (rankings feathers, loom, full regulons CSV) live on
  TRUBA scratch, not committed to git (size).
