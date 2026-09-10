# Section 2 — Manuscript patches

Apply to `results/manuscript/04_jad/revision/Manuscript.md` as a **separate,
reviewed commit** (or into the bioRxiv branch). Line numbers are from the
2026-08-29 revision. Each patch is scoped to what Section 2 (2A/2B/2C) proved.

---

## Patch 1 — Limitations: region-stratified modelling (was "not done")

**Location:** Discussion / *Limitations*, ~line 718-721.

**Remove:**
> The absence of region-stratified pseudobulk modelling means some
> transcriptional differences may reflect regional composition rather than
> disease-state transitions.

**Replace with:**
> Region-stratified linear mixed-effects modelling of pseudobulk expression
> across all ten SEA-AD regions (`expr ~ substate + age + sex + (1|donor) +
> (1|region)`) confirmed that elevated *IKZF1* in DAM (+0.57 log1p-CPM, 95% CI
> 0.53-0.61) and LateAD-DAM (+0.55, 95% CI 0.44-0.65) remained significant
> after modelling regional variance (region intraclass correlation approximately 0.04). The
> substate effect survived a likelihood-ratio test against a no-substate model
> (P < 1e-90), every leave-one-region-out refit (11/11), region as a fixed
> effect, and a negative-binomial GLMM on raw counts (Supplementary Fig. S-new).

---

## Patch 2 — Limitations: LateAD-DAM sparsity / outlier robustness (new sentence)

**Location:** Discussion / *Limitations*, after the region sentence above.

**Add:**
> The LateAD-DAM regulon peak, though based on a sparse state (median 4
> cells/donor), was robust to donor sampling: across leave-one-donor-out
> refits (57 LateAD-DAM donors) LateAD-DAM remained the highest-AUCell
> substate in every case, capping every LateAD-DAM donor at a single cell
> preserved this ranking in 1,000/1,000 resamples, and the IKZF1(+)
> AUCell-pseudotime correlation was stable under a donor jackknife
> (rho 0.30-0.33) and bootstrap (95% CI 0.27-0.35).

---

## Patch 3 — External replication now has a second (state-level) line

**Location:** Results, *Multi-evidence prioritisation of IKZF1*, ~line 427-433;
and Discussion *Limitations* ~line 716-718.

**In Results**, after the GSE95587 sentence ("...adjusted P value = .004."),
**add:**
> In an independent open-access single-nucleus cohort (Grubman et al. 2019,
> entorhinal cortex; GSE138852), the IKZF1(+) regulon signature was enriched
> in the reactive microglial subcluster relative to the homeostatic subcluster
> (Cliff's delta +0.72; Mann-Whitney P = 4e-17; random-gene-set permutation
> P = 0.001) but was not elevated by crude AD-versus-control status, mirroring
> the state-linked (rather than diagnosis-linked) association seen in SEA-AD.

**In Limitations**, change:
> analysis 6 (GSE95587) is the sole external replication, in a single cohort
> and brain region.

**to:**
> external replication rests on one bulk cohort (GSE95587) plus one
> open-access single-nucleus cohort (Grubman 2019), the latter supporting the
> state association but not a donor-diagnosis difference; credentialed
> single-nucleus cohorts (Mathys 2019/2023, Sun 2023) were not accessed.

---

## Patch 4 — Multi-evidence section: reclassify the six analyses as a block

**Location:** Results, *Multi-evidence prioritisation of IKZF1*, ~line 415-434.

The paragraph already labels analyses 1-5 as non-independent — keep that.
**Append** one sentence at the end of the paragraph:
> Because analyses 1-5 share the same SEA-AD RNA data, they constitute a single
> non-independent evidence block; the independent lines are the SEA-AD ATAC
> assay (item 5, orthogonal), the GSE95587 bulk cohort, the Grubman 2019
> single-nucleus cohort, and the Ballasch 2023 western blot. Robustness checks
> (regional LMM; donor leave-one-out and downsampling) do not add independent
> evidence but establish that the internal block is not an artefact of regional
> composition or donor sparsity.

Update the **Figure 3D / Figure 4 evidence table** to the grouping in
`results/phase7/epistemic/evidence_summary.csv` (five blocks: internal RNA,
internal ATAC, internal robustness, external, literature).

---

## Patch 5 — Abstract: make non-independence explicit

**Location:** Abstract (~line 25-40).

If not already present, add a clause to the sentence listing the convergent
evidence, e.g. after "...multi-modality convergence...":
> (the internal SEA-AD analyses are non-independent; external support is one
> bulk and one single-nucleus cohort)

---

## Patch 6 — Future priorities: item (d) is now done

**Location:** Discussion, *Immediate experimental priorities*, ~line 760-762.

**Remove** item (d):
> and (d) region-stratified pseudobulk modelling across the ten SEA-AD brain
> regions to determine whether the *IKZF1* trajectory association is
> cortex-wide or region-specific.

(now reported in Results / Limitations per Patch 1). Renumber or drop.

---

## NOT changed by Section 2 (flagged for other sections)

- Line 719-722 states the *IKZF1* finding is "conditional on the cisTarget
  v10/HOCOMOCO v11 combination" and that BHLHE40/41 "could not be evaluated
  because HOCOMOCO v11 lacks their atypical E-box motifs". The rest of the
  manuscript (Methods, audit) describes **cisTarget v10** rankings, not
  HOCOMOCO v11. This inconsistency is Section 1 territory (pySCENIC / JASPAR
  2026 re-run) — do not fix here, but resolve it there.
- SLIT2->ROBO2 single-region caveat (line 723-725) — Section 3.
- Virtual screening / MD reframing (line 726-728) — Section 4.
