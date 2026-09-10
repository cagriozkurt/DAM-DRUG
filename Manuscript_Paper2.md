<!--
DRAFT — Paper #2. Not for submission. Built on branch paper2-robustness-glue.
The accepted JAD paper (Paper #1) is frozen at git tag v1.0-jad-accepted and is
NOT modified by this work. This document is the follow-up study.
Numbers are drawn from results/phase7/*/CONCLUSION.md; regenerate with
code/phase7_robustness/.
-->

# De Novo Cereblon Molecular-Glue Design and Multi-Region Mixed-Effects Validation of IKZF1 in Late-Stage Alzheimer's Disease Microglia

*(provisional title)*

**Çağrı Özkurt**¹ *(author list and affiliations to be finalised)*

¹ Hacettepe University, Ankara, Türkiye. Correspondence: cagriozkurt@hacettepe.edu.tr

---

## Abstract

**Background.** We previously applied an integrative single-nucleus framework to
the SEA-AD atlas and nominated the transcription factor *IKZF1* (Ikaros) as a
candidate regulator associated with a terminal disease-associated microglial
state (LateAD-DAM) in Alzheimer's disease (AD) [Paper #1, *J Alzheimers Dis*
2026]. That discovery rested on non-independent internal analyses of one cohort
and pointed to an orthosteric target pocket that structural analysis judged
undruggable.

**Objective.** To test whether the *IKZF1*–late-state association is robust to
regional sampling, state sparsity, and an independent cohort; and to define a
central-nervous-system (CNS)-compliant chemical strategy for engaging IKZF1
through targeted protein degradation rather than orthosteric inhibition.

**Methods.** (i) Pseudobulk expression was modelled across all ten SEA-AD brain
regions with linear mixed-effects models (`expr ~ substate + age + sex +
(1|donor) + (1|region)`), leave-one-region-out refits, and a negative-binomial
GLMM. (ii) The LateAD-DAM *IKZF1(+)* regulon peak and its pseudotime correlation
were stress-tested by leave-one-donor-out, iterative per-donor downsampling
(1,000 draws), and donor bootstrap. (iii) The discovery *IKZF1(+)* signature was
projected onto an independent open-access snRNA-seq cohort (Grubman et al. 2019,
entorhinal cortex; GSE138852). (iv) The two apparent virtual-screening hits from
Paper #1 (tafamidis, diflunisal) were re-analysed as MM-GBSA negative controls.
(v) A scaffold-based library of candidate cereblon (CRBN) molecular glues was
enumerated on the lenalidomide isoindolinone–glutarimide warhead by BRICS
fragmentation of CNS-approved chemical space, filtered for CNS druglikeness, and
docked into the CRBN–IKZF1(ZF2) ternary interface (PDB 8RQC).

**Results.** Elevated *IKZF1* in DAM (+0.57 log1p-CPM, 95% CI 0.53–0.61) and
LateAD-DAM (+0.55, 95% CI 0.44–0.65) survived explicit modelling of regional
variance (region intraclass correlation ≈ 0.04), a likelihood-ratio test
against a no-substate model (P < 10⁻⁹⁰), all eleven leave-one-region-out
refits, region as a fixed effect, and a negative-binomial GLMM on raw counts.
The LateAD-DAM regulon peak remained the highest-activity substate in 57/57
leave-one-donor-out refits and in 1,000/1,000 draws that capped every
LateAD-DAM donor at a single cell; the *IKZF1(+)* AUCell–pseudotime correlation
was stable under a donor jackknife (ρ 0.30–0.33) and bootstrap (95% CI
0.27–0.35). In the independent cohort, the *IKZF1(+)* signature was enriched in
the reactive microglial subcluster relative to the homeostatic subcluster
(Cliff's δ = +0.72; Mann–Whitney P = 4 × 10⁻¹⁷; random-gene-set permutation
P = 10⁻³) but was not elevated by crude AD-versus-control status, mirroring the
state-linked (not diagnosis-linked) association seen in SEA-AD. Re-analysed as
negative controls, tafamidis and diflunisal each ranked first by MM-GBSA ΔG
within their target's compound set yet left the binding site within 100 ns of
explicit-solvent molecular dynamics (ligand core-RMSD 25.17 Å and 83.42 Å),
whereas lower-ranked compounds were retained — an implicit-solvent,
shallow-pocket artefact. Of 2,750 anchor-preserving generated glues, none met a
strict CNS filter (topological polar surface area [TPSA] < 90 Å², CNS-MPO ≥ 4):
the bifunctional glutarimide warhead alone (TPSA 92.5 Å², cLogP ≈ 0) precludes
it. 184 candidates passed a pre-registered relaxed filter; the top 25 docked
into the 8RQC ternary interface with Vina scores of −5.9 to −7.0 kcal/mol, all
more favourable than the bare anchor (−5.2 kcal/mol).

**Conclusions.** The *IKZF1*-associated late-disease microglial state withstands
regional-confounding, sparsity, and cross-cohort challenges, strengthening it as
a target hypothesis. The undruggable orthosteric pocket can be bypassed by a
cereblon molecular-glue strategy, but the glutarimide warhead places such
molecules at the boundary of CNS druglikeness — a design constraint that any
brain-penetrant degrader programme must confront. All candidates here are
unvalidated computational scaffolds requiring biophysical and cellular testing.

**Keywords.** Alzheimer's disease; disease-associated microglia; IKZF1/Ikaros;
linear mixed-effects models; targeted protein degradation; cereblon; molecular
glue; CNS multiparameter optimisation.

---

## 1. Introduction

Disease-associated microglia (DAM) accumulate around amyloid pathology and
undergo a stereotyped transcriptional progression from a homeostatic state
through intermediate reactive states to terminal, lipid- and
stress-associated phenotypes. In an integrative analysis of 236,002 microglial
nuclei from 84 SEA-AD donors, we nominated *IKZF1* (Ikaros) as a transcription
factor whose regulon activity peaks in a terminal state we termed LateAD-DAM
(AUCell 0.153) and is the only one of 46 retained regulons positively
correlated with diffusion pseudotime (Spearman ρ = +0.309) [Paper #1]. That
work also found the *IKZF1* orthosteric zinc-finger surface undruggable
(fpocket drug_score 0.001) and its two apparent small-molecule hits, tafamidis
and diflunisal, unstable in molecular dynamics.

Two limitations block translation of that discovery. First, the single-cell
evidence was **non-independent**: pseudobulk differential expression, regulon
inference, pseudotime correlation, and in-silico perturbation all draw on the
same SEA-AD RNA from the same 84 donors, and no region-stratified model was
fitted despite well-established regional heterogeneity of microglial states.
The late-state peak also rests on a sparse population (median 4 cells per donor,
66/84 donors), raising a survivorship concern. Second, the nominated target is
a transcription factor with no classical binding pocket; conventional
orthosteric medicinal chemistry is not a viable path.

Here we address both bottlenecks as an independent follow-up study. For the
genomic bottleneck we apply mixed-effects modelling across all ten SEA-AD
regions, donor-resampling stress tests, and projection onto an independent
snRNA-seq cohort, and we re-cast the evidence architecture to separate
non-independent internal analyses from genuine external replication. For the
chemical bottleneck we first use the two failed virtual-screening hits as an
explicit case study of implicit-solvent scoring pitfalls, then pivot to
targeted protein degradation: we enumerate a scaffold-based library of
candidate cereblon (CRBN) molecular glues on the lenalidomide
isoindolinone–glutarimide warhead, grow it toward the IKZF1 ZF2 β-hairpin
degron using the CRBN–IKZF1 ternary structure (PDB 8RQC), and filter for CNS
druglikeness. The accepted discovery paper is unchanged; this study builds on
top of it.

---

## 2. Results

### 2.1 The IKZF1-associated late state is not a regional-composition artefact

We built donor × substate × brain-region pseudobulk profiles from 182,187
microglia (the DAM-IRM hybrid excluded) spanning all ten SEA-AD regions
(AnG, DFC, FI, HIP, ITG, LEC, MEC, MTG, STG, V1C) and fitted, for *IKZF1* and
six companion transcription factors, the linear mixed-effects model
`log1p(CPM) ~ substate + age_z + sex + (1|donor) + (1|region)`
(reference substate = Homeostatic; **Figure 1A**).

Relative to homeostatic microglia, *IKZF1* was elevated in DAM
(β = +0.57 log1p-CPM, 95% CI 0.53–0.61, P = 2 × 10⁻¹⁴⁴) and in LateAD-DAM
(β = +0.55, 95% CI 0.44–0.65, P = 3 × 10⁻²⁴), and reduced in the branch state
LAM (β = −0.28). Adding the `(1|region)` term improved fit (likelihood-ratio
test P = 9 × 10⁻¹¹) but regional variance was small — the region intraclass
correlation for *IKZF1* was ≈ 0.04, against a donor intraclass correlation of
0.32. The substate effect was overwhelming after accounting for region: a
likelihood-ratio test of the full model against a model without substate gave
P < 10⁻⁹⁰ in both the ≥ 10-cell and unconstrained pseudobulk variants.

Three sensitivity analyses agreed (**Figure 1B–C**). (i) In eleven
leave-one-region-out refits, the DAM and LateAD-DAM *IKZF1* coefficients
remained positive with 95% confidence intervals excluding zero in **11/11**
cases. (ii) With region as a fixed rather than random effect, the estimates
were unchanged (DAM +0.57 [0.53, 0.61]; LateAD-DAM +0.55 [0.44, 0.65]).
(iii) A negative-binomial GLMM on raw counts with a library-size offset gave
*IKZF1* log-rate-ratios of +0.51 for both DAM (P = 5 × 10⁻²²³) and LateAD-DAM
(P = 1 × 10⁻³⁸). The companion factors *IRF8*, *BHLHE41*, and *RUNX1* behaved
identically; *PPARG* was consistently down-regulated; *SPI1* and *CEBPB* were
DAM-associated but **not** LateAD-DAM-associated (0/11 leave-one-region-out),
underscoring the distinct late-stage profile of *IKZF1*.

### 2.2 The LateAD-DAM regulon peak is not driven by donor sparsity or outliers

In the 100,000-cell pySCENIC subsample, LateAD-DAM comprised 637 cells from 57
donors (median 4 cells/donor) with a mean *IKZF1(+)* AUCell of 0.1527 — the
highest of the six substates — and a gap of 0.0326 above homeostatic microglia
(**Figure 2A**).

- **Leave-one-donor-out** (57 LateAD-DAM donors): LateAD-DAM remained the
  highest-AUCell substate in **57/57** refits; the homeostatic gap never more
  than halved and stayed within 0.0312–0.0332.
- **Iterative downsampling** (1,000 draws per cap): capping every LateAD-DAM
  donor at k = 1, 2, 3, or 4 cells preserved a positive gap and top-substate
  ranking in **1,000/1,000** draws at every k; at k = 1 the mean gap was 0.030
  (95% CI 0.027–0.034) (**Figure 2B**).
- **Pseudotime correlation**: the *IKZF1(+)* AUCell–diffusion-pseudotime
  Spearman ρ (baseline +0.309) ranged 0.302–0.327 across an 84-donor jackknife
  (sign stable in all) and had a donor-bootstrap 95% CI of 0.27–0.35, entirely
  positive (**Figure 2C**).

### 2.3 Independent single-nucleus replication of the state association

We projected the discovery *IKZF1(+)* target-gene signature (1,271 genes;
1,072 present) onto 449 microglial nuclei from an independent open-access
cohort (Grubman et al. 2019, entorhinal cortex snRNA-seq; GSE138852), which
carries author-defined microglial subclusters and per-nucleus AD/control
labels (**Figure 3**).

The signature was strongly enriched in the reactive subcluster (m2) relative to
the homeostatic subcluster (m1): median score 0.163 vs 0.099, Mann–Whitney
P = 4.1 × 10⁻¹⁷, Cliff's δ = +0.72, and a random-gene-set permutation P = 10⁻³
(1,000 size-matched gene sets). It was **not** elevated by crude
AD-versus-control status (Cliff's δ = −0.18; not significant), and raw *IKZF1*
transcript did not differ by diagnosis. This is the pattern predicted by the
SEA-AD analysis, where the association is with microglial state (pseudotime,
substate) rather than with donor diagnosis, and it argues against the SEA-AD
signal being a disease-status batch artefact. Credentialed single-nucleus
cohorts (Mathys 2019/2023, Sun 2023) were not accessed; this is a single
open-access replication.

### 2.4 A revised evidence architecture for IKZF1

We re-classified the evidence for *IKZF1* into five blocks by epistemic status
(**Table 1**; **Figure 3D of Paper #1** superseded):

| block | lines | independence |
|---|---|---|
| Internal SEA-AD RNA | pseudobulk DGE; RcisTarget regulon; diffusion pseudotime; CellOracle KO | one non-independent block (same 84 donors, same assay) |
| Internal SEA-AD chromatin | ATAC motif enrichment (IKZF1 3.90×) | orthogonal assay, same donors |
| Internal robustness (this study) | 10-region LMM; donor LOO + downsampling | confirmatory, not independent evidence |
| External replication | GSE95587 bulk (log₂FC +0.64, P_adj = 0.004, n = 117); Grubman 2019 snRNA (δ +0.72, state-level) | independent cohorts |
| Literature triangulation | Ballasch 2023 western blot (IKZF1 protein ↑, AD hippocampus) | independent |

The robustness analyses in §2.1–2.2 do not add an independent line of evidence;
they establish that the internal block is not an artefact of regional
composition or donor sparsity.

### 2.5 Implicit-solvent MM-GBSA on shallow TF pockets: a negative-control case study

In Paper #1, structure-based virtual screening surfaced tafamidis
(CHEMBL2103837, nominally → IRF8) and diflunisal (CHEMBL898, nominally → PPARG).
Both ranked **first by MM-GBSA ΔG within their target's compound set**
(tafamidis ΔG = −9.5 kcal/mol, rank 1/5; diflunisal ΔG = −2.8 kcal/mol,
rank 1/10 — in a set where every ΔG ≈ 0), yet both left the binding site within
100 ns of explicit-solvent molecular dynamics (ligand core-RMSD 25.17 ± 3.56 Å
and 83.42 ± 2.89 Å over the final 20 ns), while the second- and third-ranked
IRF8 compounds (ΔG −6.9 and −7.1 kcal/mol) were retained (**Figure 4**).

The IRF8 pocket is small (≈ 500 Å³), weakly hydrophobic (fpocket
hydrophobicity 11.8) and lies on an AlphaFold2 model flagged low-confidence at
the pocket; the PPARG hit sits outside the canonical ligand-binding sub-site.
Implicit-solvent MM-GBSA removes the desolvation penalty that, in explicit
water, ejects a ligand from such a shallow or solvent-exposed surface groove,
and the docking composite then promotes it. We therefore present tafamidis and
diflunisal not as repurposing candidates but as **MM-GBSA sensitivity
benchmarks / negative controls** that quantify this failure mode; neither
result is evidence about IRF8 or PPARG biology.

### 2.6 A CNS-constrained cereblon molecular-glue library for the IKZF1 ZF2 degron

Because the *IKZF1* orthosteric zinc-finger pocket is undruggable
(drug_score 0.001), we redirected chemical design to the cereblon (CRBN)
ternary interface, where thalidomide-type agents recruit C2H2 zinc-finger
degrons — including the IKZF1 ZF2 β-hairpin — to the CRL4^CRBN E3 ligase
(**Figure 5A**; PDB 8RQC, CRBN chains A/D, IKZF1 ZF2 residues 144–170).

We fixed the lenalidomide **isoindolinone–glutarimide** warhead
(`O=C1CCC(N2Cc3cccc(N)c3C2=O)C(=O)N1`; MW 259.3) and set the growth exit vector
at its 4-amino position (the CELMoD linker-attachment point; the 8RQC reference
ligand grows from the analogous 4-position). We BRICS-decomposed 1,962
CNS-approved and approved-drug SMILES into **561** mono-attachment fragments
(MW 40–250, ≤ 2 rings, PAINS-free, no reactive groups) and joined each to the
warhead nitrogen through five linker chemistries (amide, reverse-amide,
sulfonamide, acetamide spacer, urea, carbamate), verifying that every product
retained the intact glutarimide + isoindolinone. This yielded **2,750** unique
anchor-preserving molecules (**Figure 5B**).

### 2.7 The glutarimide warhead sets a CNS-druglikeness ceiling

**No molecule met the strict CNS filter** specified a priori (TPSA < 90 Å²
*and* CNS-MPO ≥ 4.0), and this is structural rather than a sampling failure:
the isoindolinone–glutarimide warhead alone has TPSA 92.5 Å² and cLogP ≈ 0, so
TPSA < 90 is unreachable for any lenalidomide-based glue and the CNS-MPO
polar-surface term is capped. This places cereblon-glue chemical space at the
**boundary of CNS druglikeness by construction** — a constraint any
brain-penetrant degrader programme must plan around, e.g. by minimising
appended polarity or by warhead engineering.

Under a pre-registered relaxed filter (TPSA < 120 Å², CNS-MPO(5-check
proxy) ≥ 3.5, MW < 450, cLogP 2–4, PAINS-free), **184** candidates survived
(MW 371–449, cLogP 2.0–2.9, QED 0.49–0.80; **Table 2**). As a
positional-plausibility check — not an affinity prediction — the top 25 were
docked into the 8RQC ternary interface with the exact grid box from the
Paper #1 screen. All 25 docked inside the box at Vina −5.9 to −7.0 kcal/mol,
each 0.7–1.8 kcal/mol better than the bare lenalidomide anchor (−5.2), i.e. the
grown fragment adds interface contacts as intended, with no anomalous scores
(**Figure 5C–D**). The best-scoring scaffolds were simple N-aroyl derivatives
(e.g. 4-methyl-, 2-methyl-, 2-chloro-benzamide).

---

## 3. Discussion

This study hardens and re-frames the *IKZF1* microglial hypothesis from our
accepted discovery paper. The late-disease association is not an artefact of
which brain regions contribute which microglial states: it survives an explicit
`(1|region)` random effect, region as a fixed effect, every leave-one-region-out
refit, and a negative-binomial model on raw counts, with regional variance
contributing only ≈ 4% of the total. It is not an artefact of the sparse
LateAD-DAM population: the regulon peak and its pseudotime correlation withstand
removal of any single donor and downsampling to one cell per donor. And the
*state* association — though not a crude AD-versus-control difference —
replicates in an independent snRNA-seq cohort and brain region. The distinction
matters: a signal that tracks microglial activation state rather than donor
diagnosis is harder to explain as a cohort or batch confound and is the
biologically expected behaviour of a state-transition regulator.

We are deliberate about what is *not* independent. Four of the original six
"convergent" analyses share one RNA dataset; the robustness work here is
confirmatory, not corroborating. The genuinely external support remains two
cohorts (one bulk, one single-nucleus, the latter at state level only) plus one
published western blot. Credentialed single-nucleus cohorts and a
replication-scale trajectory analysis are the obvious next steps.

On chemistry, the tafamidis/diflunisal case study is a cautionary,
generalisable result: within-target MM-GBSA rank-1 status meant nothing here,
because implicit solvent over-stabilised ligands on shallow, solvent-exposed
transcription-factor surfaces that explicit-solvent dynamics immediately
rejected. Any TF-directed virtual screen that stops at implicit-solvent
rescoring risks the same false positives.

The pivot to cereblon molecular glues is mechanistically sound — the IKZF1 ZF2
degron is a validated CRBN neosubstrate motif — but our generation exercise
exposes a hard constraint: the bifunctional glutarimide warhead alone already
exceeds common CNS polar-surface limits, so brain-penetrant IKZF1 degraders
occupy a narrow and unforgiving property window. The 184 relaxed-filter
scaffolds and their ternary-interface docking are a starting point for
synthesis and biophysical triage, not a candidate set. Essential next
experiments include cereblon TR-FRET or competitive-binding assays for warhead
engagement, a ternary-complex assay (e.g. AlphaLISA) for IKZF1 recruitment,
and *IKZF1* degradation plus LateAD-DAM phenotype assays in human iPSC-derived
microglia.

**Limitations.** One external single-nucleus cohort, open-access only; no
replication-scale PAGA/pseudotime; the CNS-MPO of record is a five-check proxy,
not the full Wager pKa-dependent score; docking is rigid-receptor and is a
positional check, not an affinity estimate; no molecular dynamics, no
free-energy perturbation, and no experimental data are presented for the
generated glues.

---

## 4. Methods (condensed)

All analysis code is in `code/phase7_robustness/` (branch
`paper2-robustness-glue`); the accepted Paper #1 code and manuscript are frozen
at git tag `v1.0-jad-accepted` and are unchanged. Random seed 42 throughout.
Environments: `envs/damdrug.yml` (Python 3.11, scanpy 1.11.5, anndata 0.12,
loompy 3.0.8, statsmodels; R 4.6 with lme4 / lmerTest / glmmTMB) for the
genomic analyses; a cheminformatics environment with RDKit 2026.03.1,
AutoDock Vina 1.2.7, and Meeko 0.7.1 for the chemistry.

**Multi-region pseudobulk and LMMs (§2.1).** Raw counts from the trajectory
object (`results/phase1/trajectory/microglia_trajectory.h5ad`; 236,002 nuclei,
36,601 genes, all 10 regions) were mapped to six substates, the DAM-IRM hybrid
removed, and summed by donor × substate × region. Two variants were retained:
all bins, and bins with ≥ 10 cells. CPM-normalised, log1p-transformed
expression was modelled with `lmerTest::lmer` (ML), with Satterthwaite p-values
and Wald confidence intervals; a no-substate model and a no-region model were
compared by likelihood-ratio test. Sensitivity: leave-one-region-out refits,
region as a fixed effect, and `lme4::glmer.nb` on raw counts with
`offset(log(total_UMI))`.

**Donor robustness (§2.2).** `IKZF1(+)` AUCell scores
(`results/phase2/GRN/scenic_auc_aggregated.loom`, 100,000-cell subsample) were
joined to donor / substate / diffusion-pseudotime labels. Leave-one-donor-out
and per-donor downsampling (k ∈ {1,2,3,4}, 1,000 draws) recomputed the
per-substate mean AUCell and the LateAD-DAM − homeostatic gap. The
AUCell–pseudotime Spearman ρ was recomputed under an 84-donor jackknife and a
1,000-replicate donor bootstrap.

**External replication (§2.3).** GSE138852 counts and covariates were retrieved
from GEO. Microglia (n = 449) were normalised (`scanpy`, target 10⁴, log1p) and
scored with `sc.tl.score_genes` on the *IKZF1(+)* target set from
`results/phase2/GRN/regulons_aggregated.csv`. Reactive (m2) versus homeostatic
(m1) subclusters and AD versus control were compared by Mann–Whitney U and
Cliff's δ; a null of 1,000 size-matched random gene sets gave a permutation
p-value for each contrast.

**Negative-control analysis (§2.5).** MM-GBSA ΔG, Vina/CNN scores, fpocket
descriptors, and the published 100 ns MD core-RMSD values were assembled into
one scorecard (`results/phase7/negctrl/`). No new simulation was run; MD RMSD
values derive from the archived TRUBA trajectories reported in Paper #1
(`gen_seed = −1`).

**Glue generation and filtering (§2.6–2.7).** The reference ligand QFC was read
from PDB 8RQC. The lenalidomide warhead was BRICS-grown with 561 fragments
(from `data/compounds/tier1_cns_approved.csv` + `tier2_approved.csv`) across
five linkers via RDKit `molzip`; products were sanitised, deduplicated by
InChIKey, and embedded (ETKDGv3, MMFF). Descriptors and CNS-MPO (repo 5-check
proxy and approximate Wager 6-parameter) were computed with RDKit. Gates:
strict — MW < 450, cLogP 2–4, TPSA < 90, CNS-MPO(proxy) ≥ 4.0, PAINS-free;
relaxed (pre-registered fallback) — TPSA < 120, CNS-MPO(proxy) ≥ 3.5, others
unchanged. The top 25 relaxed passers were docked with AutoDock Vina into the
8RQC receptor (`data/docking/receptors/IKZF1_8RQC_CRBN_prep.pdbqt`) using the
Paper #1 grid box (centre 0.566, −2.235, 4.760; 20 Å³; exhaustiveness 16),
ligands prepared with Meeko.

---

## 5. Data and code availability

- **Paper #1** (accepted, *J Alzheimers Dis*): code frozen at git tag
  `v1.0-jad-accepted`; SEA-AD atlas at the Allen Institute AWS Open Data
  registry; GSE95587 at GEO.
- **This study:** all code at `code/phase7_robustness/`; result tables,
  figures, and per-analysis `CONCLUSION.md` files at `results/phase7/`;
  GSE138852 at GEO. A versioned Zenodo archive of `results/phase7/` (including
  the generated library SDF and docked poses, which are excluded from git as
  binaries) will be deposited on submission.

---

## 6. Figures and Tables (planned)

- **Figure 1.** Multi-region mixed-effects modelling. (A) *IKZF1* substate
  coefficients (± 95% CI) with `(1|donor) + (1|region)`, small-multiples for
  seven TFs. (B) Leave-one-region-out coefficient stability. (C) NB-GLMM and
  region-as-fixed-effect comparison. Source: `results/phase7/lmm/`.
- **Figure 2.** Donor robustness of the LateAD-DAM peak. (A) Leave-one-donor-out
  homeostatic gap. (B) Downsampling distributions (violin, k = 1–4). (C) ρ
  jackknife and donor bootstrap. Source: `results/phase7/downsample/`.
- **Figure 3.** Independent replication (Grubman 2019). Signature score by
  diagnosis and by microglial subcluster; permutation null. Source:
  `results/phase7/external/`.
- **Figure 4.** MM-GBSA negative controls: ΔG versus 100 ns MD core-RMSD for
  tafamidis, diflunisal, and retained comparators. Source:
  `results/phase7/negctrl/`.
- **Figure 5.** Cereblon glue design. (A) 8RQC ternary interface with warhead
  and exit vector. (B) Generation scheme. (C) 2D grid of top scaffolds.
  (D) Ternary-interface docking scores versus the bare anchor. Source:
  `results/phase7/glue_design/`.
- **Table 1.** Five-block evidence architecture for *IKZF1*
  (`results/phase7/epistemic/evidence_summary.csv`).
- **Table 2.** Top relaxed-filter glue candidates with properties, CNS-MPO, QED,
  and ternary-interface Vina score
  (`results/phase7/glue_design/glue_candidates_top.csv`).

---

## 7. References (key; to be completed)

1. Özkurt Ç. Single-cell gene networks nominate IKZF1 as an Alzheimer's
   microglial regulator. *J Alzheimers Dis* 2026 (accepted; ID ALZ-26-0852).
2. Gazestani V, et al. Early Alzheimer's disease pathology in human cortex …
   (SEA-AD). *Cell* 2023.
3. Grubman A, et al. A single-cell atlas of entorhinal cortex from individuals
   with Alzheimer's disease. *Nat Neurosci* 2019;22:2087–2097. (GSE138852)
4. Ballasch I, et al. [IKZF1 in Alzheimer's disease microglia — western blot].
   2023. *(cited as ref 30 in Paper #1)*
5. Petzold G, et al. Structural basis of lenalidomide-induced CK1α degradation
   by the CRL4^CRBN ubiquitin ligase. *Nature* 2016;532:127–130.
6. Watson ER, et al. Molecular glue CELMoD compounds … IKZF2 degradation
   (PDB 8RQC). *Science* 2024.
7. Wager TT, Hou X, Verhoest PR, Villalobos A. Central nervous system
   multiparameter optimization desirability. *ACS Chem Neurosci* 2016;7:767–775.
8. Degen J, Wegscheid-Gerlach C, Zaliani A, Rarey M. On the art of compiling and
   using 'drug-like' chemical fragment spaces (BRICS). *ChemMedChem*
   2008;3:1503–1507.
9. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0.
   *J Chem Inf Model* 2021;61:3891–3898.
10. Bates D, Mächler M, Bolker B, Walker S. Fitting linear mixed-effects models
    using lme4. *J Stat Softw* 2015;67:1–48.
