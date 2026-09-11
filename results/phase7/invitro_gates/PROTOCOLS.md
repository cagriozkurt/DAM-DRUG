# Section 4 — Concrete in vitro testing gates for translating leads

**Date:** 2026-09-11
**Status:** protocol definition only — no wet-lab data. This document specifies
exactly what must be run, on what material, with what pass/fail criterion,
before any computational finding in this study (§2.5–2.8) is treated as more
than a hypothesis. Each protocol below is scoped to resolve a *specific*
uncertainty the computational work could not close on its own.

---

## Gate 1 — Recombinant IRF8-DBD / PPARG-LBD biophysical confirmation

**Motivates:** §2.5 (tafamidis/diflunisal negative-control finding). The
computational result (MM-GBSA rank-1 within target set, then MD egress) is a
*prediction* that these compounds do not engage their nominal targets with
meaningful affinity. That prediction has not been tested against real protein.

**Protocol.**
- **Constructs.** IRF8 DNA-binding domain (residues 1–120, His6-SUMO fusion,
  *E. coli* BL21(DE3), Ni-NTA + SUMO-protease cleavage + SEC polish) and PPARG
  ligand-binding domain (residues 206–477, His6-tag, baculovirus/Sf9
  expression, Ni-NTA + SEC), each ≥95% pure by SDS-PAGE, monodisperse by
  analytical SEC.
- **Surface plasmon resonance (SPR).** Immobilize His-tagged protein on a
  Ni-NTA or anti-His-capture chip (Biacore T200 or equivalent) to
  200–400 RU. Flow tafamidis and diflunisal (DMSO stock, 1% final DMSO,
  8-point 2-fold dilution series from 200 µM to 1.6 µM) plus two positive
  controls per target (a literature-validated IRF8-DBD DNA-duplex binder as
  the IRF8 positive control channel; rosiglitazone as the PPARG-LBD positive
  control) and one DMSO-only blank. Single-cycle kinetics, 60 s
  association/120 s dissociation, 25 °C, running buffer 10 mM HEPES pH 7.4,
  150 mM NaCl, 0.05% P20, 2% DMSO-matched.
- **Thermal-shift assay (TSA) as an orthogonal readout.** DSF (Protein
  Thermal Shift dye or SYPRO Orange), 5 µM protein, 50 µM compound (1%
  DMSO), 25→95 °C at 1 °C/min, n=4 technical replicates per condition across
  2 independent protein preparations.

**Pass/fail gate.** A compound is confirmed as a genuine binder only if BOTH
(a) SPR gives a fittable 1:1 or two-state binding isotherm with K_D < 50 µM
and a Hill/steady-state R² > 0.9, AND (b) TSA gives a thermal shift ΔTm ≥
+1.5 °C (mean ± SEM across replicates, two-tailed t-test vs. DMSO P < 0.05).
The pre-registered expectation, consistent with the computational negative-
control finding (§2.5), is that **tafamidis and diflunisal fail this gate**
against IRF8-DBD and PPARG-LBD respectively; a pass would falsify the MM-GBSA/
MD interpretation and require re-opening the repurposing hypothesis.

---

## Gate 2 — TR-FRET / competitive binding for CRBN-track glue candidates

**Motivates:** §2.6–2.8 (glue generation, CNS-gate deviation, and the
uninterpretable docking-based selectivity result). Two independent questions
need real-protein answers: does a candidate actually engage cereblon through
the intended glutarimide pocket, and is engagement selective over the
molecular-glue-relevant liability (loss of normal CRL4^CRBN neosubstrate
recruitment / off-pathway degradation), which docking could not establish.

**Protocol — CRBN engagement (TR-FRET).**
- **Reagents.** Recombinant CRBN–DDB1 heterodimer (His6-CRBN + untagged
  DDB1, baculovirus co-expression, Ni-NTA + SEC), Tb-anti-His antibody
  donor, and a fluorescein- or Cy5-labelled thalidomide/lenalidomide-class
  tracer ligand (commercially available CRBN TR-FRET tracer, e.g. the
  Eurofins/Revvity CRBN TR-FRET kit tracer, or an in-house-synthesized
  fluorescein-glutarimide tracer built on the same anchor used in §2.6).
- **Assay.** 384-well, 10 nM CRBN-DDB1, 5 nM Tb-antibody, 20 nM tracer,
  11-point 3-fold compound dilution series (top 100 µM, 1% DMSO final),
  4 h room-temperature equilibration, TR-FRET ratio (665 nm/620 nm) read on
  a PHERAstar or equivalent. Positive control: unlabelled lenalidomide
  (expected IC50 in the low-µM range for tracer displacement); negative
  control: the des-anchor fragment alone (no glutarimide).
- **Analysis.** 4-parameter logistic fit per compound (n=3 independent
  plates), IC50 with 95% CI.

**Protocol — ternary-complex / neosubstrate-recruitment counter-screen.**
- **AlphaLISA or TR-FRET ternary assay**: recombinant CRBN–DDB1 (donor bead
  or Tb-tag) + recombinant IKZF1 ZF2 domain (residues 141–174,
  acceptor-bead- or fluorophore-tagged). Titrate each Gate-2-positive
  candidate (from the TR-FRET engagement assay above); a true molecular glue
  gives a bell-shaped (hook-effect) ternary-complex signal, distinguishing
  productive ternary-complex formation from simple CRBN occupancy without
  IKZF1 recruitment.
- **Selectivity counter-screen (replaces the uninterpretable docking SI,
  §2.8).** Run the same TR-FRET engagement assay against a panel of the
  native CRL4^CRBN neosubstrates most relevant to on-mechanism toxicity
  (GSPT1, SALL4, ZBTB16 ZF domains, commercially available or in-house
  recombinant), plus a DRD2/5-HT2A radioligand-competition-binding panel
  (e.g. via a commercial CEREP/Eurofins SafetyScreen44-style panel) and a
  manual-patch or automated-patch (e.g. Qube/SyncroPatch) hERG assay in a
  stably transfected HEK293-hERG line.

**Pass/fail gates.**
1. CRBN engagement: tracer-displacement IC50 < 10 µM.
2. Ternary complex: a detectable bell-shaped AlphaLISA/TR-FRET signal with
   IKZF1-ZF2, EC50 < 30 µM at the ternary-signal peak.
3. Neosubstrate selectivity: ≥10-fold weaker engagement (by the same
   TR-FRET readout) against GSPT1/SALL4/ZBTB16 than against IKZF1-ZF2.
4. hERG patch-clamp IC50 > 10 µM (the standard early-safety threshold,
   independent of and superseding the QSAR/docking triage in §2.8, which we
   have already flagged as insufficient on its own).
A candidate advances only if it clears all four; §2.8's QSAR/docking triage
is a pre-filter to prioritise which of the 184 candidates enter this gate
first (QSAR-clean, i.e. the 140/184 predicted hERG-negative subset), not a
substitute for it.

---

## Gate 3 — CRISPR-mediated *IKZF1* knockout in iPSC-derived microglia

**Motivates:** the core discovery hypothesis itself (Paper #1 + §2.1–2.4):
that *IKZF1* is not merely correlated with but causally required for the
LateAD-DAM transition.

**Protocol.**
- **Cell source.** A well-characterized human iPSC line with an existing
  microglial differentiation protocol (e.g. Abud et al. 2017 / McQuade et al.
  2018-style iMGL protocol), from a control (non-carrier) genetic background,
  in an institutionally approved iPSC line already in use or obtainable via
  a repository (e.g. WTSIi series, KOLF2.1J).
- **Knockout.** CRISPR-Cas9 ribonucleoprotein (RNP) nucleofection with 2–3
  independent sgRNAs targeting early *IKZF1* exons (validated cutting
  efficiency by TIDE/ICE sequencing, ≥70% indel rate), single-cell clonal
  expansion, and validation by (a) Sanger/amplicon sequencing of the edited
  locus, (b) Western blot confirming loss of Ikaros protein, and (c) a
  scrambled-sgRNA isogenic control clone carried through the same cloning
  process. Minimum 2 independent KO clones and 2 independent control clones,
  each differentiated to microglia in ≥3 independent differentiation batches.
- **Disease-state induction.** Expose KO and control iMGLs to the
  stimulus/stimuli established (in Paper #1 or the field) to drive a
  DAM-like/LateAD-DAM-like transcriptional shift in vitro (e.g. amyloid-β
  oligomers + apoptotic-neuron debris co-culture, or the specific stimulus
  cocktail used for the in-silico CellOracle perturbation in Paper #1, for
  direct comparability).
- **Readout.** Bulk RNA-seq (minimum n=3 biological replicates per
  genotype × condition) scored against the discovery LateAD-DAM/DAM
  pseudobulk signature (the same gene signature used in §2.1 and the
  Grubman projection in §2.3) via `sc.tl.score_genes` or GSEA; orthogonal
  validation by qPCR of the top 10 LateAD-DAM marker genes and by flow
  cytometry / immunocytochemistry for 2–3 protein-level DAM markers.

**Pass/fail gate.** *IKZF1*-KO iMGLs must show a **≥50% blunting** of the
LateAD-DAM/DAM signature score shift (KO-stimulated vs. KO-unstimulated,
relative to control-stimulated vs. control-unstimulated), consistent across
≥2 independent KO clones and reproducible in ≥3 independent differentiation
batches (mixed-effects model, clone and batch as random effects, genotype ×
stimulus interaction P < 0.05). This is the single experiment that would
most directly test — and could falsify — the central causal claim of Paper #1
and this study.

---

## Summary table

| Gate | Question resolved | Pass criterion | Blocks |
|---|---|---|---|
| 1. IRF8-DBD / PPARG-LBD SPR+TSA | Do tafamidis/diflunisal really bind? | K_D<50 µM (SPR) AND ΔTm≥+1.5°C (TSA) | Confirms/refutes §2.5 |
| 2. CRBN TR-FRET + ternary + selectivity panel | Do glue candidates engage CRBN/IKZF1 selectively? | IC50<10 µM; ternary EC50<30 µM; ≥10x neosubstrate selectivity; hERG patch IC50>10 µM | Confirms/refutes §2.6–2.8; supersedes the docking-based SI |
| 3. CRISPR *IKZF1* KO in iMGLs | Is IKZF1 causally required for LateAD-DAM? | ≥50% signature blunting, P<0.05, ≥2 clones × ≥3 batches | Tests the core hypothesis (Paper #1 + §2.1–2.4) |

None of these experiments have been run. This document exists so that any
future wet-lab collaborator (or reviewer) has an unambiguous, falsifiable
specification rather than a vague "further validation is needed" statement.
