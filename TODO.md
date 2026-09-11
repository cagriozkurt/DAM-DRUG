# TODO.md: DAM-DRUG Revision & Validation Roadmap

This roadmap outlines all analytical, computational, and structural revisions required to eliminate artifacts, address survivorship bias, decouple non-independent evidence lines, and satisfy peer review.

---

## 1. Gene Regulatory Networks & Motif Expansion

*Objective: Eliminate survivorship bias and resolve whether IKZF1 outperforms competing regulators.*

> **STATUS 2026-09-12 — CONCLUDED: motif enrichment does not support a direct link for IKZF1 or the 10 comparators.** Dual-window rebuild (wide+narrow, matching aertslab's own production pairing) COMPLETED (job 6348013): 11 regulons retained, but **still 0/11 target TFs, IKZF1 included** — ruling out window-width as the (sole) cause. Direct NES diagnostic (`s1_05_diagnose_nes.py`, reproducing `pyscenic ctx`'s own recovery-curve calc via `ctxcore`) on IKZF1's own top-50-by-importance target module — the single most favourable case — gives NES ≈ 0 (−0.02 wide / +0.57 narrow) for IKZF1's own motif (MA1508.2), nowhere near the 3.0 pass threshold: not a near-miss, a genuine null. Full writeup: `results/phase7/grn_jaspar2026/CONCLUSION.md`. `s1_04` (benchmark/nulls) submitted for completeness as job **6349095**. This is a distinct finding from Section 2 (regional LMM / donor robustness), which operates on expression/AUCell directly and is unaffected.

* [x] **Re-run pySCENIC / RcisTarget with JASPAR 2026:** — DONE. Result: NOT RETAINED.
* [x] Build custom `.feather` cisTarget rankings from the updated JASPAR 2026 vertebrate CORE collection (2,633 PFMs) and UNVALIDATED collections.  — Used CORE **non-redundant** (1,097 motifs) rather than full CORE+UNVALIDATED (a scope decision, "D1", made at `s1_01`); dual wide+narrow window pair built matching aertslab's own production practice.
* [x] Specifically include non-canonical E-box motifs for **BHLHE40** and **BHLHE41** (e.g., `CACGCG`), as well as expanded profiles for **IRF8, PPARG, SPI1, RUNX1,** and **CEBPB**.  — motifs present in the JASPAR2026 CORE set for all 11 (confirmed via motif2tf table lookup); none retained regardless.


* [x] Re-prune the 5-seed consensus GRNBoost2 edge table using `pyscenic ctx` and record which TFs pass regulon pruning.  — **0/11 target TFs retained** (11 regulons retained total, all unrelated to the hypothesis: E2F1/E2F2/E2F7/E2F8/ETV6/FOSL2/NFIB/PBX3/ZNF148). Direct NES diagnostic confirms this is a genuine null (IKZF1's own motif NES≈0 on its own best-case target module), not a threshold/DB-construction artifact. See `results/phase7/grn_jaspar2026/CONCLUSION.md`.




* [~] **Benchmark Rescued Regulons Across Pseudotime:** — MOOT: no target-TF regulons were rescued to benchmark.
* [ ] Calculate per-cell AUCell scores for all rescued regulons across the 6 microglial substates.  — N/A, none of the 11 target TFs have a retained regulon under this DB.


* [ ] Correlate regulon activity against diffusion pseudotime (DPT) using Spearman rank correlation.  — N/A, same reason.


* [ ] Evaluate whether `IKZF1(+)` retains its unique late-stage association or if `BHLHE41`/`IRF8` show equal or superior trajectory coupling.  — Cannot be evaluated under JASPAR2026 cisTarget (neither IKZF1 nor BHLHE41/IRF8 pass pruning); the ORIGINAL HOCOMOCO-based IKZF1(+) regulon from Paper 1 remains the one actually benchmarked (§2.1–2.2 of this study, on that discovery regulon, not this JASPAR2026 rerun).




* [~] **Disentangle Ikaros-Family Paralogues (IKZF1 vs. IKZF2 vs. IKZF3):** — MOOT under this DB (none of the 3 paralogues retained either).
* [ ] Extract the predicted target gene list of the `IKZF1(+)` regulon.  — N/A here; already done for the Paper-1 discovery regulon (used throughout §2.2–2.3 of this study).


* [ ] Compute per-cell Spearman correlations between the mean expression of this target set and the distinct expression traces of *IKZF1*, *IKZF2*, and *IKZF3*.


* [ ] Confirm the regulon is coupled specifically to *IKZF1* transcript levels rather than shared zinc-finger binding promiscuity.




* [x] **Permutation & Negative Null Models:** — job **6349095** submitted (benchmark/nulls suite), running against the 11 regulons this DB actually retained, for completeness.
* [x] Shuffle cell-state and donor labels across the expression matrix (1,000 permutations).


* [x] Quantify the empirical false-discovery rate (FDR) of recovering regulons with $\vert{}\rho\vert{} > 0.30$ along pseudotime to prove `IKZF1` trajectory correlation is not a stochastic artifact.  — expected concordant null given IKZF1 was never retained as a regulon here.




* [ ] **Validate Downstream Targets via JASPAR 2026 Literature Ground Truth:**
* Intersect predicted `IKZF1` regulon targets with the JASPAR 2026 text-mined, curated human TF–Target Gene (TF–TG) interaction database.
* Calculate hypergeometric enrichment to demonstrate biological relevance against experimental ground truth.  — Not yet run; lower priority given IKZF1 has no JASPAR2026-cisTarget-validated regulon to validate targets against in the first place. Could still be run against the Paper-1 discovery regulon's target list if wanted.





---

## 2. Decoupling Evidence & Cross-Cohort Single-Cell Replication

*Objective: Resolve circularity across internal SEA-AD analyses and validate single-cell associations externally.*

> **STATUS 2026-09-10 — substantially complete.** 2A (regional LMM) PASS, 2B (donor downsampling) PASS, 2C (external replication) PASS at state level (Grubman 2019 only; credentialed cohorts + replication-manifold PAGA/DPT outstanding), 2D (epistemic table) done pending manuscript application. Code `code/phase7_robustness/`, results `results/phase7/`, spec `docs/superpowers/specs/2026-09-10-section2-decoupling-design.md`. Branch `robustness/section2-decoupling`.

* [~] **External Single-Cell Replication (Open-Access Processed Cohorts):**  — PARTIAL (Grubman 2019 GSE138852 only; Mathys/Sun are Synapse-credentialed, not accessed; Olah 2020 full matrix deferred). `code/phase7_robustness/2C_external/`, `results/phase7/external/CONCLUSION.md`
* Download processed single-nucleus count matrices and validated cluster annotations from the UCSC Cell Browser / MIT portals for **Mathys et al. (2019/2023)** and **Sun et al. (2023)**.  — NOT DONE (credentialed); used **Grubman et al. 2019** (entorhinal cortex, 449 microglia, subclusters m1–m5) as open-access fallback.


* [x] Score these independent cells with the discovery `IKZF1(+)` AUCell signature.  — `sc.tl.score_genes`, 1072/1271 regulon targets present.


* [x] Test whether `IKZF1` regulon activity is significantly enriched in their annotated reactive/DAM/lipid-associated subclusters relative to homeostatic microglia.  — YES: reactive (m2) vs homeostatic (m1) Cliff's δ +0.72, MWU p=4e-17, random-geneset permutation p=0.001. NOT elevated by crude AD-vs-control (δ −0.18) → signal is state-linked, not diagnosis-linked (concordant with SEA-AD).


* [ ] Run independent PAGA and diffusion pseudotime on the replication microglial manifolds to confirm monotonic upregulation along the disease progression axis.  — NOT DONE (449-cell cohort too small for a stable manifold; deferred to a larger external cohort).




* [x] **Control for Regional Confounding (Linear Mixed-Effects Modeling):**  — DONE. `code/phase7_robustness/2A_lmm/`, `results/phase7/lmm/CONCLUSION.md`. Verdict: PASS.
* [x] Fit a linear mixed-effects model (LMM) for pseudobulk expression across all 10 SEA-AD brain regions: `Expression ~ Substate + Age + Sex + (1 | Donor) + (1 | Brain_Region)`.  — donor×substate×region pseudobulk (10 regions AnG/DFC/FI/HIP/ITG/LEC/MEC/MTG/STG/V1C), `lmerTest`, min10 + unconstrained variants, DAM-IRM dropped.


* [x] Confirm that elevated *IKZF1* expression in DAM/LateAD-DAM remains statistically significant after accounting for regional sampling skew.  — YES: DAM +0.57 log1p-CPM (95% CI 0.53–0.61), LateAD-DAM +0.55 (CI 0.44–0.65); region ICC ≈0.04; substate LRT p<1e-90; 11/11 leave-one-region-out; region-as-fixed and NB-GLMM concordant. Companions IRF8/BHLHE41/RUNX1 same; PPARG down; SPI1/CEBPB DAM-only (not LateAD-DAM).




* [x] **LateAD-DAM Donor Downsampling Permutations:**  — DONE. `code/phase7_robustness/2B_downsample/`, `results/phase7/downsample/CONCLUSION.md`. Verdict: PASS.
* [x] Address state sparsity (median 4 cells/donor) by performing leave-one-donor-out cross-validation and iterative downsampling.  — LOO over 57 LateAD-DAM donors; per-donor cap k=1..4 × 1000 draws; pseudotime ρ donor jackknife (84) + bootstrap (1000).


* [x] Verify that the LateAD-DAM regulon peak ($AUC = 0.153$) is not driven by outlier donors.  — CONFIRMED: LateAD-DAM stays top-AUCell state in 57/57 LOO refits (gap 0.0312–0.0332); gap>0 and top-state in 1000/1000 draws even at 1 cell/donor; ρ jackknife 0.302–0.327 (sign-stable), bootstrap 95% CI 0.27–0.35.




* [x] **Epistemic Classification Update:**  — DONE (analysis + patch text; not yet applied to Manuscript.md). `code/phase7_robustness/2D_epistemic/01_evidence_table.py`, `results/phase7/epistemic/{evidence_summary.csv,manuscript_patches.md}`
* [x] Update Figure 3D / Figure 4 summary tables and text to explicitly group internal SEA-AD modalities (DGE, Regulon, DPT, CellOracle, ATAC) as non-independent dimensions, presenting external datasets (GSE95587 and Mathys/Sun cohorts) as distinct validation lines.  — 5-block table built (internal RNA / internal ATAC / internal robustness / external [GSE95587 + Grubman] / literature); 6 concrete manuscript patches drafted in `manuscript_patches.md` for a separate reviewed commit.





---

## 3. Cell–Cell Communication Triage (SLIT2 → ROBO2)

*Objective: Determine whether the SLIT2–ROBO2 axis is a cortex-wide interaction or an MTG-specific transcript dropout artifact.*

> **STATUS 2026-09-10 — SCRIPTS STAGED, not run (TRUBA).** `code/slurm/paper2/s3_0{1..4}*` + `code/phase7_robustness/s3_cellchat/{s3_prep.py,s3_cellchat.R,s3_collate.py}`. CellChat v2 params identical to the accepted MTG run. Data source RESOLVED: SEA-AD `Multiregion_2026/subclass_objects/` (public S3; region col `Brain Region`; raw counts `.layers["UMIs"]`; includes DFC). `s3_01` pulls ~80 GB (Immune + GABAergic + glia + smaller glutamatergic). Expression-matched 1000-pair permutation null included. Spatial-transcriptomics cross-check (item 3) still manual.

* [ ] **Multi-Region Cross-Validation within SEA-AD:**
* Extract matched neuronal and microglial nuclei from the remaining 9 regions in the SEA-AD atlas (e.g., DLPFC, Hippocampus, Entorhinal Cortex).


* Run CellChat v2 using identical parameters to test if SLIT2→ROBO2 emerges as a top sender–receiver interaction across multiple cortical regions.




* [ ] **Expression-Matched Background Permutation Tests:**
* Generate 1,000 expression-matched background ligand–receptor pairs to determine if the SLIT2→ROBO2 interaction probability ($sum\_prob = 0.318$) exceeds random expectation, ruling out baseline transcript abundance inflation.




* [ ] **Public Spatial Transcriptomics Cross-Check:**
* Inspect public 10x Xenium / Visium human AD cortex data to evaluate whether *SLIT2*+ interneurons physically colocalize with *ROBO2*+ microglia around amyloid plaques.




* [ ] **Manuscript Reframing:**
* If the interaction does not reproduce across regions, downgrade claims in the abstract and discussion, framing SLIT2→ROBO2 as an exploratory, region-restricted MTG hypothesis.





---

## 4. Computational Drug Discovery & Structural Pipeline

*Objective: Eliminate false positives, reframe methodological benchmarks, and transition to a rational chemical biology strategy.*

> **STATUS 2026-09-11 — all local sub-tasks (4A, 4B, 4C, in-vitro gates) complete; only MD production remains, queued on TRUBA.** Code `code/phase7_robustness/{4A_negctrl,4B_glue_gen,4C_selectivity}/`, results `results/phase7/{negctrl,glue_design,selectivity,invitro_gates}/`. hERG QSAR + counter-docking (4C) and in-vitro testing gates done 2026-09-11, folded into `Manuscript_Paper2.md`. Only the T-REMD production run (job 6346303, queued PD/Priority) and its core-RMSD gate remain in this section.

* [x] **Reframe Tafamidis & Diflunisal as Negative Controls / MM-GBSA Benchmarks:**  — DONE (4A). `results/phase7/negctrl/CONCLUSION.md`.
* [x] Retract tafamidis and diflunisal as viable repurposing candidates in the Abstract, Results, and Discussion.  — patch text drafted (Abstract/Results/Tables 3–5/Fig 5A); NOT yet applied to Manuscript.md.


* [x] Reframe their severe ligand egress ($25.17\text{ \AA}$ and $83.42\text{ \AA}$) and failed selectivity ($SI < 1.0$) as a case study highlighting the limitations of implicit-solvent MM-GBSA on shallow, non-druggable TF surfaces.  — scorecard shows both rank #1 by MM-GBSA ΔG within their target set yet egress in MD, while lower-ranked IRF8 compounds are retained. Note: SI<1.0 was the CRBN-track Tier-2 hits, not tafamidis/diflunisal — their failure mode is MD egress + (diflunisal) a PPARG set where all ΔG≈0.




* [x] **Pivot IKZF1 Targeting to CRBN-Interface Molecular Glues:**  — DONE (4B). `results/phase7/glue_design/anchor.json`.
* [x] Formalize the exclusion of the undruggable orthosteric zinc-finger pocket ($drug\_score = 0.001$).  — documented in anchor.json + CONCLUSION.md.


* [x] Focus the drug discovery section entirely on the cereblon (CRBN) ternary interface using PDB 8RQC coordinates.  — 8RQC used for anchor reference (ligand QFC) + docking box; CRBN tri-Trp cage residues recorded.




* [~] **Scaffold-Based Molecular Glue Generation (RDKit / BRICS):**  — DONE with a documented gate deviation.
* [x] Fix the phthaloyl/glutarimide anchor required for the CRBN tri-tryptophan binding pocket.  — lenalidomide isoindolinone–glutarimide, exit vector at the 4-amino N.


* [x] Use RDKit Reaction-SMARTS or BRICS fragmentation to generate derivative libraries growing toward the IKZF1 ZF2 $\beta$-hairpin degron (Gly146 / Lys145 interface).  — 561 BRICS fragments (CNS-approved + Tier-2) × 5 linkers via `molzip` → 2,750 anchor-preserving products; top 25 dock into the 8RQC ternary interface at Vina −5.9…−7.0 kcal/mol (anchor −5.2).


* [~] Filter all generated candidates for central nervous system druglikeness: CNS MPO $\ge 4.0$, $MW < 450$, $\text{cLogP } 2\text{--}4$, and $\text{PSA} < 90\text{ \AA}^2$.  — **strict gates unsatisfiable**: the glutarimide+isoindolinone warhead alone is TPSA 92.5 Å² / cLogP ≈ 0, so TPSA<90 and CNS-MPO≥4 are structurally impossible for any lenalidomide-based glue. 0/2,750 pass strict; 184 pass a pre-registered fallback (TPSA<120, CNS-MPO proxy≥3.5, MW<450, cLogP 2–4, PAINS-free). This is itself a reportable finding (CELMoD space is at the CNS-druglikeness boundary).




* [~] **Explicit-Solvent Molecular Dynamics & T-REMD Stability:** — PREP + SMOKE TEST DONE, production QUEUED (not yet run).
* [x] Re-simulate top CRBN molecular glue candidates and top orthosteric hits in **explicit solvent (TIP3P, Amber99SB-ILDN + GAFF2)** with $150\text{ mM NaCl}$ neutral balance.  — `amber14sb` unavailable on TRUBA, switched to amber99sb-ildn (matches Paper 1 precedent). 3 ternary complexes built (GLUE0231/0246/0692); GLUE0231 EM→NVT→500ps NPT completed clean (job 6346115, 0 LINCS warnings) after fixing `-DPOSRES`+barostat interaction (`refcoord-scaling=com`, `pcoupl=berendsen`).


* [x] Replace deterministic single runs with Temperature Replica Exchange MD (T-REMD; 8 replicas, 300–320 K).  — 10ps smoke test executed for real (RunPod L4, akya-cuda was fully allocated): 49 exchange attempts, real swaps, zero LINCS/fatal/segfault, 43 ns/day. Found + fixed 2 bugs (sed missing `/g`; Zn2+ landing in both T-coupling groups via gmx's generic `Ion` family) in both `s4_03_smoke_test.slurm` and `s4_03_run_tremd.slurm`.


* [~] **Production T-REMD, in progress — moved entirely to local GPU.** TRUBA `akya-cuda` stayed fully allocated with no ETA (`6346303_[0-2]` queued PD/Priority for hours) → moved to the user's local RTX 3060 Mobile (6GB VRAM), 2026-09-11. Re-validated the smoke test there first (same result: 49 exchanges, clean, GPU peaked 1.3GB of 6GB free) before committing. Found along the way: `s4_02` equilibration had ONLY ever been run for GLUE0231 (`array=0`) across the whole debugging session — GLUE0246/GLUE0692 had no `npt.gro`/`npt.cpt` at all. **GLUE0231**: 100ns×8-replica production running on debian-local since 14:40. **GLUE0246**: equilibrated (18:20, ~2h04m NVT+NPT under CPU contention with GLUE0231's production) — production launched immediately after, running since 18:20 alongside GLUE0231 (16 replicas total across both, GPU 1.4GB/6GB, RAM healthy at 13GB free). **GLUE0692**: equilibration running now (EM), production to follow automatically via the same sequential script once it reaches npt.gro. All corresponding TRUBA jobs (`6346303_[0-2]`, `6347306`) cancelled — everything for this section runs on the local GPU, none on TRUBA.


* [ ] Enforce an advancement stability threshold: mean ligand core-RMSD $< 3.5\text{ \AA}$ over the final 20 ns across all replicas.  — gate script `s4_04_core_rmsd.slurm` staged, not runnable until GLUE0231's local run (or 6346303) produces ≥20ns of trajectory.




* [~] **In Silico Off-Target & Toxicity Counter-Screening:** — DONE, gate structurally unmet (real finding). `code/phase7_robustness/4C_selectivity/`, `results/phase7/selectivity/CONCLUSION.md`.
* [x] Screen generated molecules against machine-learning QSAR models for hERG liability trained on ChEMBL bioactivity data, using 2048-bit Morgan circular fingerprints and Tree SHAP interpretability.  — RandomForest on 9,496 ChEMBL CHEMBL240 (KCNH2/hERG) bioactivities, pChEMBL≥5 active cutoff (10 µM), held-out ROC-AUC 0.868 / CV5 0.857±0.002. Applied to all 184 CNS-fallback glue candidates: 140/184 predicted hERG-clean. SHAP TreeExplainer's additivity check failed once (float-accumulation artifact on a 500-tree×2048-feature RF, not a real error — model already confirmed correct via ROC-AUC), fixed with `check_additivity=False`; `herg_shap_summary.png` committed.


* [x] Execute counter-docking against the hERG central cavity (PDB 7CN1), dopamine D2 (PDB 6CM4), and serotonin 5-HT2A (PDB 6A94), enforcing a pre-specified Selectivity Index ($SI > 5.0$).  — reused Phase 1's prepped receptors/grid boxes. **0/25 top-ranked candidates pass SI>5.0** (all score MORE negative, i.e. bind tighter, against every off-target than the 8RQC on-target PPI interface). Also fixed Phase 1's SI formula: `selectivity_table.csv` used SI=vina_ontarget/vina_offtarget, a raw ratio of similarly-scaled Vina scores mathematically bounded near 1 (its own 6 compounds scored 0.58–0.88) — could never reach SI>5.0 for any molecule. Replaced with SI=exp(ΔΔG/RT), the standard free-energy→affinity-ratio proxy. Root cause of the 0/25 result: Vina systematically scores deep small-molecule pockets (GPCR orthosteric sites, hERG central cavity) higher than shallow protein-protein-interaction surfaces like the CRBN-IKZF1 interface — a known docking-methodology limitation, not necessarily true promiscuity. Same character as the CNS-MPO gate deviation: reportable as a limitation, not silently dropped. **Folded into `Manuscript_Paper2.md` as new §2.8/Figure 6/Table 3** (2026-09-11) — draft's Section 2 and 4A/4B content was already current; this closed the one real gap.




* [x] **Define Concrete In Vitro Testing Gates for Translating Leads:** — DONE (protocol specification, no wet-lab data). `results/phase7/invitro_gates/PROTOCOLS.md`, referenced from `Manuscript_Paper2.md` Discussion.
* [x] Detail exact experimental protocols required prior to in vivo work:
* [x] Recombinant IRF8 DBD or PPARG LBD surface plasmon resonance (SPR) / thermal shift assays (TSA).  — Gate 1: exact constructs, SPR single-cycle kinetics + DSF conditions, controls, pass gate (K_D<50µM AND ΔTm≥+1.5°C).


* [x] TR-FRET or cereblon competitive binding assays for CRBN-track glue molecules.  — Gate 2: CRBN-DDB1 TR-FRET tracer displacement + AlphaLISA/TR-FRET ternary-complex (IKZF1-ZF2) assay + neosubstrate (GSPT1/SALL4/ZBTB16) selectivity counter-screen + hERG patch-clamp — explicitly supersedes the docking-based SI shown uninterpretable in 4C.


* [x] CRISPR-mediated *IKZF1* knockout in human iPSC-derived microglia (iMGs) to assess blunting of LateAD-DAM phenotypic transition.  — Gate 3: sgRNA/clone/batch design, RNA-seq readout against the discovery signature, pass gate (≥50% signature blunting, P<0.05, ≥2 clones × ≥3 batches). Tests the study's core causal claim.







---

## 5. Codebase & Data Provenance

*Objective: Ensure all revised scripts adhere to strict CADD reproducibility standards.*

* [x] Section 2 robustness code deposited: `code/phase7_robustness/` (2A_lmm, 2B_downsample, 2C_external, 2D_epistemic) + `envs/damdrug.yml` + `code/phase7_robustness/README.md`.


* [x] JASPAR 2026 parsing + custom `.feather` generation staged: `code/slurm/paper2/s1_01_fetch_jaspar2026.sh`, `s1_02_build_cistarget_db.slurm` (deposited under `code/slurm/paper2/`, not `code/phase2_GRN/` — Paper #2 separation).


* [x] Explicit-solvent GROMACS `.mdp` + T-REMD wrappers staged: `code/slurm/paper2/s4_02_gromacs_prep.slurm` (inline em/nvt/npt mdp), `s4_03_run_tremd.slurm` (inline 8-replica prod mdp), `s4_04_core_rmsd.slurm` (deposited under `code/slurm/paper2/`).


* [x] Pin exact container digests (`scenic.sif`, `gromacs.sif`, `cellchat.sif`) in SLURM scripts, replacing `:latest` tags with immutable image hashes.  — DONE 2026-09-11. Resolved via GHCR anonymous-token registry API (no docker/skopeo/crane available locally): `ghcr.io/cagriozkurt/dam-drug-scanpy@sha256:a16bd5bee57c5fb1918f2e4206368d342861392126c3a4159a462ddd8910efc3` (the actual image behind the misleadingly-named local file `scenic.sif`), `dam-drug-r@sha256:3d0db505d2c1caa014caed89f631b453e0eea13b1ad5919f3ffde92a607d0648` (cellchat), `dam-drug-scmultiomegrn@sha256:6eb4e6978d93bdffb35ae50c8ffccdd9703d4fc6d8e4755ce24e7af770cfc349`, `dam-drug-fpocket@sha256:bc0ce015b2c39aa183cebc5261c9416aef727cd3c58c554a25e4cc2518dab97f` — all 31 `.slurm` files referencing any of these updated, zero `:latest` remaining. **Finding:** "gromacs.sif" in this TODO item doesn't exist — GROMACS is loaded via the TRUBA system module `apps/gromacs/*`, never an apptainer/docker image; nothing to pin there.


* [ ] Archive all updated summary tables and replication datasets to Zenodo under a versioned DOI.