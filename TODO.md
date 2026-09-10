# TODO.md: DAM-DRUG Revision & Validation Roadmap

This roadmap outlines all analytical, computational, and structural revisions required to eliminate artifacts, address survivorship bias, decouple non-independent evidence lines, and satisfy peer review.

---

## 1. Gene Regulatory Networks & Motif Expansion

*Objective: Eliminate survivorship bias and resolve whether IKZF1 outperforms competing regulators.*

> **STATUS 2026-09-10 — SCRIPTS STAGED, not run (TRUBA).** `code/slurm/paper2/s1_0{1..4}*` + `code/phase7_robustness/s1_grn/s1_benchmark.py`. Decisions: JASPAR 2026 CORE vertebrate only, genome-wide cisTarget DB. FIXMEs (JASPAR URL, hg38 region BED, curated TF-TG set) in `code/slurm/paper2/README.md`. Covers all 5 bullets below (permutation FDR = 1000-perm unrestricted + donor-block; paralogue test; hypergeometric).

* [ ] **Re-run pySCENIC / RcisTarget with JASPAR 2026:**
* Build custom `.feather` cisTarget rankings from the updated JASPAR 2026 vertebrate CORE collection (2,633 PFMs) and UNVALIDATED collections.
* Specifically include non-canonical E-box motifs for **BHLHE40** and **BHLHE41** (e.g., `CACGCG`), as well as expanded profiles for **IRF8, PPARG, SPI1, RUNX1,** and **CEBPB**.


* Re-prune the 5-seed consensus GRNBoost2 edge table using `pyscenic ctx` and record which TFs pass regulon pruning.




* [ ] **Benchmark Rescued Regulons Across Pseudotime:**
* Calculate per-cell AUCell scores for all rescued regulons across the 6 microglial substates.


* Correlate regulon activity against diffusion pseudotime (DPT) using Spearman rank correlation.


* Evaluate whether `IKZF1(+)` retains its unique late-stage association or if `BHLHE41`/`IRF8` show equal or superior trajectory coupling.




* [ ] **Disentangle Ikaros-Family Paralogues (IKZF1 vs. IKZF2 vs. IKZF3):**
* Extract the predicted target gene list of the `IKZF1(+)` regulon.


* Compute per-cell Spearman correlations between the mean expression of this target set and the distinct expression traces of *IKZF1*, *IKZF2*, and *IKZF3*.


* Confirm the regulon is coupled specifically to *IKZF1* transcript levels rather than shared zinc-finger binding promiscuity.




* [ ] **Permutation & Negative Null Models:**
* Shuffle cell-state and donor labels across the expression matrix (1,000 permutations).


* Quantify the empirical false-discovery rate (FDR) of recovering regulons with $\vert{}\rho\vert{} > 0.30$ along pseudotime to prove `IKZF1` trajectory correlation is not a stochastic artifact.




* [ ] **Validate Downstream Targets via JASPAR 2026 Literature Ground Truth:**
* Intersect predicted `IKZF1` regulon targets with the JASPAR 2026 text-mined, curated human TF–Target Gene (TF–TG) interaction database.
* Calculate hypergeometric enrichment to demonstrate biological relevance against experimental ground truth.





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

> **STATUS 2026-09-10 — local sub-tasks (4A + 4B) complete.** Branch `robustness/section4-drug`. Code `code/phase7_robustness/{4A_negctrl,4B_glue_gen}/`, results `results/phase7/{negctrl,glue_design}/`, spec `docs/superpowers/specs/2026-09-10-section4-drug-reframe-glue-gen-design.md`. Explicit-solvent MD/T-REMD, hERG QSAR, counter-docking, in-vitro gates deferred (TRUBA / later pass). Manuscript patches drafted, not applied.

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




* [ ] **Explicit-Solvent Molecular Dynamics & T-REMD Stability:**
* Re-simulate top CRBN molecular glue candidates and top orthosteric hits in **explicit solvent (TIP3P, Amber14SB + GAFF2)** with $150\text{ mM NaCl}$ neutral balance.


* Replace deterministic single runs with Temperature Replica Exchange MD (T-REMD; 8 replicas, 300–320 K).


* Enforce an advancement stability threshold: mean ligand core-RMSD $< 3.5\text{ \AA}$ over the final 20 ns across all replicas.




* [ ] **In Silico Off-Target & Toxicity Counter-Screening:**
* Screen generated molecules against machine-learning QSAR models for hERG liability trained on ChEMBL bioactivity data, using 2048-bit Morgan circular fingerprints and Tree SHAP interpretability.


* Execute counter-docking against the hERG central cavity (PDB 7CN1), dopamine D2 (PDB 6CM4), and serotonin 5-HT2A (PDB 6A94), enforcing a pre-specified Selectivity Index ($SI > 5.0$).




* [ ] **Define Concrete In Vitro Testing Gates for Translating Leads:**
* Detail exact experimental protocols required prior to in vivo work:
* Recombinant IRF8 DBD or PPARG LBD surface plasmon resonance (SPR) / thermal shift assays (TSA).


* TR-FRET or cereblon competitive binding assays for CRBN-track glue molecules.


* CRISPR-mediated *IKZF1* knockout in human iPSC-derived microglia (iMGs) to assess blunting of LateAD-DAM phenotypic transition.







---

## 5. Codebase & Data Provenance

*Objective: Ensure all revised scripts adhere to strict CADD reproducibility standards.*

* [x] Section 2 robustness code deposited: `code/phase7_robustness/` (2A_lmm, 2B_downsample, 2C_external, 2D_epistemic) + `envs/damdrug.yml` + `code/phase7_robustness/README.md`.


* [x] JASPAR 2026 parsing + custom `.feather` generation staged: `code/slurm/paper2/s1_01_fetch_jaspar2026.sh`, `s1_02_build_cistarget_db.slurm` (deposited under `code/slurm/paper2/`, not `code/phase2_GRN/` — Paper #2 separation).


* [x] Explicit-solvent GROMACS `.mdp` + T-REMD wrappers staged: `code/slurm/paper2/s4_02_gromacs_prep.slurm` (inline em/nvt/npt mdp), `s4_03_run_tremd.slurm` (inline 8-replica prod mdp), `s4_04_core_rmsd.slurm` (deposited under `code/slurm/paper2/`).


* [ ] Pin exact container digests (`scenic.sif`, `gromacs.sif`, `cellchat.sif`) in SLURM scripts, replacing `:latest` tags with immutable image hashes.


* [ ] Archive all updated summary tables and replication datasets to Zenodo under a versioned DOI.