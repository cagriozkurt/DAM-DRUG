# TRUBA HPC Staging for Paper #2 — Design Spec

**Date:** 2026-09-10
**Branch:** `paper2-robustness-glue`
**Scope:** SLURM + prep scripts for the three cluster-bound work packages
(TODO.md Sections 1, 3, and the Section 4 MD tail). Scripts are authored here
and submitted on TRUBA by the user; nothing runs locally.
**New location:** `code/slurm/paper2/` (keeps the accepted-paper `code/slurm/01–36`
numbering frozen).

## Conventions (from existing `code/slurm/`)

- Header: `#SBATCH --output=logs/%x-%j.out`, account via `$SBATCH_ACCOUNT`,
  `set -eo pipefail`, `PROJDIR=${DAM_DRUG_DIR:-$(pwd)}`.
- CPU partition `barbun`; GPU partition `akya-cuda`.
- Containers in `$PROJDIR/containers/*.sif`, `apptainer exec --bind $PROJDIR:$PROJDIR`.
  SCENIC image runs `conda run -n scenic`; R image uses `Rscript --vanilla`.
- **Container digests:** the reproducibility audit wants `:latest` replaced with
  immutable `@sha256:` digests. Decision D3 below.

---

## Work Package 1 — GRN overhaul (pySCENIC + JASPAR 2026)

*Objective (TODO §1): rebuild cisTarget rankings from JASPAR 2026, re-prune the
5-seed GRNBoost2 consensus, and benchmark whether IKZF1 keeps its unique
late-stage trajectory coupling once BHLHE40/41 non-canonical E-boxes and
expanded IRF8/PPARG/SPI1/RUNX1/CEBPB profiles are represented.*

Allocation: 1× `barbun`, 20 cores, 64–128 GB (DB build); ctx/aucell as in
`08_run_ctx_aucell.slurm`.

Inputs already on TRUBA: `containers/scenic.sif`, the 5 per-seed adjacency
matrices + `adj_matrix_aggregated.tsv`, `microglia_raw.loom`,
`scenic_auc_aggregated.loom` (old), `results/phase1/trajectory/pseudotime.csv`.

### Scripts
| file | does | notes |
|---|---|---|
| `s1_01_fetch_jaspar2026.sh` | JASPAR 2026 CORE vertebrate **non-redundant** (1,019 matrices; transfac + jaspar formats, URLs verified) → per-motif transfac files + `jaspar2026_motif2tf.tbl` (dimers split → 1,096 rows) | D1 = CORE vertebrate only |
| `s1_02_build_cistarget_db.slurm` | `cbust` (aertslab binary + source-build fallback), hg38.fa, the **canonical aertslab v10 region BED** (`hg38-limited-upstream10000-tss-downstream10000-full-transcript.bed`, 92,636 GENE#N segments — exact v10 region definition), pure-python FASTA extract (no bedtools on TRUBA, tested), transfac→cbust `.cb` (tested), then `create_cistarget_motif_databases.py` via `conda run -n scenic` | heaviest step (1–3 d). D2 = genome-wide. No methodological deviation from v10 in the region set now. |
| `s1_03_pyscenic_ctx_jaspar.slurm` | `pyscenic ctx` on `adj_matrix_aggregated.tsv` with the new feather + a JASPAR motif2tf table; then `pyscenic aucell` on `microglia_raw.loom` | mirrors `08_run_ctx_aucell.slurm` |
| `s1_00_fetch_curated_targets.sh` | ChEA3 Literature+ENCODE+ReMap ChIP-seq ∪ DoRothEA A/B/C (OmniPath TSV) → `data/references/curated_ikzf1_targets.txt` (3,238 symbols, verified 2026-09-10) | login node, no R |
| `s1_04_benchmark_nulls.slurm` + `s1_benchmark.py` | (a) per-cell AUCell for every rescued regulon across the 6 substates; (b) Spearman vs diffusion pseudotime; (c) rank IKZF1 vs BHLHE41/IRF8/etc.; (d) 1,000-permutation pseudotime shuffle (unrestricted + donor-block) → empirical FDR for \|ρ\|>0.30; (e) IKZF1 vs IKZF2 vs IKZF3 paralogue test; (f) hypergeometric enrichment of the rescued IKZF1 regulon targets vs `curated_ikzf1_targets.txt` | pure Python, `scenic.sif` |

Outputs → `results/phase7/grn_jaspar2026/`.

---

## Work Package 2 — Multi-region CellChat v2 (SLIT2 → ROBO2 triage)

*Objective (TODO §3): test whether SLIT2→ROBO2 (interneuron → microglia) is
MTG-specific or pan-cortical, with an expression-matched permutation null.*

Allocation: 1× `barbun`, 20–40 cores, 128–350 GB; array over regions.

**Data gap:** only `SEAAD_MTG_*` and the microglia-only multi-regional h5ad are
on the SSD. CellChat needs **neurons + microglia** from the other 9 regions.
SEA-AD releases these as per-region or per-class objects on the Allen AWS
bucket — **exact object names / URLs are unknown to me** (D4).

### Scripts
| file | does | notes |
|---|---|---|
| `s3_01_download_regions.slurm` | pull the 9 non-MTG region h5ads to `data/raw/SEA-AD/regions/` | modelled on `01_download_seaad.slurm`; **D4** URLs |
| `s3_02_prep_regions.slurm` + `s3_prep.py` | per region: subset to neuronal + microglial supertypes, extract raw counts + metadata to `counts_raw.h5` / `gene_names.csv` / `cell_meta.csv` (same layout `01_prep_mtg_for_cellchat.py` produces), seed 42 subsample matching the MTG run | reuse existing prep logic |
| `s3_03_cellchat_array.slurm` + `s3_cellchat.R` | array 0–8: CellChat v2 with the **identical** parameters to `02_cellchat_nichechat.R`; save the per-region CellChat object + the ranked L–R interaction table | `dam-drug-r.sif` |
| `s3_04_collate_permute.slurm` + `s3_collate.py` | (a) is SLIT2→ROBO2 a top sender–receiver pair in each region? rank + sum_prob per region; (b) 1,000 expression-matched background L–R pairs per region → permutation p for the observed `sum_prob = 0.318`; (c) verdict table: pan-cortical vs MTG-restricted | drives the manuscript reframing decision |

Outputs → `results/phase7/cellchat_multiregion/`.
Spatial transcriptomics cross-check (TODO §3 item 3) is **out of scope** here —
separate manual task.

---

## Work Package 3 — Explicit-solvent MD / T-REMD of glue ternary complexes

*Objective (TODO §4): kinetic stability of the top generated glues in the
CRBN–IKZF1(ZF2) ternary complex; advancement gate mean core-RMSD < 3.5 Å over
the final 20 ns.*

Allocation: `akya-cuda` GPU, array; T-REMD = 8 replicas 300–320 K per complex.

Inputs: `data/structures/pdb/8RQC_CRBN_ZF2.pdb` (receptor + ZF2 + 4 Zn),
`results/phase7/glue_design/glue_top.sdf`,
`results/phase7/glue_design/docked/GLUE*.pdbqt` (top poses — **not in git**,
regenerate with `06_dock_glues.py` on TRUBA first), `containers/` GROMACS
(system module `apps/gromacs/2024.1-oneapi2024` per prior runs).

### Scripts
| file | does | notes |
|---|---|---|
| `s4_01_build_ternary.py` (local or TRUBA) | for the top 2–3 glues: place the docked glue pose into `8RQC_CRBN_ZF2.pdb`, keep the 4 Zn, protonate (pH 7.4), GAFF2-parameterise the ligand (acpype), assemble the complex | RDKit + parmed/acpype; **D5**: how many complexes (2 or 3) |
| `s4_02_gromacs_prep.slurm` | solvate TIP3P dodecahedron, 150 mM NaCl neutralising, Amber99SB-ILDN + GAFF2, EM → NVT → NPT | **Zn RESOLVED** — harmonic distance restraints via `s4_md/s4_zn_restraints.py` (`[ intermolecular_interactions ]` funct 6, k=10000 kJ/mol/nm², r0 0.23 nm Zn–S / 0.20 nm Zn–N; coordination auto-detected). Verified on 8RQC: CRBN C4 + IKZF1-ZF2 C2H2. |
| `s4_03_run.slurm` | **D6**: (a) plain 100 ns × 1 per complex (`akya-cuda`, ~1–2 days) *or* (b) full T-REMD, 8 replicas 300–320 K, REST/`plumed` or GROMACS `-replex` | array; checkpoint-resume like `24_run_md100ns.slurm` |
| `s4_04_core_rmsd.slurm` + `s4_rmsd.py` | ligand core-RMSD on the 10 % lowest-RMSF residues, last 20 ns, per replica; PASS if mean < 3.5 Å across all replicas | mirrors `25_core_rmsd_analysis.slurm` |

Outputs → `results/phase7/glue_md/`.

---

## Decisions (resolved 2026-09-10)

- **D1 — JASPAR 2026 collections:** **CORE vertebrate, non-redundant** — 1,019
  matrices (the "2,633" figure is CORE across all taxa / all formats). URLs
  verified 2026-09-10. All target TFs present (IKZF1 MA1508.2, BHLHE40
  MA0464.3, BHLHE41 MA0636.1, IRF8, SPI1, RUNX1, CEBPB, PPARG) — this resolves
  the manuscript's BHLHE40/41 atypical-E-box coverage gap.
- **D2 — cisTarget DB scope:** **full genome-wide** `create_cisTarget_databases`
  over the **canonical aertslab v10 region set**
  (`hg38-limited-upstream10000-tss-downstream10000-full-transcript.bed`, 92,636
  GENE#N segments — the exact definition behind
  `hg38_10kbp_up_10kbp_down_full_tx_v10_clust`). Full methodological symmetry
  with the accepted paper's DB. `cbust`, hg38.fa, the BED, and FASTA extraction
  are all handled in `s1_02` — no remaining FIXME.
- **D3 — container digests:** keep `:latest`, add a Section-5 TODO to pin
  `@sha256:` at deposit time (default).
- **D4 — SEA-AD non-MTG data — RESOLVED (2026-09-10):** SEA-AD "Multiregion
  2026" release, `s3://sea-ad-single-cell-profiling/Multiregion_2026/subclass_objects/`
  (public). One h5ad per subclass, all 10 regions each; region col `Brain
  Region`; raw counts in `.layers["UMIs"]`. `DFC/RNAseq/` folder is empty but
  DFC is present in the subclass objects. `s3_01_download_subclass_objects.slurm`
  fetches the Immune + GABAergic + glia + smaller-glutamatergic set (~80 GB;
  `all` arg adds the 3 giant IT classes).
- **D5 — number of ternary complexes for MD:** **3** top glues (default).
- **D-Zn — RESOLVED:** harmonic distance restraints (not cationic-dummy).
  `s4_md/s4_zn_restraints.py` auto-detects each Zn's coordinators from the
  structure and writes `[ intermolecular_interactions ]` funct-6 bonds
  (k = 10000 kJ/mol/nm²; r0 = 0.230 nm Zn–S(Cys), 0.200 nm Zn–N(His)).
  Both the CRBN Zn (4×Cys, C4) and the IKZF1 ZF2 Zn (2×Cys + 2×His, C2H2) are
  restrained. Only one ternary copy (8RQC chains A + B + 2 Zn) is built.
- **D6 — MD protocol:** **full T-REMD**, 8 replicas, 300–320 K per complex;
  advancement gate = mean ligand core-RMSD < 3.5 Å over the final 20 ns across
  all replicas.

## Deliverables
`code/slurm/paper2/` (12–14 scripts + `README.md` with submit order and the
D1–D6 choices recorded), analysis `.py`/`.R` alongside or under
`code/phase7_robustness/{s1_grn,s3_cellchat,s4_md}/`.
