# Paper #2 — TRUBA HPC scripts

Cluster-bound work for the robustness & de novo glue study (branch
`paper2-robustness-glue`). Authored locally, submitted on TRUBA. The accepted
JAD paper's `code/slurm/01–36` are frozen; new work is numbered `s1_*`, `s3_*`,
`s4_*` here.

Spec: `docs/superpowers/specs/2026-09-10-truba-staging-paper2.md`
(decisions D1–D6 recorded there).

```
export DAM_DRUG_DIR=/arf/scratch/mozkurt/DAM-DRUG
export SBATCH_ACCOUNT=<your_account>          # not set in non-interactive shells
cd $DAM_DRUG_DIR
```

**TRUBA facts (verified 2026-09-10 via `ssh truba`, node arf-ui2):**
- Project dir `/arf/scratch/mozkurt/DAM-DRUG` is a plain copy, **not a git
  checkout** — sync new scripts with `rsync`/`scp` from the local repo.
- Containers present: `containers/{scenic.sif,dam-drug-r.sif,dam-drug-scmultiomegrn.sif,fpocket-env.sif}`
  (no gromacs.sif — MD uses the module).
- GROMACS modules: `apps/gromacs/2024.1-oneapi2024` (CPU, oneAPI — used by the
  accepted paper), `apps/gromacs/2023.2-cuda` (**use this for the GPU T-REMD**),
  `apps/gromacs/2023.3`.
- `rclone` is at `/usr/bin/rclone` (handy for the large S3 pulls); `aws` CLI is
  not installed. UI node has outbound internet (curl/wget work).
- `data/raw/SEA-AD/` currently holds only the microglia + MTG-RNA objects —
  every WP needs its downloads run first.

Analysis code these wrappers call: `code/phase7_robustness/{s1_grn,s3_cellchat,s4_md}/`.
Outputs: `results/phase7/{grn_jaspar2026,cellchat_multiregion,glue_md}/`.

---

## WP1 — GRN overhaul (pySCENIC + JASPAR 2026 CORE vertebrate, genome-wide DB)

| order | script | partition | ~time | notes |
|---|---|---|---|---|
| 1 | `s1_01_fetch_jaspar2026.sh` | login node | mins | JASPAR 2026 CORE vertebrate PFMs + motif2tf; **check FIXME URLs** |
| 2 | `s1_02_build_cistarget_db.slurm` | barbun 20c 128G | **1–3 days** | genome-wide `create_cisTarget_databases`; needs `cbust`, hg38.fa, region BED (**FIXME**) |
| 3 | `s1_03_pyscenic_ctx_jaspar.slurm` | barbun 20c 180G | 6–24 h | `pyscenic ctx` + `aucell` with the new feather |
| 4 | `s1_04_benchmark_nulls.slurm` | barbun 20c 128G | 2–6 h | AUCell×substate, pseudotime ρ, IKZF1 vs BHLHE41/IRF8 ranking, 1,000-permutation FDR, IKZF1/2/3 paralogue test, JASPAR TF–TG hypergeometric |

## WP2 — Multi-region CellChat v2 (SLIT2 → ROBO2 triage)

**Data source resolved (2026-09-10):** SEA-AD "Multiregion 2026" release,
`s3://sea-ad-single-cell-profiling/Multiregion_2026/subclass_objects/`
(public, `--no-sign-request`). One h5ad per **subclass**, each spanning all 10
regions; region column `Brain Region`; raw counts in `.layers["UMIs"]`.
The per-region folders exist but `DFC/RNAseq/` is **empty** — the subclass
objects are the reliable source and include DFC. `s3_01` fetches the
receiver (Immune) + all GABAergic senders + glia/vascular + smaller
glutamatergic context (~80 GB); the 3 giant IT classes are optional
(`sbatch s3_01_download_subclass_objects.slurm all`).

| order | script | partition | ~time | notes |
|---|---|---|---|---|
| 1 | `s3_01_download_subclass_objects.slurm` | login/transfer | hours | ~80 GB (required set); `all` arg adds L23/L4/L5/L6-IT (~+90 GB) |
| 2 | `s3_02_prep_regions.slurm` | barbun 20c 180G, array 0–8 | 2–6 h | per region: subset each subclass obj by `Brain Region`, `.layers["UMIs"]`, coarsen, 5k/type cap, concat → prep/ |
| 3 | `s3_03_cellchat_array.slurm` | barbun 40c 350G, array 0–8 | 4–12 h/region | CellChat v2, identical params to `code/phase2_LR/02_cellchat_nichechat.R` |
| 4 | `s3_04_collate_permute.slurm` | barbun 20c 64G | 1–2 h | SLIT2→ROBO2 rank + sum_prob per region; 1,000 expr-matched background L–R pairs → permutation p; pan-cortical vs MTG-restricted verdict |

## WP3 — T-REMD of glue ternary complexes (3 top glues)

| order | script | partition | ~time | notes |
|---|---|---|---|---|
| 0 | (rerun `code/phase7_robustness/4B_glue_gen/06_dock_glues.py` on TRUBA) | — | mins | regenerate `docked/GLUE*.pdbqt` (not in git) |
| 1 | `code/phase7_robustness/s4_md/s4_01_build_ternary.py` | login/barbun | mins | 3 ternary complexes: **one copy** (chains A CRBN + B IKZF1-ZF2 + 2 Zn), glue docked pose, GAFF2 via acpype |
| 2 | `s4_02_gromacs_prep.slurm` (calls `s4_md/s4_zn_restraints.py`) | barbun 20c 32G, array 0–2 | 1–3 h | TIP3P dodecahedron, 150 mM NaCl, Amber14SB+GAFF2; **Zn²⁺ harmonic distance restraints** (auto-detected coordination, `[ intermolecular_interactions ]` funct 6, k=10000 kJ/mol/nm², r0 0.23 nm Zn–S(Cys) / 0.20 nm Zn–N(His)); EM→NVT→NPT |
| 3 | `s4_03_run_tremd.slurm` | akya-cuda, array 0–2 (×8 replicas) | 1–3 days/complex | T-REMD 8 replicas 300–320 K, GROMACS `-replex`; checkpoint-resume |
| 4 | `s4_04_core_rmsd.slurm` | barbun 20c 32G | <1 h | ligand core-RMSD on 10% lowest-RMSF residues, last 20 ns, per replica; **PASS if mean < 3.5 Å across all 8** |

---

## Known FIXMEs (fill on TRUBA)

- `s1_01`: JASPAR 2026 bulk-download URL (check `https://jaspar.elixir.no/downloads/`).
- `s1_02`: hg38 gene-region BED matching `10kbp_up_10kbp_down_full_tx` (aertslab
  `create_cisTarget_databases` wiki), path to `cbust` binary, hg38.fa location.
- ~~`s3_01`: SEA-AD S3 object keys~~ **RESOLVED** — `Multiregion_2026/subclass_objects/`
  (see the WP2 section above). Confirm subclass names against a fresh
  `aws s3 ls s3://sea-ad-single-cell-profiling/Multiregion_2026/subclass_objects/ --no-sign-request`
  before submitting (taxonomy may add types).
- ~~`s4_02`: Zn model choice~~ **RESOLVED** — harmonic distance restraints
  (`s4_md/s4_zn_restraints.py`), k=10000 kJ/mol/nm². Verified on the 8RQC
  structure: CRBN Zn → 4×Cys-SG (C4 site); IKZF1 ZF2 Zn → 2×Cys-SG + 2×His-NE2
  (C2H2 site). Actual construct numbering is Cys147/150, His163/167 (auto-detected).
- All: replace container `:latest` with `@sha256:` digests at Zenodo deposit
  (TODO §5).
