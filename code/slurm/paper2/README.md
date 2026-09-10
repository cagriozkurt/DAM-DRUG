# Paper #2 — TRUBA HPC scripts

Cluster-bound work for the robustness & de novo glue study (branch
`paper2-robustness-glue`). Authored locally, submitted on TRUBA. The accepted
JAD paper's `code/slurm/01–36` are frozen; new work is numbered `s1_*`, `s3_*`,
`s4_*` here.

Spec: `docs/superpowers/specs/2026-09-10-truba-staging-paper2.md`
(decisions D1–D6 recorded there).

```
export DAM_DRUG_DIR=/arf/scratch/mozkurt/DAM-DRUG
export SBATCH_ACCOUNT=<your_account>
cd $DAM_DRUG_DIR
```

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

| order | script | partition | ~time | notes |
|---|---|---|---|---|
| 1 | `s3_01_download_regions.slurm` | login/transfer | hours | **FIXME**: Allen AWS S3 URLs for the 9 non-MTG SEA-AD objects (neurons + microglia) |
| 2 | `s3_02_prep_regions.slurm` | barbun 20c 180G | 2–6 h | per region → `counts_raw.h5` / `gene_names.csv` / `cell_meta.csv` |
| 3 | `s3_03_cellchat_array.slurm` | barbun 40c 350G, array 0–8 | 4–12 h/region | CellChat v2, identical params to `code/phase2_LR/02_cellchat_nichechat.R` |
| 4 | `s3_04_collate_permute.slurm` | barbun 20c 64G | 1–2 h | SLIT2→ROBO2 rank + sum_prob per region; 1,000 expr-matched background L–R pairs → permutation p; pan-cortical vs MTG-restricted verdict |

## WP3 — T-REMD of glue ternary complexes (3 top glues)

| order | script | partition | ~time | notes |
|---|---|---|---|---|
| 0 | (rerun `code/phase7_robustness/4B_glue_gen/06_dock_glues.py` on TRUBA) | — | mins | regenerate `docked/GLUE*.pdbqt` (not in git) |
| 1 | `code/phase7_robustness/s4_md/s4_01_build_ternary.py` | login/barbun | mins | 3 ternary complexes: glue pose into `8RQC_CRBN_ZF2.pdb`, keep 4 Zn, GAFF2 via acpype |
| 2 | `s4_02_gromacs_prep.slurm` | barbun 20c 32G, array 0–2 | 1–3 h | TIP3P dodecahedron, 150 mM NaCl, Amber14SB+GAFF2, Zn restraints, EM→NVT→NPT |
| 3 | `s4_03_run_tremd.slurm` | akya-cuda, array 0–2 (×8 replicas) | 1–3 days/complex | T-REMD 8 replicas 300–320 K, GROMACS `-replex`; checkpoint-resume |
| 4 | `s4_04_core_rmsd.slurm` | barbun 20c 32G | <1 h | ligand core-RMSD on 10% lowest-RMSF residues, last 20 ns, per replica; **PASS if mean < 3.5 Å across all 8** |

---

## Known FIXMEs (fill on TRUBA)

- `s1_01`: JASPAR 2026 bulk-download URL (check `https://jaspar.elixir.no/downloads/`).
- `s1_02`: hg38 gene-region BED matching `10kbp_up_10kbp_down_full_tx` (aertslab
  `create_cisTarget_databases` wiki), path to `cbust` binary, hg38.fa location.
- `s3_01`: Allen SEA-AD AWS S3 object keys for AnG, DFC, FI, HIP, ITG, LEC, MEC,
  STG, V1C (neurons + microglia; see `https://registry.opendata.aws/allen-sea-ad-atlas`).
- `s4_02`: Zn model choice — bonded/cationic-dummy vs restrained.
- All: replace container `:latest` with `@sha256:` digests at Zenodo deposit
  (TODO §5).
