# DAM-DRUG

**Integrative Single-Cell Analysis of Alzheimer's Disease Microglia Identifies IKZF1 as a Late-Stage Neuroinflammatory Regulator and Prioritises Tafamidis and Diflunisal as Candidate Repurposable Therapeutics**

Çağrı Özkurt, Pelin Kelicen-Uğur — Hacettepe University

> Özkurt Ç, Kelicen-Uğur P. (manuscript in preparation).

---

## Overview

Starting from 236,002 microglial nuclei across 84 donors (SEA-AD atlas), this pipeline spans six phases:

1. **QC & substate annotation** — scVI batch correction; five ordered substates (HM → IRM → DAM → LAM → LateAD-DAM)
2. **GRN inference** — pySCENIC 5-seed GRNBoost2 consensus (1.46 M edges; 46 validated regulons); scMultiomeGRN multi-omic cross-validation
3. **Cell–cell communication** — CellChat (3,718 L-R pairs; 15 MTG cell types)
4. **Structure-based virtual screening** — AutoDock Vina + GNINA CNN rescoring + MM-GBSA (1,677 FDA-approved CNS drugs; 6 TF targets)
5. **Multi-modal validation** — ATAC chromatin enrichment; pseudobulk DESeq2; bulk replication (GSE95587); CellOracle in-silico KO
6. **Figures & tables** — all manuscript outputs

**Key findings:**
- IKZF1 is the primary late-disease microglial transcriptional amplifier, validated across six independent modalities
- Tafamidis (→IRF8, ΔG = −9.5 kcal/mol) and diflunisal (→PPARG, ΔG = −2.8 kcal/mol) are the top repurposing candidates
- SLIT2→ROBO2 from inhibitory interneurons is the dominant extrinsic signal to AD microglia

---

## Repository structure

```
DAM-DRUG/
├── code/
│   ├── phase1_QC/          # Data acquisition, QC, scVI, DGE, trajectory
│   ├── phase2_GRN/         # pySCENIC, AUCell, scMultiomeGRN, ATAC validation
│   ├── phase2_LR/          # CellChat ligand-receptor (R)
│   ├── phase3_structure/   # Structure fetch, domain trimming, fpocket
│   ├── phase4_docking/     # Docking prep, consensus scoring, MM-GBSA
│   ├── phase5_validation/  # CellOracle perturbation
│   ├── phase6_figures/     # All figure generation scripts
│   ├── tables/             # Table generation scripts
│   └── slurm/              # SLURM job scripts for TRUBA HPC (steps 14–42)
├── envs/
│   ├── scanpy_env.yml      # Two-env spec: scenic (pySCENIC) + scanpy_env (analysis/figures)
│   ├── scmultiomegrn.yml   # scMultiomeGRN, CellOracle, velocyto, PyTorch
│   ├── cellchat_r.yml      # R 4.3.3, CellChat 2.2, Bioconductor deps
│   └── mmpbsa.yml          # OpenBabel, GROMACS tools, gmx_MMPBSA, meeko, pydeseq2
├── data/
│   ├── raw/                # Downloaded h5ad files (not in repo — see Data download)
│   ├── compounds/          # ChEMBL CNS compound library (CSV + PDBQT)
│   ├── structures/         # AF2/PDB protein structures (trimmed, prepared)
│   └── resources/          # HOCOMOCO v11 motif databases
└── results/
    ├── figures/            # All main and supplementary figures (PDF + PNG)
    ├── tables/             # Summary tables (CSV)
    ├── phase1/             # Trajectory, DGE, pseudobulk outputs
    ├── phase2/             # GRN adjacency, AUCell, ATAC, CellChat outputs
    ├── phase4/             # Docking scores, MM-GBSA summary
    └── phase5/             # CellOracle, GSE95587 replication outputs
```

---

## Reproducing figures and tables (no HPC required)

All main and supplementary figures, and all manuscript tables, can be reproduced locally from pre-computed result files already present in `results/`. No HPC, no GROMACS, no GNINA, no large data downloads required.

**Setup (one environment, ~5 minutes):**

```bash
conda env create -f envs/scanpy_env.yml -n scanpy_env
conda activate scanpy_env
export DAM_DRUG_DIR=/path/to/DAM-DRUG
```

**Run all figures:**

```bash
cd $DAM_DRUG_DIR
python code/phase6_figures/fig2_targets.py
python code/phase6_figures/fig3_docking.py
python code/phase6_figures/fig4_validation.py
python code/phase6_figures/fig_cellchat.py
python code/phase6_figures/supp_fig_s1_regulon_heatmap.py
python code/phase6_figures/supp_fig_s2_celloracle.py
python code/phase6_figures/supp_fig_s3_qc.py
python code/phase6_figures/supp_fig_s5_md_rmsd.py
python code/phase6_figures/supp_fig_s6_af2_quality.py
python code/phase6_figures/supp_fig_s7_bhlhe_coexpr.py
python code/phase6_figures/supp_fig_s8_slit2_robo2_expr.py
```

> `fig1_atlas.py` requires the trajectory h5ad file (~3.3 GB, from SEA-AD). Download it first with `bash code/phase1_QC/01_data_acquisition.sh` (requires AWS CLI).

**Run all tables:**

```bash
python code/tables/make_table2_markers.py
python code/tables/make_table3_target_scorecard.py
python code/tables/make_table4_5_docking.py
```

> `make_table1_cohort.py` requires the SEA-AD h5ad.

Outputs are written to `results/figures/` and `results/tables/` and can be compared directly to the committed PDFs and CSVs in the repository.

**What requires HPC:** The upstream computations (pySCENIC GRN inference, AutoDock Vina/GNINA docking, GROMACS MD, CellOracle, pseudobulk DESeq2) were run on the TRUBA ARF cluster (TÜBİTAK ULAKBİM). GNINA has no macOS binary. Full pipeline reproduction requires a Linux HPC with SLURM; see the [Execution order](#execution-order) section below for complete instructions.

---

## Prerequisites

### Hardware

| Step | Minimum | Used in this study |
|---|---|---|
| Phases 1, 3, 6, tables | 32 GB RAM, 8 cores | MacBook Pro (local) |
| pySCENIC (Phase 2) | 64 GB RAM, 16 cores | TRUBA ARF node, 128 GB RAM |
| scMultiomeGRN (Phase 2) | 64 GB RAM, GPU optional | TRUBA ARF node |
| Docking (Phase 4) | 32 GB RAM, 20 cores | TRUBA ARF node |
| MM-GBSA / MD (Phase 4) | 64 GB RAM, 20 cores | TRUBA ARF node, GROMACS 2024.1 (CPU) |
| CellOracle (Phase 5) | 32 GB RAM, 8 cores | TRUBA ARF node |

> Steps marked **HPC** in the execution table below require a SLURM cluster. All other steps run on a standard laptop/workstation.

### External binaries (not conda-managed)

Install these before running Phases 3–4:

| Tool | Version | Install |
|---|---|---|
| AWS CLI | any | `brew install awscli` or https://aws.amazon.com/cli/ |
| synapseclient | any | `pip install synapseclient` |
| AutoDock Vina | 1.2 | https://github.com/ccsb-scripps/AutoDock-Vina/releases — place binary at `~/apps/vina/vina` or set `VINA` env var |
| GNINA | 1.0 | https://github.com/gnina/gnina/releases — place binary at `~/apps/gnina/gnina` or set `GNINA` env var |
| OpenBabel | 3.1.1 | `conda install -c conda-forge openbabel=3.1.1` (in mmpbsa env) |
| fpocket | 4.x | https://github.com/Discngine/fpocket — build from source; place `fpocket` in `PATH` |
| GROMACS | 2024.1 | System-level install or HPC module (`module load apps/gromacs/2024.1-oneapi2024`) |

---

## Environment setup

All environments are specified in `envs/`. Set the project root environment variable first — every script reads it:

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG   # set this before running any script
```

### Environment 1 — scenic (pySCENIC GRN inference)

Used by: Phase 2 GRN steps 05, 25, 26, 28

```bash
conda env create -f envs/scenic.yml
conda activate scenic
```

### Environment 2 — scanpy_env (main analysis + figures)

Used by: Phases 1, 2 (aggregate), 3, 5 (partial), 6, tables

```bash
conda env create -f envs/scanpy_env.yml
conda activate scanpy_env
```

### Environment 3 — scMultiomeGRN (multi-omic GRN + CellOracle)

Used by: Phase 2 scMultiomeGRN (steps 30–33), Phase 2 ATAC (step 33), Phase 5 CellOracle (step 35)

```bash
conda env create -f envs/scmultiomegrn.yml
conda activate scMultiomeGRN
# Install build-time deps first, then velocyto, then celloracle:
pip install "cython" "numpy<2"
pip install velocyto==0.17.17 --no-build-isolation
pip install celloracle==0.20.0
# PyTorch with CUDA 12.1 (for TRUBA V100; omit --index-url for CPU-only):
pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 \
    --index-url https://download.pytorch.org/whl/cu121
pip install torch_geometric lightning
```

### Environment 4 — cellchat_r (CellChat in R)

Used by: Phase 2 LR steps 36–39

```bash
conda env create -f envs/cellchat_r.yml
conda activate cellchat_r
# Install CellChat from GitHub (exact version used: 2.2.0.9001):
Rscript -e 'devtools::install_github("jinworks/CellChat@v2.2.0")'
```

### Environment 5 — mmpbsa (docking + MM-GBSA + pydeseq2)

Used by: Phase 4 docking/MM-GBSA (steps 13–23), Phase 5 bulk DESeq2 (steps 19, 24)

```bash
conda env create -f envs/mmpbsa.yml
conda activate mmpbsa
export PATH="$HOME/.local/bin:$PATH"
# gmx_MMPBSA and meeko are pip-installed to ~/.local/bin:
pip install gmx_MMPBSA==1.6.4 meeko==0.7.1 pydeseq2==0.5.4 acpype==2023.10.27
```

---

## Data download

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG
bash code/phase1_QC/01_data_acquisition.sh
```

This script downloads:

| Dataset | File | Size | Required for |
|---|---|---|---|
| SEA-AD microglia (primary) | `SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad` | ~3.3 GB | All phases |
| SEA-AD MTG RNAseq | `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad` | ~36 GB | CellChat (Phase 2 LR) |
| Mathys et al. 2019 | via Synapse syn18485175 | ~1 GB | Supplementary only |

> The MTG ATACseq file (`SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad`, ~18 GB) is commented out in the script. Download it manually before running step 33 (ATAC DA validation) if needed.

**GSE95587** (bulk replication, Phase 5) is downloaded automatically by `code/slurm/19_deseq2_gse95587.slurm` via `GEOparse`.

**Synapse authentication** (Mathys download only):
```bash
pip install synapseclient
synapse login   # enter Synapse credentials; one-time
```

**AWS CLI** is required for SEA-AD downloads (public bucket, no credentials needed):
```bash
# macOS
brew install awscli
# Linux
pip install awscli
```

---

## Execution order

Scripts are globally numbered 01–42. Numbers in `code/slurm/` correspond to the same step run on HPC.

**Legend:** `[local]` = runs on laptop/workstation | `[HPC]` = requires SLURM cluster | `[optional]` = not required for main figures

### Phase 1 — QC and cell state annotation

| Step | Script | Environment | Notes |
|---|---|---|---|
| 01 | `code/phase1_QC/01_data_acquisition.sh` | system | Downloads SEA-AD; see Data download above |
| 02a | `code/phase1_QC/02a_explore_microglia_object.py` | scanpy_env | [local] Explore pre-release microglia h5ad (already QC'd); marker and batch checks |
| 02b | `code/phase1_QC/02b_QC_microglia_isolation.py` | scanpy_env | [local] QC and isolation from full MTG RNAseq + Mathys 2019; outputs microglia h5ad |
| 03 | `code/phase1_QC/03_batch_correction_scVI.py` | scanpy_env | [local] scVI latent space; UMAP |
| 04 | `code/phase1_QC/04_DGE_microglial_states.py` | scanpy_env | [local] Wilcoxon rank-sum per substate |
| 24 | `code/slurm/24_pseudobulk_deseq2.slurm` | mmpbsa | [HPC] Pseudobulk DESeq2 validation |
| 27 | `code/slurm/27_trajectory_paga.slurm` | scanpy_env | [HPC] PAGA + diffusion pseudotime (run after step 26) |

```bash
# Steps 02–04: local
conda activate scanpy_env
export DAM_DRUG_DIR=/path/to/DAM-DRUG
python code/phase1_QC/02a_explore_microglia_object.py
python code/phase1_QC/02b_QC_microglia_isolation.py
python code/phase1_QC/03_batch_correction_scVI.py
python code/phase1_QC/04_DGE_microglial_states.py

# Step 24 (pseudobulk): HPC
sbatch code/slurm/24_pseudobulk_deseq2.slurm

# Step 27 (trajectory): HPC, run after AUCell (step 26)
sbatch code/slurm/27_trajectory_paga.slurm
```

### Phase 2 — Gene regulatory network inference

#### pySCENIC (multi-seed GRNBoost2)

| Step | Script | Environment | Notes |
|---|---|---|---|
| 05 | `code/phase2_GRN/05_pySCENIC_GRN.py` | scenic | [local] Prepare pySCENIC inputs |
| 25 | `code/slurm/25_grn_multiseed.slurm` | scenic | [HPC] 5-seed GRNBoost2 + RcisTarget (array job) |
| 26 | `code/slurm/26_run_ctx_aucell.slurm` | scenic | [HPC] AUCell scoring across all cells |
| 25b | `code/phase2_GRN/25_aggregate_grn.py` | scanpy_env | [local] Aggregate multi-seed GRN; filter regulons |
| 28 | `code/slurm/28_regulon_pseudotime_corr.slurm` | scenic | [HPC] Regulon–pseudotime Spearman correlation |
| 28b | `code/phase2_GRN/28_regulon_pseudotime_correlation.py` | scanpy_env | [local] Post-process correlation results |

```bash
# Step 05: local — prepare inputs
conda activate scenic
python code/phase2_GRN/05_pySCENIC_GRN.py

# Steps 25, 26, 28: HPC
sbatch code/slurm/25_grn_multiseed.slurm
# wait for completion, then:
sbatch code/slurm/26_run_ctx_aucell.slurm
sbatch code/slurm/28_regulon_pseudotime_corr.slurm

# Aggregate results: local
conda activate scanpy_env
python code/phase2_GRN/25_aggregate_grn.py
python code/phase2_GRN/28_regulon_pseudotime_correlation.py
```

#### scMultiomeGRN (multi-omic cross-validation)

| Step | Script | Environment | Notes |
|---|---|---|---|
| 30 | `code/slurm/30_setup_scMultiomeGRN.slurm` | scMultiomeGRN | [HPC] Install scMultiomeGRN package |
| 31 | `code/phase2_GRN/31_prep_scMultiomeGRN_input.py` | scMultiomeGRN | [local] Prepare MTX input files |
| 31b | `code/slurm/31_scMultiomeGRN_env.slurm` | scMultiomeGRN | [HPC] Validate env + run inference |
| 32 | `code/slurm/32_scMultiomeGRN_infer.slurm` | scMultiomeGRN | [HPC] scMultiomeGRN inference |

#### ATAC chromatin validation

| Step | Script | Environment | Notes |
|---|---|---|---|
| 32b | `code/phase2_GRN/32_atac_da_validation.py` | scMultiomeGRN | [local] DA peak analysis |
| 33 | `code/slurm/33_atac_da_validation.slurm` | scMultiomeGRN | [HPC] Full ATAC DA on cluster |

> Download `SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad` (~18 GB) before running step 33.

### Phase 2 — Cell–cell communication (CellChat)

| Step | Script | Environment | Notes |
|---|---|---|---|
| 36 | `code/phase2_LR/36_prep_mtg_for_cellchat.py` | scanpy_env | [local] Extract MTG microglia subset |
| 36b | `code/slurm/36_install_r_cellchat.slurm` | cellchat_r | [HPC] Install R packages (one-time) |
| 37 | `code/slurm/37_prep_mtg.slurm` | cellchat_r | [HPC] Prepare MTG Seurat object |
| 38 | `code/slurm/38_cellchat.slurm` | cellchat_r | [HPC] CellChat inference (all cell types) |
| 39 | `code/slurm/39_cellchat_final.slurm` | cellchat_r | [HPC] Finalize and export CellChat object |

### Phase 3 — Protein structure preparation

All steps run locally with the `mmpbsa` or `scanpy_env` environment.

| Step | Script | Environment | Notes |
|---|---|---|---|
| 06 | `code/phase3_structure/06_fetch_structures.py` | scanpy_env | Download AF2 + PDB structures for 6 TF targets |
| 07 | `code/phase3_structure/07_trim_domains.py` | scanpy_env | Trim to DNA-binding domain only |
| 08 | `code/phase3_structure/08_prepare_structures.py` | mmpbsa | PDBFixer → PDBQT (obabel, Gasteiger charges) |
| 09 | `code/phase3_structure/09_run_fpocket.sh` | system | fpocket pocket detection (requires `fpocket` in PATH) |
| 10 | `code/phase3_structure/10_parse_fpocket.py` | scanpy_env | Parse fpocket scores; select docking pockets |
| 11 | `code/phase3_structure/11_druggability_summary.py` | scanpy_env | Druggability summary table |

```bash
conda activate scanpy_env
python code/phase3_structure/06_fetch_structures.py
python code/phase3_structure/07_trim_domains.py

conda activate mmpbsa
export PATH="$HOME/.local/bin:$PATH"
python code/phase3_structure/08_prepare_structures.py

bash code/phase3_structure/09_run_fpocket.sh

conda activate scanpy_env
python code/phase3_structure/10_parse_fpocket.py
python code/phase3_structure/11_druggability_summary.py
```

### Phase 4 — Virtual screening

ChEMBL compound library (1,677 FDA-approved CNS drugs) is in `data/compounds/`. PDBQT preparation and docking require HPC for throughput.

| Step | Script | Environment | Notes |
|---|---|---|---|
| 12 | `code/phase4_docking/12_docking_prep.py` | mmpbsa | [local] Compute Vina grid boxes; receptor PDBQT |
| 13 | `code/phase4_docking/13_prep_ligands.py` | mmpbsa | [local/HPC] SMILES → 3D SDF (obabel) → PDBQT (meeko) |
| 14 | `code/slurm/14_run_vina.slurm` | mmpbsa | [HPC] AutoDock Vina screen (20 parallel workers) |
| 15 | `code/slurm/15_run_gnina.slurm` | mmpbsa | [HPC] GNINA CNN rescoring of Vina poses |
| 16 | `code/phase4_docking/16_consensus_filter.py` | scanpy_env | [local] Consensus filter (Vina ≤ −8.5 AND CNN ≥ 0.75) |
| 17 | `code/phase4_docking/17_prep_mmpbsa.py` | mmpbsa | [local] Prepare MM-GBSA workdirs (GAFF2 topology) |
| 17b | `code/slurm/17_run_mmpbsa.slurm` | mmpbsa | [HPC] gmx_MMPBSA 1 ns NPT MM-GBSA (array job) |
| 17c | `code/slurm/17_collect_mmpbsa.py` | mmpbsa | [local] Aggregate FINAL_RESULTS_MMPBSA.dat → CSV |
| 18 | `code/slurm/18_run_md100ns.slurm` | mmpbsa | [HPC] 100 ns GROMACS 2024.1 MD stability (CPU) |
| 18b | `code/slurm/18b_core_rmsd_analysis.slurm` | mmpbsa | [HPC] Core backbone + ligand RMSD analysis |

**Tier-2 library (approved drugs beyond CNS filter) — optional:**

| Step | Script | Notes |
|---|---|---|
| 20 | `code/slurm/20_fetch_tier2.py` + `code/slurm/20_prep_tier2_library.slurm` | Fetch ChEMBL approved drugs; prep PDBQT |
| 21 | `code/slurm/21_run_docking_tier2.slurm` | Vina screen of tier-2 |
| 22 | `code/slurm/22_run_mmpbsa_tier2.slurm` | MM-GBSA on tier-2 hits |
| 23 | `code/slurm/23_run_selectivity_docking.slurm` | Selectivity docking (off-target receptors) |

```bash
# Local prep steps
conda activate mmpbsa && export PATH="$HOME/.local/bin:$PATH"
export DAM_DRUG_DIR=/path/to/DAM-DRUG
python code/phase4_docking/12_docking_prep.py
python code/phase4_docking/13_prep_ligands.py

# HPC docking jobs
sbatch code/slurm/14_run_vina.slurm
# wait for completion, then:
sbatch code/slurm/15_run_gnina.slurm

# Local consensus filter
conda activate scanpy_env
python code/phase4_docking/16_consensus_filter.py

# MM-GBSA prep (local) + scoring (HPC)
conda activate mmpbsa && export PATH="$HOME/.local/bin:$PATH"
python code/phase4_docking/17_prep_mmpbsa.py
sbatch code/slurm/17_run_mmpbsa.slurm
# wait, then collect:
python code/slurm/17_collect_mmpbsa.py

# 100 ns MD: HPC only
sbatch code/slurm/18_run_md100ns.slurm
sbatch code/slurm/18b_core_rmsd_analysis.slurm
```

### Phase 5 — Validation

| Step | Script | Environment | Notes |
|---|---|---|---|
| 19 | `code/slurm/19_deseq2_gse95587.slurm` | mmpbsa | [HPC] Bulk replication: GSE95587 DESeq2 |
| 35 | `code/phase5_validation/35_celloracle_perturbation.py` | scMultiomeGRN | [local/HPC] CellOracle in-silico TF KO |
| 35b | `code/slurm/35_celloracle_perturbation.slurm` | scMultiomeGRN | [HPC] CellOracle on cluster |

```bash
# Bulk replication: HPC
sbatch code/slurm/19_deseq2_gse95587.slurm

# CellOracle: local (requires >32 GB RAM) or HPC
conda activate scMultiomeGRN
export DAM_DRUG_DIR=/path/to/DAM-DRUG
python code/phase5_validation/35_celloracle_perturbation.py
# or: sbatch code/slurm/35_celloracle_perturbation.slurm
```

### Phase 6 — Figures

All figure scripts run locally with `scanpy_env`. Inputs are pre-computed CSVs in `results/`; large h5ad files are not required for most figures.

```bash
conda activate scanpy_env
export DAM_DRUG_DIR=/path/to/DAM-DRUG

python code/phase6_figures/fig1_atlas.py          # Fig 1: microglial atlas (requires trajectory h5ad)
python code/phase6_figures/fig2_targets.py        # Fig 2: target prioritisation
python code/phase6_figures/fig3_docking.py        # Fig 3: docking results
python code/phase6_figures/fig4_validation.py     # Fig 4: multi-modal validation
python code/phase6_figures/fig_cellchat.py        # CellChat panel
python code/phase6_figures/graphical_abstract.py  # Graphical abstract

# Supplementary figures
python code/phase6_figures/supp_fig_s1_regulon_heatmap.py
python code/phase6_figures/supp_fig_s2_celloracle.py
python code/phase6_figures/supp_fig_s3_qc.py
python code/phase6_figures/supp_fig_s4_docking_roc.py
python code/phase6_figures/supp_fig_s5_md_rmsd.py
python code/phase6_figures/supp_fig_s6_af2_quality.py
python code/phase6_figures/supp_fig_s7_bhlhe_coexpr.py
python code/phase6_figures/supp_fig_s8_slit2_robo2_expr.py
```

> **Fig 2 note:** `regulon_auc_by_state_aggregated.csv` and `regulon_pseudotime_corr.csv` must be synced from TRUBA after completing steps 25–28 before running `fig2_targets.py` on a local machine.

### Tables

```bash
conda activate scanpy_env
export DAM_DRUG_DIR=/path/to/DAM-DRUG

python code/tables/make_table1_cohort.py          # Cohort demographics
python code/tables/make_table2_markers.py         # Substate marker genes
python code/tables/make_table3_target_scorecard.py  # TF target scorecard
python code/tables/make_table4_5_docking.py       # Docking + MM-GBSA results
```

---

## Software versions

### Python environments

| Package | Version | Environment |
|---|---|---|
| Python | 3.10 | all |
| scanpy | 1.11.5 | scanpy_env |
| anndata | 0.11.4 | scanpy_env |
| numpy | 2.2.6 | scanpy_env |
| pandas | 2.3.3 | scanpy_env |
| scipy | 1.15.2 | scanpy_env |
| matplotlib | 3.10.8 | scanpy_env |
| seaborn | 0.13.2 | scanpy_env |
| umap-learn | 0.5.11 | scanpy_env |
| leidenalg | 0.11.0 | scanpy_env |
| statsmodels | 0.14.6 | scanpy_env |
| pyscenic | 0.12.1 | scenic |
| arboreto | 0.1.6 | scenic |
| ctxcore | 0.2.0 | scenic |
| loompy | 3.0.8 | scenic |
| dask | 2023.12.0 | scenic |
| numpy | 1.26.4 | scenic |
| scvi-tools | 1.4.2 | scenic (base) |
| scrublet | 0.2.3 | scenic (base) |
| celloracle | 0.20.0 | scMultiomeGRN |
| velocyto | 0.17.17 | scMultiomeGRN |
| anndata | 0.10.8 | scMultiomeGRN |
| torch | 2.4.0 | scMultiomeGRN |
| pydeseq2 | 0.5.4 | mmpbsa |
| openbabel | 3.1.1 | mmpbsa |
| meeko | 0.7.1 | mmpbsa (pip, ~/.local) |
| gmx_MMPBSA | 1.6.4 | mmpbsa (pip, ~/.local) |
| acpype | 2023.10.27 | mmpbsa |
| Python | 3.12 | mmpbsa |

### R environment

| Package | Version |
|---|---|
| R | 4.3.3 |
| CellChat | 2.2.0.9001 |
| NMF | 0.28 |
| ComplexHeatmap | 2.18.0 |
| limma | 3.58.1 |
| MAST | 1.28.0 |
| igraph | 2.1.4 |
| ggplot2 | 3.5.2 |

### External tools

| Tool | Version |
|---|---|
| AutoDock Vina | 1.2 |
| GNINA | 1.0 |
| fpocket | 4.x |
| GROMACS | 2024.1 (gmx_mpi; TRUBA system module) |
| MEME Suite | 5.4.1 |

---

## Reproducibility notes

- All stochastic steps use fixed seeds. See manuscript §Random seeds for full table.
- **GROMACS MD** (`gen_seed = -1`): velocity generation uses a system-generated random seed. RMSD values in the manuscript are not bit-reproducible without the archived trajectory files (stored on TRUBA ARF scratch). The trajectory files can be made available on request.
- **pySCENIC** uses 5 seeds (42, 1, 2, 3, 4) with a minimum consensus of 3/5. Results are consensus-stable across seeds.
- All scripts read the project root from `$DAM_DRUG_DIR` (with `os.environ.get("DAM_DRUG_DIR", "<original_path>")` fallback). Set this variable before running any script.
- Environment specs with exact version pins are in `envs/`. Versions were recovered from container images (`scanpy.sif`, `cellchat.sif`) and miniconda backup archives (`~/miniconda_backup.tar`) used during the original analysis on TRUBA.

---

## Data availability

| Dataset | Source | License |
|---|---|---|
| SEA-AD microglial atlas | [Allen Brain Cell Atlas](https://portal.brain-map.org/explore/Seattle-Alzheimers-Disease) — S3: `s3://sea-ad-single-cell-profiling/` | CC BY 4.0 |
| SEA-AD MTG (RNAseq + ATACseq) | Allen Brain Cell Atlas (same S3 bucket) | CC BY 4.0 |
| Mathys et al. 2019 | [Synapse syn18485175](https://www.synapse.org/#!Synapse:syn18485175) | CC BY 4.0 |
| GSE95587 (bulk replication) | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95587) | Public |
| ChEMBL compound library | [ChEMBL](https://www.ebi.ac.uk/chembl/) (max_phase = 4, CNS filter) | CC BY-SA 3.0 |
| AlphaFold2 structures | [EBI AlphaFold Database v4](https://alphafold.ebi.ac.uk/) | CC BY 4.0 |
| PDB structures | [RCSB PDB](https://www.rcsb.org/) | Open |
| HOCOMOCO v11 motifs | [HOCOMOCO](https://hocomoco11.autosome.org/) | CC BY 4.0 |

---

## Citation

```
Özkurt Ç, Kelicen-Uğur P. Integrative Single-Cell Analysis of Alzheimer's Disease Microglia
Identifies IKZF1 as a Late-Stage Neuroinflammatory Regulator and Prioritises Tafamidis and
Diflunisal as Candidate Repurposable Therapeutics. (manuscript in preparation)
```

---

## License

Code is released under the MIT License. See `LICENSE`.
