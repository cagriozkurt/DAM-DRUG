# DAM-DRUG

**Integrative Single-Cell Analysis of Alzheimer's Disease Microglia Identifies IKZF1 as a Late-Stage Neuroinflammatory Regulator and Generates Low-Confidence Computational Hypotheses for Tafamidis and Diflunisal as Repurposing Candidates**

Çağrı Özkurt — Hacettepe University

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19535749.svg)](https://doi.org/10.5281/zenodo.19535749)

> Özkurt Ç. (manuscript in preparation).

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
│   ├── phase1_QC/          # Data acquisition, QC, DGE, trajectory
│   ├── phase2_GRN/         # pySCENIC, AUCell, scMultiomeGRN, ATAC validation
│   ├── phase2_LR/          # CellChat ligand-receptor (R)
│   ├── phase3_structure/   # Structure fetch, domain trimming, fpocket
│   ├── phase4_docking/     # Docking prep, consensus scoring, MM-GBSA
│   ├── phase5_validation/  # CellOracle perturbation
│   ├── phase6_figures/     # All figure generation scripts
│   ├── tables/             # Table generation scripts
│   └── slurm/              # SLURM job scripts for TRUBA HPC (steps 14–42)
├── docker/                 # Dockerfiles for reproducible environments
├── envs/                   # Conda environment YAMLs
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

## Part 1 — Reproduce figures and tables (no HPC required)

All main and supplementary figures and manuscript tables can be reproduced from pre-computed result files already present in `results/`. No HPC, no GROMACS, no GNINA, and no large data downloads are required.

### Step 1 — Set up environment

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG
apptainer pull "$DAM_DRUG_DIR/containers/scenic.sif" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest
alias pyrun="apptainer exec --bind \"\$DAM_DRUG_DIR:\$DAM_DRUG_DIR\" \"\$DAM_DRUG_DIR/containers/scenic.sif\" conda run -n scanpy_env python"
alias scenicrun="apptainer exec --bind \"\$DAM_DRUG_DIR:\$DAM_DRUG_DIR\" --env LD_LIBRARY_PATH=/opt/conda/envs/scenic/lib \"\$DAM_DRUG_DIR/containers/scenic.sif\" conda run -n scenic python"
```

> `pyrun` (`scanpy_env`) and `scenicrun` (`scenic`) are the aliases used for all commands below.

### Step 2 — Run figures

```bash
cd $DAM_DRUG_DIR

pyrun code/phase6_figures/02_fig2_targets.py
pyrun code/phase6_figures/03_fig3_docking.py
pyrun code/phase6_figures/04_fig4_validation.py
pyrun code/phase6_figures/05_fig_cellchat.py
pyrun code/phase6_figures/07_supp_fig_s1_regulon_heatmap.py
pyrun code/phase6_figures/08_supp_fig_s2_celloracle.py
pyrun code/phase6_figures/10_supp_fig_s5_md_rmsd.py
pyrun code/phase6_figures/11_supp_fig_s6_af2_quality.py
pyrun code/phase6_figures/12_supp_fig_s7_bhlhe_coexpr.py
```

Outputs are written to `results/figures/` as PDF and PNG.

**Four figures require additional data not included in the repo:**

| Figure | Missing file | How to obtain |
|--------|-------------|---------------|
| `01_fig1_atlas.py`, `12_supp_fig_s7_bhlhe_coexpr.py` | SEA-AD h5ad (~3.3–36 GB) | `bash code/phase1_QC/01_data_acquisition.sh` (requires AWS CLI) |
| `01_fig1_atlas.py`, `09_supp_fig_s3_qc.py` | `results/phase1/trajectory/microglia_trajectory.h5ad` | Run `sbatch code/slurm/05_trajectory_paga.slurm` on HPC (128 GB RAM, ~8 h) after downloading the SEA-AD microglia h5ad |
| `13_supp_fig_s8_slit2_robo2_expr.py` | `results/phase2/LR/prep/counts_raw.h5` (873 MB) | `wget -P results/phase2/LR/prep/ https://zenodo.org/records/19535759/files/counts_raw.h5` |

### Step 3 — Run tables

```bash
pyrun code/tables/make_table2_markers.py
pyrun code/tables/make_table3_target_scorecard.py
pyrun code/tables/make_table4_5_docking.py
```

> `make_table1_cohort.py` requires the SEA-AD microglia h5ad.

---

## Part 2 — Full pipeline reproduction (HPC required)

The upstream computations (pySCENIC GRN inference, AutoDock Vina/GNINA docking, GROMACS MD, CellOracle, pseudobulk DESeq2) were run on the TRUBA ARF cluster (TÜBİTAK ULAKBİM). Full pipeline reproduction requires a Linux HPC with SLURM, GROMACS 2024.1, and GPU access for GNINA.

### Prerequisites

#### Hardware

| Phase | Minimum | Used in this study |
|-------|---------|-------------------|
| Phases 1, 3, 6, tables | 32 GB RAM, 8 cores | MacBook Pro (local) |
| pySCENIC (Phase 2) | 64 GB RAM, 16 cores | TRUBA ARF node, 128 GB RAM |
| scMultiomeGRN (Phase 2) | 64 GB RAM, GPU optional | TRUBA ARF node |
| Docking (Phase 4) | 32 GB RAM, 20 cores | TRUBA ARF node |
| MM-GBSA / MD (Phase 4) | 64 GB RAM, 20 cores | TRUBA ARF node, GROMACS 2024.1 (CPU) |
| CellOracle (Phase 5) | 32 GB RAM, 8 cores | TRUBA ARF node |

#### External binaries

| Tool | Version | Install |
|------|---------|---------|
| AWS CLI | any | `brew install awscli` or https://aws.amazon.com/cli/ |
| AutoDock Vina | 1.2 | https://github.com/ccsb-scripps/AutoDock-Vina/releases — place at `$DAM_DRUG_DIR/bin/vina` or set `VINA` env var |
| GNINA | 1.3.2 | https://github.com/gnina/gnina/releases — place at `$DAM_DRUG_DIR/bin/gnina` or set `GNINA` env var; requires CUDA/cuDNN 9 |
| OpenBabel | 3.1.1 | `conda install -c conda-forge openbabel=3.1.1` (in mmpbsa env) |
| fpocket | 4.x | Included in `fpocket-env.sif` (compiled from source in `docker/Dockerfile.fpocket`) |
| GROMACS | 2024.1 | HPC module: `module load apps/gromacs/2024.1-oneapi2024` |

#### Apptainer images

Most SLURM scripts use Apptainer images pulled automatically on first run. Pull them manually if preferred:

```bash
apptainer pull "$DAM_DRUG_DIR/containers/scenic.sif"                docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest
apptainer pull "$DAM_DRUG_DIR/containers/dam-drug-r.sif"            docker://ghcr.io/cagriozkurt/dam-drug-r:latest
apptainer pull "$DAM_DRUG_DIR/containers/dam-drug-scmultiomegrn.sif" docker://ghcr.io/cagriozkurt/dam-drug-scmultiomegrn:latest
apptainer pull "$DAM_DRUG_DIR/containers/fpocket-env.sif"           docker://ghcr.io/cagriozkurt/dam-drug-fpocket:latest
```

| Image | Covers |
|-------|--------|
| `scenic.sif` | analysis + figures (scanpy, pySCENIC, pydeseq2) |
| `dam-drug-r.sif` | CellChat/R |
| `dam-drug-scmultiomegrn.sif` | CellOracle + PyTorch GRN (phases 2 & 5) |
| `fpocket-env.sif` | PDBFixer + fpocket pocket detection (phase 3) |

> GPU scripts pass `--nv` to expose the host GPU driver. SLURM scripts handle this automatically.

#### Conda environment (pipeline-only, not containerized)

`mmpbsa` requires the GROMACS system module and cannot be containerized portably:

```bash
conda env create -f envs/mmpbsa.yml
conda activate mmpbsa
pip install gmx_MMPBSA==1.6.4 meeko==0.7.1 acpype==2023.10.27
```

### Data download

```bash
bash code/phase1_QC/01_data_acquisition.sh
```

| Dataset | File | Size | Required for |
|---------|------|------|-------------|
| SEA-AD microglia (primary) | `SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad` | ~3.3 GB | All phases |
| SEA-AD MTG RNAseq | `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad` | ~36 GB | CellChat (Phase 2 LR) |

> The MTG ATACseq file (`SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad`, ~18 GB) is commented out in the script. Download it manually before running step 33 (ATAC DA validation).

**GSE95587** (bulk replication, Phase 5) must be downloaded manually before running `32_deseq2_gse95587.slurm`:

```bash
mkdir -p $DAM_DRUG_DIR/data/raw/Mayo/suppl
cd $DAM_DRUG_DIR/data/raw/Mayo/suppl

# Series matrix (~2 MB)
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95587/suppl/GSE95587_series_matrix.txt.gz"

# Per-sample count files (~500 MB total) — one TSV.gz per GSM accession
wget -r -nd -np -A "GSM*.tsv.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95587/suppl/"
```

**pySCENIC resource files** must be present before running `06_run_pySCENIC.slurm`. The SLURM script downloads them automatically if missing, but on clusters without compute-node internet access download them on the login node first:

```bash
mkdir -p $DAM_DRUG_DIR/data/resources/motifs

# Rankings database (~1.5 GB)
wget -P $DAM_DRUG_DIR/data/resources/motifs/ \
  https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# Motif-to-TF table (~50 MB)
wget -P $DAM_DRUG_DIR/data/resources/motifs/ \
  https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

### Execution order

Scripts are numbered 01–42. Numbers in `code/slurm/` correspond to the same step on HPC.

**Legend:** `[HPC]` = requires SLURM cluster | `[optional]` = not required for main figures

> All snippets use Apptainer (`pyrun` / `scenicrun`). `conda activate mmpbsa` is the only exception — GROMACS cannot be containerized portably.

#### Phase 1 — QC and cell state annotation

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 01 | `code/phase1_QC/01_data_acquisition.sh` | system | Downloads SEA-AD datasets |
| 02 | `code/slurm/02_run_explore.slurm` | scanpy_env | [HPC] Explore pre-release microglia h5ad; marker and batch checks |
| 03 | `code/slurm/03_run_DGE.slurm` | scanpy_env | [HPC] Wilcoxon rank-sum DGE per substate |
| 04 | `code/slurm/04_pseudobulk_deseq2.slurm` | mmpbsa | [HPC] Pseudobulk DESeq2 validation |
| 05 | `code/slurm/05_trajectory_paga.slurm` | scanpy_env | [HPC] PAGA + diffusion pseudotime (128 GB, ~8 h) |

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG
sbatch code/slurm/02_run_explore.slurm
sbatch code/slurm/03_run_DGE.slurm
sbatch code/slurm/04_pseudobulk_deseq2.slurm
sbatch code/slurm/05_trajectory_paga.slurm   # after step 26
```

#### Phase 2 — Gene regulatory network inference

**pySCENIC (multi-seed GRNBoost2):**

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 06 | `code/slurm/06_run_pySCENIC.slurm` | scenic | [HPC] Download resources + prepare pySCENIC inputs (`01_pySCENIC_GRN.py`) |
| 07 | `code/slurm/07_grn_multiseed.slurm` | scenic | [HPC] 5-seed GRNBoost2 + RcisTarget (array job) |
| 08 | `code/slurm/08_run_ctx_aucell.slurm` | scenic | [HPC] AUCell scoring across all cells |
| 02 | `code/phase2_GRN/02_aggregate_grn.py` | scanpy_env | Aggregate multi-seed GRN; filter regulons |
| 09 | `code/slurm/09_regulon_pseudotime_corr.slurm` | scanpy_env | [HPC] Regulon–pseudotime Spearman correlation |
| 03 | `code/phase2_GRN/03_regulon_pseudotime_correlation.py` | scanpy_env | Post-process correlation results |

> Do not run `01_pySCENIC_GRN.py` directly — resource files in `data/resources/motifs/` must exist first (see download commands above). `06_run_pySCENIC.slurm` will also auto-download them if the compute node has internet access.

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG
sbatch code/slurm/06_run_pySCENIC.slurm

sbatch code/slurm/07_grn_multiseed.slurm
sbatch code/slurm/08_run_ctx_aucell.slurm   # after step 25

pyrun code/phase2_GRN/02_aggregate_grn.py

sbatch code/slurm/09_regulon_pseudotime_corr.slurm
```

**scMultiomeGRN (multi-omic cross-validation):**

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 11 | `code/slurm/11_setup_scMultiomeGRN.slurm` | scMultiomeGRN | [HPC] Download scMultiomeGRN tool source from Zenodo |
| 04 | `code/phase2_GRN/04_prep_scMultiomeGRN_input.py` | scMultiomeGRN | Prepare MTX input files |
| 12 | `code/slurm/12_scMultiomeGRN_infer.slurm` | scMultiomeGRN | [HPC] scMultiomeGRN inference |

```bash
sbatch code/slurm/11_setup_scMultiomeGRN.slurm

apptainer exec --bind "$DAM_DRUG_DIR:$DAM_DRUG_DIR" \
    "$DAM_DRUG_DIR/containers/dam-drug-scmultiomegrn.sif" \
    conda run -n scMultiomeGRN python \
    "$DAM_DRUG_DIR/code/phase2_GRN/04_prep_scMultiomeGRN_input.py"

sbatch code/slurm/12_scMultiomeGRN_infer.slurm
```

**ATAC chromatin validation:**

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 05 | `code/phase2_GRN/05_atac_da_validation.py` | scMultiomeGRN | DA peak analysis |
| 13 | `code/slurm/13_atac_da_validation.slurm` | scMultiomeGRN | [HPC] Full ATAC DA on cluster |

> Download `SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad` (~18 GB) before running step 5/13.

```bash
# local (small subset):
apptainer exec --bind "$DAM_DRUG_DIR:$DAM_DRUG_DIR" \
    "$DAM_DRUG_DIR/containers/dam-drug-scmultiomegrn.sif" \
    conda run -n scMultiomeGRN python \
    "$DAM_DRUG_DIR/code/phase2_GRN/05_atac_da_validation.py"

# or full run on HPC:
sbatch code/slurm/13_atac_da_validation.slurm
```

#### Phase 2 — Cell–cell communication (CellChat)

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 01 | `code/phase2_LR/01_prep_mtg_for_cellchat.py` | scanpy_env | Extract MTG microglia subset |
| 15 | `code/slurm/15_prep_mtg.slurm` | cellchat_r | [HPC] Prepare MTG Seurat object |
| 16 | `code/slurm/16_cellchat.slurm` | cellchat_r | [HPC] CellChat inference (all cell types) |
| 17 | `code/slurm/17_cellchat_final.slurm` | cellchat_r | [HPC] Finalize and export CellChat object |

> CellChat is pre-installed in `ghcr.io/cagriozkurt/dam-drug-r:latest`. SLURM scripts use this image automatically.

```bash
pyrun code/phase2_LR/01_prep_mtg_for_cellchat.py

sbatch code/slurm/15_prep_mtg.slurm
sbatch code/slurm/16_cellchat.slurm
sbatch code/slurm/17_cellchat_final.slurm
```

#### Phase 3 — Protein structure preparation

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 01 | `code/phase3_structure/01_fetch_structures.py` | scanpy_env | Download AF2 + PDB structures for 6 TF targets (login node) |
| 02 | `code/phase3_structure/02_trim_domains.py` | scanpy_env | Trim to DNA-binding domain (login node) |
| 18 | `code/slurm/18_run_phase3.slurm` | fpocket-env.sif | [HPC] PDBFixer → fpocket → parse → druggability summary (steps 3–6) |

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG
# Steps 1-2: run on login node (compute nodes have no internet)
pyrun code/phase3_structure/01_fetch_structures.py
pyrun code/phase3_structure/02_trim_domains.py

# Steps 3-6: prepare + fpocket + parse + summarize
sbatch code/slurm/18_run_phase3.slurm
```

#### Phase 4 — Virtual screening

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 01 | `code/phase4_docking/01_docking_prep.py` | mmpbsa | Compute Vina grid boxes; receptor PDBQT |
| 03 | `code/phase4_docking/03_prep_ligands.py` | mmpbsa | SMILES → 3D SDF → PDBQT |
| 20 | `code/slurm/20_run_vina.slurm` | mmpbsa | [HPC] AutoDock Vina screen (20 parallel workers) |
| 21 | `code/slurm/21_run_gnina.slurm` | mmpbsa | [HPC] GNINA CNN rescoring of Vina poses |
| 06 | `code/phase4_docking/06_consensus_filter.py` | scanpy_env | Consensus filter (Vina ≤ −8.5 AND CNN ≥ 0.75) |
| 07 | `code/phase4_docking/07_prep_mmpbsa.py` | mmpbsa | Prepare MM-GBSA workdirs (GAFF2 topology) |
| 22 | `code/slurm/22_run_mmpbsa.slurm` | mmpbsa | [HPC] gmx_MMPBSA 1 ns NPT MM-GBSA (array job) |
| 23 | `code/slurm/23_collect_mmpbsa.py` | mmpbsa | Aggregate FINAL_RESULTS_MMPBSA.dat → CSV |
| 24 | `code/slurm/24_run_md100ns.slurm` | mmpbsa | [HPC] 100 ns GROMACS 2024.1 MD stability |
| 25 | `code/slurm/25_core_rmsd_analysis.slurm` | mmpbsa | [HPC] Core backbone + ligand RMSD analysis |

**Tier-2 library (approved drugs beyond CNS filter) — optional:**

| Step | Script | Notes |
|------|--------|-------|
| 26 | `code/slurm/26_fetch_tier2.py` + `27_prep_tier2_library.slurm` | Fetch ChEMBL approved drugs; prep PDBQT |
| 28 | `code/slurm/28_run_docking_tier2.slurm` | Vina screen of tier-2 |
| 30 | `code/slurm/30_run_mmpbsa_tier2.slurm` | MM-GBSA on tier-2 hits |
| 31 | `code/slurm/31_run_selectivity_docking.slurm` | Selectivity docking (off-target receptors) |

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG
conda activate mmpbsa
python code/phase4_docking/01_docking_prep.py
python code/phase4_docking/03_prep_ligands.py

sbatch code/slurm/20_run_vina.slurm
sbatch code/slurm/21_run_gnina.slurm   # after step 14

pyrun code/phase4_docking/06_consensus_filter.py

conda activate mmpbsa
python code/phase4_docking/07_prep_mmpbsa.py
sbatch code/slurm/22_run_mmpbsa.slurm
python code/slurm/23_collect_mmpbsa.py

sbatch code/slurm/24_run_md100ns.slurm
sbatch code/slurm/25_core_rmsd_analysis.slurm
```

#### Phase 5 — Validation

| Step | Script | Environment | Notes |
|------|--------|-------------|-------|
| 32 | `code/slurm/32_deseq2_gse95587.slurm` | scenic.sif | [HPC] Bulk replication: GSE95587 DESeq2 |
| 01 | `code/phase5_validation/01_celloracle_perturbation.py` | scMultiomeGRN | CellOracle in-silico TF KO |

```bash
export DAM_DRUG_DIR=/path/to/DAM-DRUG
sbatch code/slurm/32_deseq2_gse95587.slurm

apptainer exec \
    --bind "$DAM_DRUG_DIR:$DAM_DRUG_DIR" \
    "$DAM_DRUG_DIR/containers/dam-drug-scmultiomegrn.sif" \
    conda run -n scMultiomeGRN python "$DAM_DRUG_DIR/code/phase5_validation/01_celloracle_perturbation.py"
# or on HPC: sbatch code/slurm/33_celloracle_perturbation.slurm
```

#### Phase 6 — Figures and tables

See [Part 1](#part-1--reproduce-figures-and-tables-no-hpc-required) above.

---

## Software versions

### Python environments

| Package | Version | Environment |
|---------|---------|-------------|
| Python | 3.10 | scanpy_env, scenic, scMultiomeGRN |
| Python | 3.12 | mmpbsa |
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
| celloracle | 0.20.0 | scMultiomeGRN |
| velocyto | 0.17.17 | scMultiomeGRN |
| torch | 2.4.0 | scMultiomeGRN |
| pydeseq2 | 0.5.4 | mmpbsa |
| openbabel | 3.1.1 | mmpbsa |
| meeko | 0.7.1 | mmpbsa (pip) |
| gmx_MMPBSA | 1.6.4 | mmpbsa (pip) |
| acpype | 2023.10.27 | mmpbsa |

### R environment

| Package | Version |
|---------|---------|
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
|------|---------|
| AutoDock Vina | 1.2 |
| GNINA | 1.3.2 |
| fpocket | 4.x |
| GROMACS | 2024.1 (gmx_mpi; TRUBA system module) |
| MEME Suite | 5.4.1 |

---

## Reproducibility notes

- All stochastic steps use fixed seeds. See manuscript §Random seeds for full table.
- **`fig_cellchat.py` Panel A** (chord diagram) requires `pdftoppm` (poppler-utils) to rasterize the CellChat chord PDF. This system tool is not included in the container; Panel A will be blank when run inside Apptainer. The full chord diagram PDF is available at `results/phase5/celloracle/` and is referenced directly in the manuscript figure.
- **GROMACS MD** (`gen_seed = -1`): velocity generation uses a system-generated random seed. RMSD values in the manuscript are not bit-reproducible without the archived trajectory files (stored on TRUBA ARF scratch). The trajectory files can be made available on request.
- **pySCENIC** uses 5 seeds (42, 1, 2, 3, 4) with a minimum consensus of 3/5. Results are consensus-stable across seeds.
- All scripts read the project root from `$DAM_DRUG_DIR`. Set this variable before running any script.
- Environment specs with exact version pins are in `envs/`.

---

## Data availability

Intermediate data files (21 files, 2.71 GB) are deposited at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19535759.svg)](https://doi.org/10.5281/zenodo.19535759)

| Dataset | Source | License |
|---------|--------|---------|
| Intermediate data (this study) | [Zenodo: 10.5281/zenodo.19535759](https://doi.org/10.5281/zenodo.19535759) | CC BY 4.0 |
| SEA-AD microglial atlas | [Allen Brain Cell Atlas](https://brain-map.org/consortia/sea-ad) — S3: `s3://sea-ad-single-cell-profiling/` | CC BY 4.0 |
| SEA-AD MTG (RNAseq + ATACseq) | Allen Brain Cell Atlas (same S3 bucket) | CC BY 4.0 |
| GSE95587 (bulk replication) | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95587) | Public |
| ChEMBL compound library | [ChEMBL](https://www.ebi.ac.uk/chembl/) (max_phase = 4, CNS filter) | CC BY-SA 3.0 |
| AlphaFold2 structures | [EBI AlphaFold Database v4](https://alphafold.ebi.ac.uk/) | CC BY 4.0 |
| PDB structures | [RCSB PDB](https://www.rcsb.org/) | Open |
| HOCOMOCO v11 motifs | [HOCOMOCO](https://hocomoco11.autosome.org/) | CC BY 4.0 |

---

## Citation

```
Özkurt Ç. Integrative Single-Cell Analysis of Alzheimer's Disease Microglia Identifies IKZF1
as a Late-Stage Neuroinflammatory Regulator and Generates Low-Confidence Computational
Hypotheses for Tafamidis and Diflunisal as Repurposing Candidates. (manuscript in preparation)
https://doi.org/10.5281/zenodo.19535749
```

---

## License

Code is released under the MIT License. See `LICENSE`.
