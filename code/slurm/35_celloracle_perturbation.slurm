#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=celloracle_ko
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/celloracle_ko-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/celloracle_ko-%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 5.7: CellOracle in-silico TF perturbation
#
# Steps (each checkpoint-skipped if output already exists):
#   1. Load microglia_trajectory.h5ad, stratified 20K subsample, raw counts
#   2. Build TFdict from pySCENIC adj_matrix_aggregated.tsv + regulons
#   3. Oracle: PCA + KNN imputation
#   4. Fit GRN per cluster (ridge, alpha=10)
#   5. Simulate IKZF1 / IRF8 / SPI1 KO → quiver + grid plots + Markov table
#
# Prereq: scMultiomeGRN env with celloracle 0.20.0 + setuptools==60
# Fix (one-time, already done interactively):
#   conda run -n scMultiomeGRN pip install --force-reinstall setuptools==60.0.0
#
# Submit: sbatch code/slurm/35_celloracle_perturbation.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
export DAM_DRUG_DIR=$PROJDIR

mkdir -p $PROJDIR/logs

echo "=== CellOracle TF Perturbation ==="
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

module load miniconda3

# Ensure setuptools==60 (provides pkg_resources needed by gimmemotifs/celloracle)
conda run -n scMultiomeGRN pip install -q "setuptools==60.0.0" --force-reinstall 2>/dev/null || true

conda run -n scMultiomeGRN \
    python $PROJDIR/code/phase5_validation/35_celloracle_perturbation.py

echo ""
echo "=== Done: $(date) ==="
