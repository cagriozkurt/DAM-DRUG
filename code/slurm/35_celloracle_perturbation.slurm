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
# Submit: sbatch code/slurm/35_celloracle_perturbation.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-/arf/scratch/mozkurt/DAM-DRUG}
export DAM_DRUG_DIR=$PROJDIR

mkdir -p $PROJDIR/logs

echo "=== CellOracle TF Perturbation ==="
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

SIF=$HOME/containers/dam-drug-scmultiomegrn.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scmultiomegrn:latest

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scMultiomeGRN \
    python $PROJDIR/code/phase5_validation/35_celloracle_perturbation.py

echo ""
echo "=== Done: $(date) ==="
