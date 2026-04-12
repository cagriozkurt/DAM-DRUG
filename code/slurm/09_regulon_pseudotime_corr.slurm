#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=reg_pt_corr
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/reg_pt_corr-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/reg_pt_corr-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — Regulon–pseudotime Spearman correlation (item 2.7)
# Submit: sbatch code/slurm/28_regulon_pseudotime_corr.slurm
# Or run interactively: apptainer exec ... conda run -n scenic python ...
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

export DAM_DRUG_DIR=$PROJDIR

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"
mkdir -p $PROJDIR/logs

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scenic \
    python $PROJDIR/code/phase2_GRN/28_regulon_pseudotime_correlation.py

echo "Done: $(date)"
