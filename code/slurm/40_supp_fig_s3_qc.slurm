#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=supp_s3_qc
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/supp_s3_qc-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/supp_s3_qc-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=64G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 6 — Supp Fig S3: QC & Batch Correction
# Submit: sbatch code/slurm/40_supp_fig_s3_qc.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://mozkurt/dam-drug-scanpy:latest

export DAM_DRUG_DIR=$PROJDIR
mkdir -p $PROJDIR/logs $PROJDIR/results/figures

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scenic \
    python $PROJDIR/code/phase6_figures/supp_fig_s3_qc.py

echo "Done: $(date)"
