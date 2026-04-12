#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=supp_s3_qc
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=64G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 6 — Supp Fig S3: QC & Batch Correction
# Submit: sbatch code/slurm/34_supp_fig_s3_qc.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

mkdir -p "$PROJDIR/logs" "$PROJDIR/results/figures"

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  conda run -n scanpy_env python "$PROJDIR/code/phase6_figures/09_supp_fig_s3_qc.py"

echo "Done: $(date)"
