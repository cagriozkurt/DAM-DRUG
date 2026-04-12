#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=reg_pt_corr
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — Regulon–pseudotime Spearman correlation (item 2.7)
# Submit: sbatch code/slurm/09_regulon_pseudotime_corr.slurm
# Or run interactively:
#   apptainer exec containers/scenic.sif python code/phase2_GRN/03_regulon_pseudotime_correlation.py
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

export DAM_DRUG_DIR=$PROJDIR

mkdir -p "$PROJDIR/logs"
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python "$PROJDIR/code/phase2_GRN/03_regulon_pseudotime_correlation.py"

echo "Done: $(date)"
