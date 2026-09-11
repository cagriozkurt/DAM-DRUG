#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=fig1_atlas
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=64G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 6 — Figure 1: Microglial Atlas & Substate Landscape
# Submit: sbatch code/slurm/10_fig1_atlas.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy@sha256:a16bd5bee57c5fb1918f2e4206368d342861392126c3a4159a462ddd8910efc3

export DAM_DRUG_DIR=$PROJDIR
mkdir -p "$PROJDIR/logs" "$PROJDIR/results/figures"

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  conda run -n scanpy_env python "$PROJDIR/code/phase6_figures/01_fig1_atlas.py"

echo "Done: $(date)"
