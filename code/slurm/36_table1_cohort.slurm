#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=table1
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=00:15:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# Table 1 — Cohort demographics from microglia_trajectory.h5ad
# Submit: sbatch code/slurm/36_table1_cohort.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest
mkdir -p "$PROJDIR/logs" "$PROJDIR/results/tables"

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  conda run -n scanpy_env python "$PROJDIR/code/tables/make_table1_cohort.py"

echo "Done: $(date)"
