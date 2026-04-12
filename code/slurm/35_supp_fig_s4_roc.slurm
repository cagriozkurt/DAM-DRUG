#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=supp_s4_roc
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=8G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 6 — Supp Fig S4: Docking Enrichment ROC Curves
# Requires: results/phase4/roc/{target}_roc.csv (one per docking target)
# Submit: sbatch code/slurm/35_supp_fig_s4_roc.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

mkdir -p "$PROJDIR/logs" "$PROJDIR/results/figures"

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Check roc/ directory exists
ROC_DIR=$PROJDIR/results/phase4/roc
if [ ! -d "$ROC_DIR" ] || [ -z "$(ls -A "$ROC_DIR"/*.csv 2>/dev/null)" ]; then
  echo "ERROR: No ROC CSV files found in $ROC_DIR"
  echo "Run docking enrichment analysis first to generate {target}_roc.csv files."
  exit 1
fi

echo "ROC files found:"
ls "$ROC_DIR"/*.csv

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python "$PROJDIR/code/phase6_figures/supp_fig_s4_docking_roc.py"

echo "Done: $(date)"
