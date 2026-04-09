#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=supp_s4_roc
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/supp_s4_roc-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/supp_s4_roc-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=8G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 6 — Supp Fig S4: Docking Enrichment ROC Curves
# Requires: results/phase4/roc/{target}_roc.csv (one per docking target)
# Submit: sbatch code/slurm/41_supp_fig_s4_roc.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://mozkurt/dam-drug-scanpy:latest

export DAM_DRUG_DIR=$PROJDIR
mkdir -p $PROJDIR/logs $PROJDIR/results/figures

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Check roc/ directory exists
ROC_DIR=$PROJDIR/results/phase4/roc
if [ ! -d "$ROC_DIR" ] || [ -z "$(ls -A $ROC_DIR/*.csv 2>/dev/null)" ]; then
  echo "ERROR: No ROC CSV files found in $ROC_DIR"
  echo "Run docking enrichment analysis first to generate {target}_roc.csv files."
  exit 1
fi

echo "ROC files found:"
ls $ROC_DIR/*.csv

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scenic \
    python $PROJDIR/code/phase6_figures/supp_fig_s4_docking_roc.py

echo "Done: $(date)"
