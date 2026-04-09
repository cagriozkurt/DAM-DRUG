#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=fig1_atlas
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/fig1_atlas-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/fig1_atlas-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=64G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 6 — Figure 1: Microglial Atlas & Substate Landscape
# Submit: sbatch code/slurm/29_fig1_atlas.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

export DAM_DRUG_DIR=$PROJDIR
mkdir -p $PROJDIR/logs $PROJDIR/results/figures

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scenic \
    pip install -q --user python-igraph leidenalg

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scenic \
    python $PROJDIR/code/phase6_figures/fig1_atlas.py

echo "Done: $(date)"
