#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=traj_paga
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 1 — PAGA trajectory inference + diffusion pseudotime
#
# Outputs:
#   results/phase1/trajectory/pseudotime.csv
#   results/phase1/trajectory/paga_graph.pdf
#   results/phase1/trajectory/umap_paga_pseudotime.pdf
#   results/phase1/trajectory/microglia_trajectory.h5ad
#
# Submit: sbatch code/slurm/05_trajectory_paga.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest
H5AD_SRC=$PROJDIR/data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad
H5AD_TMP=/tmp/SEA-AD_Microglia.h5ad

mkdir -p "$PROJDIR/logs" "$PROJDIR/results/phase1/trajectory"

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Copy h5ad to NVMe /tmp for fast I/O
echo "Copying h5ad to /tmp..."
cp "$H5AD_SRC" "$H5AD_TMP"
echo "Copy done: $(date)"

export DAM_DRUG_DIR=$PROJDIR
export SEA_AD_H5AD=$H5AD_TMP

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  --bind /tmp:/tmp \
  "$SIF" \
  python "$PROJDIR/code/phase1_QC/04_trajectory_paga.py"

echo "Done: $(date)"
rm -f "$H5AD_TMP"
