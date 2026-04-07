#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=traj_paga
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/traj_paga-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/traj_paga-%j.err
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
# Submit: sbatch code/slurm/27_trajectory_paga.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
CONTAINER=$HOME/containers/scanpy-env.sif
H5AD_SRC=$PROJDIR/data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad
H5AD_TMP=/tmp/SEA-AD_Microglia.h5ad

mkdir -p $PROJDIR/logs $PROJDIR/results/phase1/trajectory

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Copy h5ad to NVMe /tmp for fast I/O
echo "Copying h5ad to /tmp..."
cp $H5AD_SRC $H5AD_TMP
echo "Copy done: $(date)"

export DAM_DRUG_DIR=$PROJDIR
export SEA_AD_H5AD=$H5AD_TMP

# Install missing deps to ~/.local (home is mounted; container is read-only)
apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  --bind /tmp:/tmp \
  $CONTAINER \
  conda run -n scenic \
    pip install -q --user python-igraph leidenalg

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  --bind /tmp:/tmp \
  $CONTAINER \
  conda run -n scenic \
    python $PROJDIR/code/phase1_QC/27_trajectory_paga.py

echo "Done: $(date)"
rm -f $H5AD_TMP
