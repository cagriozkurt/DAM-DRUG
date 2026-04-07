#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=table1
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/table1-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/table1-%j.err
#SBATCH --time=00:15:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# Table 1 — Cohort demographics from microglia_trajectory.h5ad
# Submit: sbatch code/slurm/42_table1_cohort.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
CONTAINER=$PROJDIR/scanpy.sif
export DAM_DRUG_DIR=$PROJDIR
mkdir -p $PROJDIR/logs $PROJDIR/results/tables

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  $CONTAINER \
  conda run -n scanpy_env \
    python $PROJDIR/code/tables/make_table1_cohort.py

echo "Done: $(date)"
