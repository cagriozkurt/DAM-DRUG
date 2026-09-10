#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_prep_regions
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=180G
#SBATCH --array=0-8
# ─────────────────────────────────────────────────────────────────────────────
# WP2.2 — Per-region CellChat input prep. Array index selects the region.
# Reads data/raw/SEA-AD/regions/SEAAD_<R>_RNAseq_final-nuclei.h5ad (from s3_01).
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
mkdir -p "$PROJDIR/logs"

REGIONS=(AnG DFC FI HIP ITG LEC MEC STG V1C)
R=${REGIONS[$SLURM_ARRAY_TASK_ID]}
H5AD=$PROJDIR/data/raw/SEA-AD/regions/SEAAD_${R}_RNAseq_final-nuclei.h5ad
[ -s "$H5AD" ] || { echo "ERROR: $H5AD missing (run s3_01, fill FIXMEs)"; exit 1; }

echo "[$R] job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID — start $(date)"
apptainer exec --bind "$PROJDIR:$PROJDIR" "$SIF" \
  conda run -n scenic python \
  "$PROJDIR/code/phase7_robustness/s3_cellchat/s3_prep.py" "$R" "$H5AD"
echo "[$R] done $(date)"
