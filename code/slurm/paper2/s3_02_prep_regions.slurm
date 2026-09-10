#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_prep_regions
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=180G
#SBATCH --array=0-8
# ─────────────────────────────────────────────────────────────────────────────
# WP2.2 — Per-region CellChat input prep from the Multiregion 2026 subclass
# objects (downloaded by s3_01). Array index selects the region.
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
mkdir -p "$PROJDIR/logs"

REGIONS=(AnG DFC FI HIP ITG LEC MEC STG V1C)
R=${REGIONS[$SLURM_ARRAY_TASK_ID]}
SC_DIR=$PROJDIR/data/raw/SEA-AD/multiregion_2026/subclass_objects
[ -d "$SC_DIR" ] && [ -n "$(ls -A "$SC_DIR" 2>/dev/null)" ] || \
    { echo "ERROR: $SC_DIR empty — run s3_01_download_subclass_objects.slurm"; exit 1; }

echo "[$R] job $SLURM_JOB_ID — start $(date)"
apptainer exec --bind "$PROJDIR:$PROJDIR" "$SIF" \
  conda run -n scenic python \
  "$PROJDIR/code/phase7_robustness/s3_cellchat/s3_prep.py" "$R" "$SC_DIR"
echo "[$R] done $(date)"
