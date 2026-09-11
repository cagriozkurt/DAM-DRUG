#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_cellchat
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=16:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G
#SBATCH --array=0-8
# ─────────────────────────────────────────────────────────────────────────────
# WP2.3 — CellChat v2 per region (array). Parameter-identical to the accepted
# MTG analysis. Reads results/phase7/cellchat_multiregion/<R>/prep/.
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/dam-drug-r.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-r@sha256:3d0db505d2c1caa014caed89f631b453e0eea13b1ad5919f3ffde92a607d0648
mkdir -p "$PROJDIR/logs"

REGIONS=(AnG DFC FI HIP ITG LEC MEC STG V1C)
R=${REGIONS[$SLURM_ARRAY_TASK_ID]}

echo "[$R] job $SLURM_JOB_ID — start $(date)"
apptainer exec --bind "$PROJDIR:$PROJDIR" "$SIF" \
  Rscript --vanilla "$PROJDIR/code/phase7_robustness/s3_cellchat/s3_cellchat.R" "$R"
echo "[$R] done $(date)"
