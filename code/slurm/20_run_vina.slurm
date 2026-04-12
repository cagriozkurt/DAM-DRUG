#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=vina-dock
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=8G
#SBATCH --array=0-5

# ── Target list (index matches SLURM_ARRAY_TASK_ID) ─────────────────────────
# Index 5 (IKZF1_8RQC_CRBN) requires 02_ikzf1_protac_prep.py to have been run first.
TARGETS=(
    "PPARG_1FM9_LBD_prep"
    "IRF8_AF2_DBD_prep"
    "MAF_4EOT_bZIP_prep"
    "RUNX1_1LJM_Runt_prep"
    "BHLHE41_AF2_bHLH_prep"
    "IKZF1_8RQC_CRBN"
)

set -eo pipefail

STEM=${TARGETS[$SLURM_ARRAY_TASK_ID]}

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

export DAM_DRUG_DIR=$PROJDIR
# Override VINA by setting it in your environment before submitting, or place the
# binary at $PROJDIR/bin/vina.
export VINA=${VINA:-$PROJDIR/bin/vina}

mkdir -p "$PROJDIR/logs" "$PROJDIR/results/phase4/vina/$STEM"

echo "Job:    $SLURM_JOB_ID  Array: $SLURM_ARRAY_JOB_ID[$SLURM_ARRAY_TASK_ID]"
echo "Node:   $SLURMD_NODENAME"
echo "Target: $STEM"
echo "Start:  $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python -u "$PROJDIR/code/phase4_docking/04_run_vina.py" \
    "$STEM" \
    --cpus 20

echo "Done:   $(date)"
