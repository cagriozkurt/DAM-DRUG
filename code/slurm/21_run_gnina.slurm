#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=gnina-rescore
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=8G
#SBATCH --array=0-5

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
# Override GNINA by setting it in your environment before submitting, or place the
# binary at $PROJDIR/bin/gnina.
export GNINA=${GNINA:-$PROJDIR/bin/gnina}

mkdir -p "$PROJDIR/logs" "$PROJDIR/results/phase4/gnina/$STEM"

echo "Job:    $SLURM_JOB_ID  Array: $SLURM_ARRAY_JOB_ID[$SLURM_ARRAY_TASK_ID]"
echo "Node:   $SLURMD_NODENAME"
echo "Target: $STEM"
echo "Start:  $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python -u "$PROJDIR/code/phase4_docking/05_run_gnina.py" \
    "$STEM" \
    --cpus 20 \
    --top 30

echo "Done:   $(date)"
