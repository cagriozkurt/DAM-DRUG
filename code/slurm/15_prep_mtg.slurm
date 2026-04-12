#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=prep_mtg
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR

SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

mkdir -p "$PROJDIR/logs"

echo "=== MTG Preprocessing (Apptainer) ==="
echo "Start: $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  conda run -n scanpy_env python "$PROJDIR/code/phase2_LR/01_prep_mtg_for_cellchat.py"

echo "Done: $(date)"
