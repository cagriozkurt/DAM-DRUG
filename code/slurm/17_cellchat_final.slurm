#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=cellchat_final
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR

SIF=$PROJDIR/containers/dam-drug-r.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-r:latest

mkdir -p "$PROJDIR/logs"

echo "=== SEA-AD CellChat Analysis (Apptainer) ==="
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  Rscript --vanilla "$PROJDIR/code/phase2_LR/02_cellchat_nichechat.R"

echo "Analysis complete: $(date)"
