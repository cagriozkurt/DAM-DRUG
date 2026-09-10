#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_cc_collate
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=64G
# ─────────────────────────────────────────────────────────────────────────────
# WP2.4 — Collate the 9 per-region CellChat runs; rank SLIT2->ROBO2; run the
# expression-matched 1,000-pair permutation null; write the pan-cortical vs
# MTG-restricted verdict.
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
mkdir -p "$PROJDIR/logs"

echo "Start $(date)"
apptainer exec --bind "$PROJDIR:$PROJDIR" "$SIF" \
  conda run -n scenic python \
  "$PROJDIR/code/phase7_robustness/s3_cellchat/s3_collate.py"
echo "Done $(date)"
