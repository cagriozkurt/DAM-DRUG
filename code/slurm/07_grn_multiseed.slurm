#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=grn_seed
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=180G
#SBATCH --array=1-4
# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — GRNBoost2 multi-seed array (seeds 1–4)
#
# seed=42 already ran via 06_run_pySCENIC.slurm → results/phase2/GRN/adj_matrix.tsv
# This array adds seeds 1–4, outputting adj_matrix_seed{1..4}.tsv
# After all 4 complete, run 02_aggregate_grn.py to average all 5 seeds.
#
# Submit:
#   sbatch code/slurm/07_grn_multiseed.slurm
# After completion:
#   apptainer exec containers/scenic.sif python code/phase2_GRN/02_aggregate_grn.py
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

SEED=$SLURM_ARRAY_TASK_ID          # 1, 2, 3, or 4
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

GRN_DIR=$PROJDIR/results/phase2/GRN
LOOM=$GRN_DIR/microglia_raw.loom
TF_LIST=$PROJDIR/data/resources/tf_lists/allTFs_hg38.txt
OUT_ADJ=$GRN_DIR/adj_matrix_seed${SEED}.tsv

mkdir -p "$PROJDIR/logs"

echo "Job: $SLURM_JOB_ID  array=$SLURM_ARRAY_TASK_ID  seed=$SEED"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

[ -f "$LOOM" ]    || { echo "ERROR: loom not found: $LOOM"; exit 1; }
[ -f "$TF_LIST" ] || { echo "ERROR: TF list not found: $TF_LIST"; exit 1; }

if [ -f "$OUT_ADJ" ] && [ -s "$OUT_ADJ" ]; then
    echo "Output already exists ($OUT_ADJ) — skipping."
    echo "Done: $(date)"
    exit 0
fi

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  --env LD_LIBRARY_PATH=/opt/conda/envs/scenic/lib \
  "$SIF" \
  pyscenic grn \
    "$LOOM" \
    "$TF_LIST" \
    --output "$OUT_ADJ" \
    --num_workers 20 \
    --seed "$SEED"

echo ""
echo "Seed $SEED done: $(date)"
echo "Output: $OUT_ADJ  ($(wc -l < "$OUT_ADJ") edges)"
