#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=grn_seed
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/grn_seed-%A_%a.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/grn_seed-%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=180G
#SBATCH --array=1-4
# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — GRNBoost2 multi-seed array (seeds 1–4)
#
# seed=42 already ran via run_pySCENIC.slurm → results/phase2/GRN/adj_matrix.tsv
# This array adds seeds 1–4, outputting adj_matrix_seed{1..4}.tsv
# After all 4 complete, run 25_aggregate_grn.py to average all 5 seeds.
#
# Submit:
#   sbatch code/slurm/25_grn_multiseed.slurm
# After completion:
#   conda run -n scenic python code/phase2_GRN/25_aggregate_grn.py
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

SEED=$SLURM_ARRAY_TASK_ID          # 1, 2, 3, or 4
PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://mozkurt/dam-drug-scanpy:latest

GRN_DIR=$PROJDIR/results/phase2/GRN
LOOM=$GRN_DIR/microglia_raw.loom
TF_LIST=$PROJDIR/data/resources/tf_lists/allTFs_hg38.txt
OUT_ADJ=$GRN_DIR/adj_matrix_seed${SEED}.tsv

echo "Job: $SLURM_JOB_ID  array=$SLURM_ARRAY_TASK_ID  seed=$SEED"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

[ -f "$LOOM" ] || { echo "ERROR: loom not found: $LOOM"; exit 1; }
[ -f "$TF_LIST" ] || { echo "ERROR: TF list not found: $TF_LIST"; exit 1; }

if [ -f "$OUT_ADJ" ] && [ -s "$OUT_ADJ" ]; then
    echo "Output already exists ($OUT_ADJ) — skipping."
    echo "Done: $(date)"
    exit 0
fi

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scenic \
    pyscenic grn \
      "$LOOM" \
      "$TF_LIST" \
      --output "$OUT_ADJ" \
      --num_workers 20 \
      --seed "$SEED"

echo ""
echo "Seed $SEED done: $(date)"
echo "Output: $OUT_ADJ  ($(wc -l < "$OUT_ADJ") edges)"
