#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=atac_da
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/atac_da-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/atac_da-%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 2.8 (revised): ATAC differential accessibility + TF motif enrichment
#
# Steps (resumable — each skipped if output already exists):
#   1. Filter ATAC to microglia, pseudobulk by donor
#   2. Wilcoxon DA test: AD vs Neurotypical donors
#   3. BED of AD-upregulated DAPs
#   4. Extract peak sequences (hg38, GetSequence.pl)
#   5. FIMO scan with 9 target TF motifs
#   6. Enrichment summary
#
# Submit: sbatch code/slurm/33_atac_da_validation.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-/arf/scratch/mozkurt/DAM-DRUG}
export DAM_DRUG_DIR=$PROJDIR

mkdir -p $PROJDIR/logs

echo "=== ATAC DA Validation ==="
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

SIF=$HOME/containers/dam-drug-scmultiomegrn.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scmultiomegrn:latest

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scMultiomeGRN \
    python $PROJDIR/code/phase2_GRN/32_atac_da_validation.py

echo ""
echo "=== Done: $(date) ==="
