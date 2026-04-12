#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=seaad-download
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=16G

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
OUTDIR=$PROJDIR/data/raw/SEA-AD
mkdir -p "$OUTDIR" "$PROJDIR/logs"

BASE=https://sea-ad-single-cell-profiling.s3.amazonaws.com

echo "Start: $(date)"

# -c = resume if interrupted; -q --show-progress = clean output
# 1. Microglia pre-release (3.0 GB) — primary analysis object
wget -c -q --show-progress -O \
  $OUTDIR/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad \
  "$BASE/Microglia-and-Immune-for-AAIC/SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC-pre-release.2025-07-24.h5ad"
echo "Microglia done: $(date)"

# 2. MTG full RNAseq (34 GB) — Phase 2 CellChat
wget -c -q --show-progress -O \
  $OUTDIR/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  "$BASE/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
echo "MTG RNAseq done: $(date)"

# 3. MTG ATACseq (18 GB) — Phase 2 scMultiomeGRN
wget -c -q --show-progress -O \
  $OUTDIR/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad \
  "$BASE/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
echo "MTG ATACseq done: $(date)"

echo "All downloads complete: $(date)"
ls -lh $OUTDIR
exit
