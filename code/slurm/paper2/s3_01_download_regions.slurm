#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_dl_regions
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
# ─────────────────────────────────────────────────────────────────────────────
# WP2.1 — Download the 9 non-MTG SEA-AD region snRNA objects (all cell types:
#         neurons + microglia needed for CellChat sender/receiver).
#
# The accepted paper only used MTG. SEA-AD's expanded release adds these regions;
# the microglia-only multi-regional object on the SSD lacks neurons.
#
# FIXME: confirm the exact S3 object keys. The MTG object lives at
#   $BASE/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad
# The expanded regions are expected under an analogous path; verify at
#   https://registry.opendata.aws/allen-sea-ad-atlas  /  the SEA-AD data portal.
# Region acronyms used in the microglia object: AnG DFC FI HIP ITG LEC MEC STG V1C
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
OUT=$PROJDIR/data/raw/SEA-AD/regions
mkdir -p "$OUT" "$PROJDIR/logs"
BASE=https://sea-ad-single-cell-profiling.s3.amazonaws.com

# FIXME: fill in the real object key per region (path + filename + date stamp).
declare -A OBJ=(
  [AnG]="FIXME/RNAseq/SEAAD_AnG_RNAseq_final-nuclei.FIXME.h5ad"
  [DFC]="FIXME/RNAseq/SEAAD_DFC_RNAseq_final-nuclei.FIXME.h5ad"
  [FI]="FIXME/RNAseq/SEAAD_FI_RNAseq_final-nuclei.FIXME.h5ad"
  [HIP]="FIXME/RNAseq/SEAAD_HIP_RNAseq_final-nuclei.FIXME.h5ad"
  [ITG]="FIXME/RNAseq/SEAAD_ITG_RNAseq_final-nuclei.FIXME.h5ad"
  [LEC]="FIXME/RNAseq/SEAAD_LEC_RNAseq_final-nuclei.FIXME.h5ad"
  [MEC]="FIXME/RNAseq/SEAAD_MEC_RNAseq_final-nuclei.FIXME.h5ad"
  [STG]="FIXME/RNAseq/SEAAD_STG_RNAseq_final-nuclei.FIXME.h5ad"
  [V1C]="FIXME/RNAseq/SEAAD_V1C_RNAseq_final-nuclei.FIXME.h5ad"
)

echo "Start $(date)"
for R in "${!OBJ[@]}"; do
    key="${OBJ[$R]}"
    if [[ "$key" == FIXME* ]]; then
        echo "SKIP $R — object key not set (FIXME)"; continue
    fi
    dst="$OUT/SEAAD_${R}_RNAseq_final-nuclei.h5ad"
    echo "[$R] $BASE/$key -> $dst"
    wget -c -q --show-progress -O "$dst" "$BASE/$key"
done
echo "Done $(date)"
ls -lh "$OUT"
