#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_dl_subclass
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20          # barbun requires 20 cores/node
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# WP2.1 — Download SEA-AD Multiregion 2026 per-SUBCLASS objects.
#
# Confirmed by anonymous S3 listing (2026-09-10):
#   bucket:  s3://sea-ad-single-cell-profiling  (public, no-sign-request)
#   subclass objects span ALL 10 regions in one file each:
#     Multiregion_2026/subclass_objects/SEAAD_<SUBCLASS>_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad
#   X = log-norm per 10k; RAW counts in .layers["UMIs"]; region col = "Brain Region".
#   The per-region folders (AnG/RNAseq/ ...) also exist but DFC/RNAseq/ is EMPTY —
#   the subclass objects are the reliable multi-region source, and include DFC.
#
# For the SLIT2->ROBO2 (interneuron -> microglia) CellChat we need:
#   receiver:  Immune (= Micro-PVM + lymphocytes)
#   senders:   all GABAergic subclasses (SLIT2+ interneurons)
#   context:   glia + a representative set of glutamatergic subclasses
# The 3 giant IT classes (L23-IT 41 GB, L4-IT 16 GB, L5-IT 15 GB) are marked
# OPTIONAL — CellChat is run on a 5,000-cell/type cap anyway.
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
OUT=$PROJDIR/data/raw/SEA-AD/multiregion_2026/subclass_objects
mkdir -p "$OUT" "$PROJDIR/logs"
BASE=https://sea-ad-single-cell-profiling.s3.amazonaws.com/Multiregion_2026/subclass_objects
SUFFIX=_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad

REQUIRED=(
  Immune                                              # receiver: microglia
  Sst Vip Pvalb Lamp5 Sncg Pax6 Chandelier Lamp5-Lhx6 Sst-Chodl   # senders: interneurons
  Astrocyte Oligodendrocyte OPC Endothelial VLMC-Perivascular      # glial/vascular context
  L5-ET L56-NP L6-CT L6b L6-IT-Car3 EC-IT CA2-4 DG Sub-CA1         # smaller glutamatergic context
)
OPTIONAL=( L23-IT L4-IT L5-IT L6-IT )   # huge; download only if you want full exc coverage

DL="${1:-required}"   # pass "all" as arg 1 to also fetch OPTIONAL
LIST=("${REQUIRED[@]}")
[ "$DL" = "all" ] && LIST+=("${OPTIONAL[@]}")

echo "Start $(date) — fetching ${#LIST[@]} subclass objects (mode: $DL)"
for SC in "${LIST[@]}"; do
    url="$BASE/SEAAD_${SC}${SUFFIX}"
    dst="$OUT/SEAAD_${SC}_multiregion.h5ad"
    if [ -s "$dst" ]; then echo "  [skip] $SC"; continue; fi
    echo "  [get ] $SC  <- $url"
    wget -c -q --show-progress -O "$dst" "$url"
done
echo "Done $(date)"
du -sh "$OUT"; ls -lh "$OUT"
