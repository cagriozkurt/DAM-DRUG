#!/usr/bin/env bash
# =============================================================================
# DAM-DRUG Phase 1 — Data Acquisition
# Correct S3 bucket: s3://sea-ad-single-cell-profiling/ (NOT allen-brain-cell-atlas)
#
# REVISED STRATEGY: Download the pre-release microglia-specific object first.
# This contains 240,651 microglia nuclei across 10 brain regions with scVI
# latent space + UMAP + Supertype labels already computed — skipping most QC.
#
# Run from project root: bash code/phase1_QC/01_data_acquisition.sh
# =============================================================================

set -euo pipefail
PROJECT_ROOT="${DAM_DRUG_DIR:-$(pwd)}"
RAW="$PROJECT_ROOT/data/raw"

if ! command -v aws &>/dev/null; then
  echo "ERROR: aws CLI not found. Install with: brew install awscli"
  exit 1
fi

# =============================================================================
# 1. SEA-AD — Microglia-and-Immune pre-release (PRIMARY — download this first)
# =============================================================================
# 240,651 microglia nuclei; 10 brain regions; 84 donors
# Contains: raw UMI counts (layers), log-normalized X, X_scVI, X_umap,
#           Supertype labels (Micro-PVM_1/2/3/4), APOE/Braak/CERAD/sex/PMI
# Size: ~3.3 GB — download in ~5 minutes on fast connection
#
# NOTE: This is a pre-release for AAIC 2025. The full official release
# covering all cell types is expected Fall 2025. Microglia are already final.

MG_DIR="$RAW/SEA-AD"
mkdir -p "$MG_DIR"

MG_FILE="$MG_DIR/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad"
MG_SHA256="f7ba373c6554ba0f04a9a67704e2046e18608bbc6a435a5f71af4c164a1a421f"
echo "=== [1/2] Downloading SEA-AD Microglia pre-release (~3.3 GB) ==="
if [ -f "$MG_FILE" ] && echo "$MG_SHA256  $MG_FILE" | sha256sum -c --quiet 2>/dev/null; then
  echo "  Already exists and checksum OK, skipping."
else
  aws s3 cp \
    "s3://sea-ad-single-cell-profiling/Microglia-and-Immune-for-AAIC/SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC-pre-release.2025-07-24.h5ad" \
    "$MG_FILE" \
    --no-sign-request \
    --only-show-errors
  echo "$MG_SHA256  $MG_FILE" | sha256sum -c || { echo "ERROR: checksum mismatch for $MG_FILE"; exit 1; }
  echo "  Done: $(du -sh "$MG_FILE" | cut -f1)"
fi

# =============================================================================
# 2. SEA-AD — MTG RNAseq full object (SECONDARY — all cell types, for context)
# =============================================================================
# Needed for: (a) confirming microglial isolation from full atlas;
#             (b) cell-cell communication (CellChat) with non-microglia
# Size: ~36 GB — download overnight
# SKIP this if storage is limited; microglia pre-release is sufficient for
# all primary DAM-DRUG analyses.

MTG_DIR="$RAW/SEA-AD"

MTG_FILE="$MTG_DIR/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
MTG_SHA256="9c1b48266d0a9aef76ad20fb8487d604158b10fad18dc0de85d5261ef06cb7c8"
echo ""
echo "=== [2/2] Downloading SEA-AD MTG full RNAseq (~36 GB) — may take 1-2 hours ==="
echo "    (Skip with Ctrl-C if storage is limited; microglia pre-release is sufficient)"
if [ -f "$MTG_FILE" ] && echo "$MTG_SHA256  $MTG_FILE" | sha256sum -c --quiet 2>/dev/null; then
  echo "  Already exists and checksum OK, skipping."
else
  aws s3 cp \
    "s3://sea-ad-single-cell-profiling/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad" \
    "$MTG_FILE" \
    --no-sign-request \
    --only-show-errors
  echo "$MTG_SHA256  $MTG_FILE" | sha256sum -c || { echo "ERROR: checksum mismatch for $MTG_FILE"; exit 1; }
  echo "  Done: $(du -sh "$MTG_FILE" | cut -f1)"
fi

# =============================================================================
# DEFERRED: MTG ATACseq (18 GB) — download before Phase 2 only
# =============================================================================
# aws s3 cp \
#   "s3://sea-ad-single-cell-profiling/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad" \
#   "$RAW/SEA-AD/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad" \
#   --no-sign-request --only-show-errors

echo ""
echo "=== File summary ==="
du -sh "$RAW"/*/* 2>/dev/null || true
echo ""
echo "=== Acquisition complete. Next step: run 02_explore_microglia_object.py ==="
