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
PROJECT_ROOT="${DAM_DRUG_DIR:-/Volumes/PortableSSD/untitled folder/DAM-DRUG}"
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

echo "=== [1/2] Downloading SEA-AD Microglia pre-release (~3.3 GB) ==="
aws s3 cp \
  "s3://sea-ad-single-cell-profiling/Microglia-and-Immune-for-AAIC/SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC-pre-release.2025-07-24.h5ad" \
  "$MG_DIR/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad" \
  --no-sign-request \
  --only-show-errors

echo "  Done: $(du -sh "$MG_DIR/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad" | cut -f1)"

# =============================================================================
# 2. SEA-AD — MTG RNAseq full object (SECONDARY — all cell types, for context)
# =============================================================================
# Needed for: (a) confirming microglial isolation from full atlas;
#             (b) cell-cell communication (CellChat) with non-microglia
# Size: ~36 GB — download overnight
# SKIP this if storage is limited; microglia pre-release is sufficient for
# all primary DAM-DRUG analyses.

MTG_DIR="$RAW/SEA-AD"

echo ""
echo "=== [2/2] Downloading SEA-AD MTG full RNAseq (~36 GB) — may take 1-2 hours ==="
echo "    (Skip with Ctrl-C if storage is limited; microglia pre-release is sufficient)"
aws s3 cp \
  "s3://sea-ad-single-cell-profiling/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad" \
  "$MTG_DIR/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad" \
  --no-sign-request \
  --only-show-errors

echo "  Done: $(du -sh "$MTG_DIR/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad" | cut -f1)"

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
