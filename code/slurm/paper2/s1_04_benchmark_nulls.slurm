#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_grn_benchmark
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# WP1.4 — Benchmark JASPAR-2026 rescued regulons + negative null models.
#   - per-regulon AUCell vs diffusion pseudotime (Spearman + BH)
#   - IKZF1 vs BHLHE41/IRF8/SPI1/RUNX1/CEBPB/PPARG trajectory-coupling ranking
#   - 1,000-permutation empirical FDR for |rho| > 0.30 (unrestricted + donor-block)
#   - IKZF1 vs IKZF2 vs IKZF3 paralogue disentangling
#   - hypergeometric enrichment of IKZF1 targets vs JASPAR 2026 curated TF-TG
#
# Needs: results/phase7/grn_jaspar2026/{scenic_auc_jaspar2026.loom,regulons_jaspar2026.csv}
#        data/resources/jaspar2026/tf_tg_curated.tsv   (FIXME: source the curated set)
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
mkdir -p "$PROJDIR/logs" "$PROJDIR/results/phase7/grn_jaspar2026"

echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME — start $(date)"
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  --env LD_LIBRARY_PATH=/opt/conda/envs/scenic/lib \
  "$SIF" \
  conda run -n scenic python "$PROJDIR/code/phase7_robustness/s1_grn/s1_benchmark.py"
echo "Done $(date)"
