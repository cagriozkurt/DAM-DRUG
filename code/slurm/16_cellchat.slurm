#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=cellchat
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# Step 2 of 2: CellChat + NicheNet R analysis
# Prereq: 14_install_r_cellchat.slurm must have completed successfully
#         15_prep_mtg.slurm must have completed (produces prep/ directory)
#
# Submit: sbatch code/slurm/16_cellchat.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR

SIF=$PROJDIR/containers/dam-drug-r.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-r:latest

mkdir -p "$PROJDIR/logs"

echo "=== CellChat + NicheNet ==="
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Verify key packages available before submitting full job
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  Rscript -e "
for (pkg in c('CellChat','Matrix','hdf5r')) {
    ok <- requireNamespace(pkg, quietly=TRUE)
    cat(sprintf('%-20s %s\n', pkg, if(ok) 'OK' else 'MISSING'))
    if (!ok) stop(paste(pkg, 'not found in dam-drug-r image'))
}
"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  Rscript --vanilla "$PROJDIR/code/phase2_LR/02_cellchat_nichechat.R"

echo ""
echo "=== Done: $(date) ==="
