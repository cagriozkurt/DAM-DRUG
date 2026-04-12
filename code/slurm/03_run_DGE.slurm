#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=DAM-DGE
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20        # barbun minimum = 20 cores/node
#SBATCH --mem=180G         # barbun has 384 GB/node; use 180 GB to leave headroom

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest
H5AD_SRC=$PROJDIR/data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad
H5AD_TMP=/tmp/SEA-AD_Microglia.h5ad

mkdir -p "$PROJDIR/logs" "$PROJDIR/results/phase1/DGE"

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Copy h5ad to local NVMe /tmp for faster I/O
echo "Copying h5ad to /tmp NVMe..."
cp "$H5AD_SRC" "$H5AD_TMP"
echo "Copy done: $(date)"

export DAM_DRUG_DIR=$PROJDIR
export SEA_AD_H5AD=$H5AD_TMP   # script will use this if set

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  --bind /tmp:/tmp \
  "$SIF" \
  conda run -n scanpy_env python "$PROJDIR/code/phase1_QC/03_DGE_microglial_states.py"

echo "Done: $(date)"

# Clean up /tmp
rm -f "$H5AD_TMP"
