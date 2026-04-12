#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=build-sif
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20        # barbun minimum
#SBATCH --mem=32G

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
mkdir -p "$PROJDIR/logs" "$PROJDIR/containers"

echo "Building SIF on $SLURMD_NODENAME: $(date)"
apptainer build --force "$PROJDIR/containers/scenic.sif" "$PROJDIR/containers/scanpy-env"
echo "Build exit code: $?"
ls -lh "$PROJDIR/containers/"
