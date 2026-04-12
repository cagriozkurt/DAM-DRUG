#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=reparse-fpocket
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=16G

# Re-run fpocket (step 4) + parse + summarize (steps 5-6).
# 04_run_fpocket.sh skips structures whose _out dir already exists.
# fpocket failed silently on SPI1/ACSL1/PIK3CA in job 5504720 — re-running to fix.
# Prepared PDB files already exist; this only processes the missing 5 structures.

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
CONTAINER=$PROJDIR/containers/fpocket-env.sif
[ -f "$CONTAINER" ] || apptainer pull "$CONTAINER" docker://ghcr.io/cagriozkurt/dam-drug-fpocket:latest

export DAM_DRUG_DIR=$PROJDIR

mkdir -p "$PROJDIR/logs"

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Step 4: fpocket (skips structures whose _out dir already exists)
echo "=== Step 4: fpocket ==="
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$CONTAINER" \
  bash "$PROJDIR/code/phase3_structure/04_run_fpocket.sh"
echo "fpocket done: $(date)"

# Step 5: parse fpocket output (now includes SPI1/ACSL1/PIK3CA in STEM_MAP)
echo "=== Step 5: parse ==="
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$CONTAINER" \
  python -u "$PROJDIR/code/phase3_structure/05_parse_fpocket.py"
echo "parse done: $(date)"

# Step 6: druggability ranking + plots (now includes SPI1/ACSL1/PIK3CA in TF_ORDER)
echo "=== Step 6: summarize ==="
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$CONTAINER" \
  python -u "$PROJDIR/code/phase3_structure/06_druggability_summary.py"
echo "summarize done: $(date)"

echo "Done: $(date)"
