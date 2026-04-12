#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=phase3-struct
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=16G

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
CONTAINER=$PROJDIR/containers/fpocket-env.sif
[ -f "$CONTAINER" ] || apptainer pull "$CONTAINER" docker://ghcr.io/cagriozkurt/dam-drug-fpocket:latest

mkdir -p "$PROJDIR/logs" \
         "$PROJDIR/data/structures/af2" \
         "$PROJDIR/data/structures/pdb" \
         "$PROJDIR/data/structures/trimmed" \
         "$PROJDIR/data/structures/prepared" \
         "$PROJDIR/results/phase3/fpocket"

export DAM_DRUG_DIR=$PROJDIR

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Steps 1-2 (fetch + trim) run on login node (no internet on compute nodes).
# This job handles steps 3-6 only (prepare + fpocket + parse + summarize).

# Step 3: PDBFixer structure preparation (add missing atoms, H, minimize)
echo "=== Step 3: prepare (PDBFixer) ==="
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$CONTAINER" \
  python -u "$PROJDIR/code/phase3_structure/03_prepare_structures.py"
echo "prepare done: $(date)"

# Step 4: run fpocket on all prepared structures
echo "=== Step 4: fpocket ==="
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$CONTAINER" \
  bash "$PROJDIR/code/phase3_structure/04_run_fpocket.sh"
echo "fpocket done: $(date)"

# Step 5: parse fpocket output
echo "=== Step 5: parse ==="
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$CONTAINER" \
  python -u "$PROJDIR/code/phase3_structure/05_parse_fpocket.py"
echo "parse done: $(date)"

# Step 6: druggability ranking + plots
echo "=== Step 6: summarize ==="
apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$CONTAINER" \
  python -u "$PROJDIR/code/phase3_structure/06_druggability_summary.py"
echo "summarize done: $(date)"

echo "Done: $(date)"
