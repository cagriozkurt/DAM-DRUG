#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_core_rmsd
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# WP3.4 — Ligand core-RMSD (10% lowest-RMSF heavy atoms, final 20 ns) per
# replica for the 3 glue ternary complexes. Advancement gate: mean < 3.5 A in
# all 8 replicas.
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif    # any env with MDAnalysis; add to the SIF if absent
mkdir -p "$PROJDIR/logs"

echo "Start $(date)"
apptainer exec --bind "$PROJDIR:$PROJDIR" "$SIF" \
  conda run -n scenic python \
  "$PROJDIR/code/phase7_robustness/s4_md/s4_core_rmsd.py"
echo "Done $(date)"
