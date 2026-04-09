#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=cellchat_final
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/cc_final-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/cc_final-%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
export DAM_DRUG_DIR=$PROJDIR

SIF=$HOME/containers/dam-drug-r.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://mozkurt/dam-drug-r:latest

echo "=== SEA-AD CellChat Analysis (Apptainer) ==="
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

apptainer exec --cleanenv --bind /arf:/arf \
    "$SIF" \
    Rscript --vanilla $PROJDIR/code/phase2_LR/37_cellchat_nichechat.R

echo "Analysis complete: $(date)"
