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

# 1. Total Isolation: Wipe R paths so it only sees the container's internals
export R_LIBS_USER=""
export R_LIBS_SITE=""
export R_LIBS=""

echo "=== SEA-AD CellChat Analysis (Apptainer) ==="
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

# 2. Execute using the verified image
# We bind /arf so the container can read the 34GB counts_raw.h5
apptainer exec --cleanenv --nv --bind /arf:/arf \
    $PROJDIR/cellchat.sif \
    /opt/conda/envs/cellchat_r/bin/Rscript --vanilla $PROJDIR/code/phase2_LR/37_cellchat_nichechat.R

echo "Analysis complete: $(date)"
