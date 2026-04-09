#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=prep_mtg
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/prep_mtg-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/prep_mtg-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
export DAM_DRUG_DIR=$PROJDIR

SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://mozkurt/dam-drug-scanpy:latest

echo "=== MTG Preprocessing (Apptainer) ==="
echo "Start: $(date)"

apptainer exec --cleanenv --bind /arf:/arf \
    "$SIF" \
    conda run -n scanpy_env python $PROJDIR/code/phase2_LR/36_prep_mtg_for_cellchat.py

echo "Done: $(date)"
