#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=phase3-struct
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/phase3-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/phase3-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=16G

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
CONTAINER=$HOME/containers/fpocket-env.sif

mkdir -p $PROJDIR/logs $PROJDIR/data/structures/{af2,pdb,trimmed,prepared}
mkdir -p $PROJDIR/results/phase3/fpocket

export DAM_DRUG_DIR=$PROJDIR

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Steps 1-2 (fetch + trim) run on login node (no internet on compute nodes).
# This job handles steps 3-6 only (prepare + fpocket + parse + summarize).

# Step 3: PDBFixer structure preparation (add missing atoms, H, minimize)
echo "=== Step 3: prepare (PDBFixer) ==="
apptainer exec \
    --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
    $CONTAINER \
    /opt/conda/bin/python -u $PROJDIR/code/phase3_structure/08_prepare_structures.py
echo "prepare done: $(date)"

# Step 4: run fpocket on all prepared structures
echo "=== Step 4: fpocket ==="
apptainer exec \
    --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
    $CONTAINER \
    bash $PROJDIR/code/phase3_structure/09_run_fpocket.sh
echo "fpocket done: $(date)"

# Step 5: parse fpocket output
echo "=== Step 5: parse ==="
apptainer exec \
    --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
    $CONTAINER \
    /opt/conda/bin/python -u $PROJDIR/code/phase3_structure/10_parse_fpocket.py
echo "parse done: $(date)"

# Step 6: druggability ranking + plots
echo "=== Step 6: summarize ==="
apptainer exec \
    --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
    $CONTAINER \
    /opt/conda/bin/python -u $PROJDIR/code/phase3_structure/11_druggability_summary.py
echo "summarize done: $(date)"

echo "Done: $(date)"
