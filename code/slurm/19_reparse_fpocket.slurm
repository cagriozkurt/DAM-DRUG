#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=reparse-fpocket
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/reparse-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/reparse-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=16G

# Re-run fpocket (step 4) + parse + summarize (steps 5-6).
# 09_run_fpocket.sh skips structures whose _out dir already exists.
# fpocket failed silently on SPI1/ACSL1/PIK3CA in job 5504720 — re-running to fix.
# Prepared PDB files already exist; this only processes the missing 5 structures.

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
CONTAINER=$HOME/containers/fpocket-env.sif

export DAM_DRUG_DIR=$PROJDIR

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Step 4: fpocket (skips structures whose _out dir already exists)
echo "=== Step 4: fpocket ==="
apptainer exec \
    --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
    $CONTAINER \
    bash $PROJDIR/code/phase3_structure/09_run_fpocket.sh
echo "fpocket done: $(date)"

# Step 5: parse fpocket output (now includes SPI1/ACSL1/PIK3CA in STEM_MAP)
echo "=== Step 5: parse ==="
apptainer exec \
    --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
    $CONTAINER \
    /opt/conda/bin/python -u $PROJDIR/code/phase3_structure/10_parse_fpocket.py
echo "parse done: $(date)"

# Step 6: druggability ranking + plots (now includes SPI1/ACSL1/PIK3CA in TF_ORDER)
echo "=== Step 6: summarize ==="
apptainer exec \
    --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
    $CONTAINER \
    /opt/conda/bin/python -u $PROJDIR/code/phase3_structure/11_druggability_summary.py
echo "summarize done: $(date)"

echo "Done: $(date)"
