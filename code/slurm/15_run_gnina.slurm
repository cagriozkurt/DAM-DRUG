#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=gnina-rescore
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/gnina-%A_%a.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/gnina-%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=8G
#SBATCH --array=0-5

TARGETS=(
    "PPARG_1FM9_LBD_prep"
    "IRF8_AF2_DBD_prep"
    "MAF_4EOT_bZIP_prep"
    "RUNX1_1LJM_Runt_prep"
    "BHLHE41_AF2_bHLH_prep"
    "IKZF1_8RQC_CRBN"
)

STEM=${TARGETS[$SLURM_ARRAY_TASK_ID]}
PROJDIR=/arf/scratch/mozkurt/DAM-DRUG

export DAM_DRUG_DIR=$PROJDIR
export GNINA=$HOME/apps/gnina/gnina

mkdir -p $PROJDIR/logs $PROJDIR/results/phase4/gnina/$STEM

echo "Job:    $SLURM_JOB_ID  Array: $SLURM_ARRAY_JOB_ID[$SLURM_ARRAY_TASK_ID]"
echo "Node:   $SLURMD_NODENAME"
echo "Target: $STEM"
echo "Start:  $(date)"

python3 -u \
    $PROJDIR/code/phase4_docking/15_run_gnina.py \
    $STEM \
    --cpus 20 \
    --top 30

echo "Done:   $(date)"
