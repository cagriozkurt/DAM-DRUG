#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=DAM-explore
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/explore-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/explore-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=64G

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://mozkurt/dam-drug-scanpy:latest
H5AD_SRC=$PROJDIR/data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad
H5AD_TMP=/tmp/SEA-AD_Microglia.h5ad

mkdir -p $PROJDIR/logs $PROJDIR/results/phase1/explore

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

cp $H5AD_SRC $H5AD_TMP

export DAM_DRUG_DIR=$PROJDIR
export SEA_AD_H5AD=$H5AD_TMP

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  --bind /tmp:/tmp \
  "$SIF" \
  conda run -n scanpy_env python $PROJDIR/code/phase1_QC/02_explore_microglia_object.py

echo "Done: $(date)"
rm -f $H5AD_TMP
exit
