#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=DAM-SCENIC
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/SCENIC-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/SCENIC-%j.err
#SBATCH --time=72:00:00          # 3-day limit; GRNBoost2+RcisTarget can be slow
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=180G

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest
H5AD_SRC=$PROJDIR/data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad
H5AD_TMP=/tmp/SEA-AD_Microglia.h5ad

mkdir -p $PROJDIR/logs $PROJDIR/results/phase2/GRN

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Copy h5ad to NVMe /tmp for fast I/O
cp $H5AD_SRC $H5AD_TMP
echo "h5ad copied to /tmp: $(date)"

export DAM_DRUG_DIR=$PROJDIR
export SEA_AD_H5AD=$H5AD_TMP
export SLURM_NTASKS=20   # passed to script for GRNBoost2 worker count

# First: download resource files if missing
MOTIF_DIR=$PROJDIR/data/resources/motifs
TF_DIR=$PROJDIR/data/resources/tf_lists
mkdir -p $MOTIF_DIR $TF_DIR

if [ ! -f "$MOTIF_DIR/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" ]; then
    echo "Downloading rankings database (~1.5 GB)..."
    wget -q -P $MOTIF_DIR \
      https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
fi

if [ ! -f "$MOTIF_DIR/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" ]; then
    echo "Downloading motif table..."
    wget -q -P $MOTIF_DIR \
      https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
fi

if [ ! -f "$TF_DIR/allTFs_hg38.txt" ]; then
    echo "Downloading TF list..."
    wget -q -P $TF_DIR \
      https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
fi

echo "Resources ready: $(date)"

apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  --bind /tmp:/tmp \
  "$SIF" \
  conda run -n scenic python $PROJDIR/code/phase2_GRN/05_pySCENIC_GRN.py

echo "Done: $(date)"
rm -f $H5AD_TMP
exit
