#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=scmgrn_env
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/scmgrn_env-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/scmgrn_env-%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# One-time setup for scMultiomeGRN (item 2.8) — RESUMABLE
#
# Steps (each skipped if already done):
#   1. Create conda env scMultiomeGRN
#   2. Extract HOCOMOCOv11 motifs
#   3. Download hg38 chromosome FASTAs
#   4. Export per-state MTX from h5ads
#
# Submit: sbatch code/slurm/31_scMultiomeGRN_env.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
TOOLDIR=$PROJDIR/tools/scMultiomeGRN/extracted/ScmultiomeGRN-main
DATADIR=$PROJDIR/tools/scMultiomeGRN/data_resource
CONTAINER=$HOME/containers/scanpy-env.sif
export ATAC_H5AD=$PROJDIR/data/raw/SEA-AD/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad
export DAM_DRUG_DIR=$PROJDIR

CKPT_DIR=$PROJDIR/results/phase2/scMultiomeGRN/.checkpoints
mkdir -p $PROJDIR/logs $CKPT_DIR

echo "=== scMultiomeGRN One-Time Setup ==="
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# ── 1. Conda env ──────────────────────────────────────────────────────────────
echo ""
echo "--- Step 1: Conda env scMultiomeGRN ---"
module load miniconda3

if [ -f "$CKPT_DIR/env.done" ]; then
    echo "SKIP (checkpoint exists)"
else
    if conda env list | grep -q "^scMultiomeGRN "; then
        echo "Env already exists."
    else
        conda create -y -n scMultiomeGRN python=3.10 \
            bedops "meme=5.4.1" perl-bioperl \
            -c bioconda -c conda-forge

        conda run -n scMultiomeGRN \
            conda install -y -c bioconda arboreto

        # PyTorch with CUDA 12.1 (V100-compatible for akya-cuda)
        conda run -n scMultiomeGRN pip install \
            torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 \
            --index-url https://download.pytorch.org/whl/cu121

        conda run -n scMultiomeGRN pip install \
            torch_geometric lightning h5py \
            pandas scikit-learn scipy tqdm gensim
    fi
    touch $CKPT_DIR/env.done
    echo "Done."
fi

# ── 2. HOCOMOCOv11 motifs ─────────────────────────────────────────────────────
echo ""
echo "--- Step 2: HOCOMOCOv11 motifs ---"
MOTIF_DIR=$DATADIR/HOCOMOCOv11

if [ -f "$CKPT_DIR/motifs.done" ]; then
    echo "SKIP (checkpoint exists): $(ls $MOTIF_DIR/*.meme 2>/dev/null | wc -l) .meme files"
else
    mkdir -p $DATADIR
    unzip -o -q $TOOLDIR/data_resource/HOCOMOCOv11.zip -d $DATADIR/
    NMEME=$(ls $MOTIF_DIR/*.meme 2>/dev/null | wc -l)
    echo "Extracted: $NMEME .meme files"
    [ "$NMEME" -gt 0 ] && touch $CKPT_DIR/motifs.done
fi

# ── 3. hg38 genome ────────────────────────────────────────────────────────────
echo ""
echo "--- Step 3: hg38 genome ---"
GENOME_DIR=$DATADIR/download_genome/hg38
mkdir -p $GENOME_DIR

if [ -f "$CKPT_DIR/genome.done" ]; then
    echo "SKIP (checkpoint exists): $(ls $GENOME_DIR/*.fa 2>/dev/null | wc -l) chromosomes"
else
    MISSING=0
    for CHR in $(seq 1 22) X Y M; do
        [ -f "$GENOME_DIR/chr${CHR}.fa" ] || MISSING=$((MISSING+1))
    done

    if [ $MISSING -gt 0 ]; then
        echo "Downloading $MISSING missing chromosome FASTAs..."
        cd $GENOME_DIR
        for CHR in $(seq 1 22) X Y M; do
            FA="chr${CHR}.fa"
            [ -f "$FA" ] && continue
            wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${FA}.gz" \
                 -O "${FA}.gz" && gunzip "${FA}.gz" &
        done
        wait
        echo "Genome ready: $(ls $GENOME_DIR/*.fa | wc -l) chromosomes"
    else
        echo "All chromosomes present."
    fi
    touch $CKPT_DIR/genome.done
fi

# ── 4. Export per-state MTX ───────────────────────────────────────────────────
echo ""
echo "--- Step 4: Export per-state MTX ---"

# Check if all 6 states are already done
STATES=("Homeostatic" "DAM" "DAM_IRM" "IRM" "LateAD_DAM" "LAM")
ALL_DONE=1
for TAG in "${STATES[@]}"; do
    MTX=$PROJDIR/results/phase2/scMultiomeGRN/input/$TAG/scrna/matrix.mtx
    AMTX=$PROJDIR/results/phase2/scMultiomeGRN/input/$TAG/atac/matrix.mtx
    if [ ! -f "$MTX" ] || [ ! -f "$AMTX" ]; then
        ALL_DONE=0
        break
    fi
done

if [ $ALL_DONE -eq 1 ]; then
    echo "SKIP — all states already exported"
else
    if [ ! -f "$ATAC_H5AD" ]; then
        echo "ERROR: ATAC h5ad not found: $ATAC_H5AD"
        exit 1
    fi
    echo "ATAC h5ad: $ATAC_H5AD"
    echo "Memory: $(free -h | awk '/^Mem/ {print $2}') total"

    apptainer exec \
      --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
      $CONTAINER \
      conda run -n scenic \
        python $PROJDIR/code/phase2_GRN/31_prep_scMultiomeGRN_input.py
fi

echo ""
echo "=== Setup complete: $(date) ==="
echo "Next: sbatch code/slurm/32_scMultiomeGRN_infer.slurm"
