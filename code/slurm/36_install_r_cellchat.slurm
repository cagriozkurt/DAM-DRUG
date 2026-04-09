#!/bin/bash
# OBSOLETE — CellChat is now pre-installed in docker://ghcr.io/cagriozkurt/dam-drug-r:latest
# Use: apptainer pull ~/containers/dam-drug-r.sif docker://ghcr.io/cagriozkurt/dam-drug-r:latest
# This script is kept for historical reference only; do not submit.
#SBATCH --account=mozkurt
#SBATCH --job-name=r_install
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/r_install-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/r_install-%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=64G
# ─────────────────────────────────────────────────────────────────────────────
# One-time: create conda env "cellchat_r" with R 4.3 + system libs,
# then install CellChat + NicheNet from GitHub.
#
# conda-forge handles libcurl, libxml2, libpng, harfbuzz, nlopt etc.
# — avoids the "unable to compile from source" failures seen with module R.
#
# Submit: sbatch code/slurm/36_install_r_cellchat.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
mkdir -p $PROJDIR/logs

module load miniconda3

# Isolate conda R from polluted user library (compiled against Intel SVML not present on nodes)
export R_LIBS_USER=""
export R_LIBS=""

echo "=== Creating cellchat_r conda env ==="
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

# Skip creation if env already exists
if conda env list | grep -q "^cellchat_r "; then
    echo "cellchat_r env already exists — skipping conda create"
else
    conda create -y -n cellchat_r -c conda-forge -c bioconda \
        r-base=4.3 \
        r-matrix r-data.table \
        r-dplyr r-tidyr r-ggplot2 r-patchwork r-ggraph r-ggalluvial \
        r-igraph r-nmf r-circlize \
        r-future r-future.apply \
        r-remotes r-devtools r-biocmanager \
        r-hdf5r \
        r-curl r-xml2 r-httr r-jsonlite \
        bioconductor-complexheatmap bioconductor-limma bioconductor-mast \
        2>&1 | tail -20
    echo "conda create done: $(date)"
fi

# Install CellChat's problematic deps via conda (pre-built, avoids cmake/zlib link issues)
echo ""
echo "=== Installing CellChat deps via conda ==="
conda install -y -n cellchat_r -c conda-forge -c bioconda \
    r-nloptr r-lme4 r-pbkrtest r-car r-rstatix r-ggpubr r-svglite \
    bioconductor-biocneighbors \
    2>&1 | tail -30
echo "conda deps done: $(date)"

# Install CellChat from GitHub (not on conda-forge)
echo ""
echo "=== Installing CellChat from GitHub ==="
conda run -n cellchat_r Rscript -e "
options(Ncpus=28L)
lib <- .libPaths()[grep('miniconda', .libPaths())][1]
cat('Installing into:', lib, '\n')

# conda-forge pins r-nmf at 0.21; CellChat requires >= 0.23 — upgrade from CRAN
if (packageVersion('NMF') < '0.23') {
    cat('Upgrading NMF from', as.character(packageVersion('NMF')), 'to CRAN version...\n')
    install.packages('NMF', lib=lib, repos='https://cloud.r-project.org', quiet=FALSE)
}

if (!requireNamespace('CellChat', quietly=TRUE)) {
    remotes::install_github('jinworks/CellChat',
        lib=lib, upgrade='never',
        build_vignettes=FALSE, quiet=FALSE)
    cat('CellChat installed\n')
} else {
    cat('CellChat already installed\n')
}
"

# NicheNet skipped: requires Seurat/tidyverse/caret/mlrMBO — not installable from
# source without full system dev libs. Use CellChat-only L-R analysis.
echo ""
echo "=== NicheNet skipped (dep tree requires Seurat — not feasible from source) ==="

# Verify
echo ""
echo "=== Verifying installs ==="
conda run -n cellchat_r Rscript -e "
cat('libPaths:\n'); cat(paste(' ', .libPaths()), sep='\n'); cat('\n')
for (pkg in c('CellChat', 'Matrix', 'ComplexHeatmap', 'hdf5r', 'igraph')) {
    ok <- requireNamespace(pkg, quietly=TRUE)
    cat(sprintf('  %-20s %s\n', pkg, if(ok) 'OK' else 'MISSING'))
}
"

echo ""
echo "=== Done: $(date) ==="
