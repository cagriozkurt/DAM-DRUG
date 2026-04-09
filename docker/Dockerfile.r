# DAM-DRUG — R / CellChat environment
#
# Contains: R 4.3.3 + all packages from envs/cellchat_r.yml
#           CellChat v2.2.0 (GitHub), NicheNet (via CellChat dep)
#
# Build (from repo root):
#   docker build -f docker/Dockerfile.r -t mozkurt/dam-drug-r:latest .
#
# Push:
#   docker push mozkurt/dam-drug-r:latest
#
# Pull on TRUBA (converts to SIF automatically):
#   apptainer pull ~/containers/dam-drug-r.sif docker://mozkurt/dam-drug-r:latest
#
# Run locally:
#   docker run --rm -v $(pwd):/work mozkurt/dam-drug-r \
#     Rscript code/phase2_LR/37_cellchat_nichechat.R

FROM rocker/r-ver:4.3.3

LABEL org.opencontainers.image.description="DAM-DRUG R/CellChat environment"
LABEL org.opencontainers.image.source="https://github.com/mozkurt/DAM-DRUG"

# System libraries required by R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    libgmp-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages — versions pinned to match envs/cellchat_r.yml
RUN Rscript -e "install.packages(c( \
    'remotes', 'BiocManager', 'devtools', \
    'Matrix', 'data.table', 'dplyr', 'tidyr', 'tidyverse', \
    'ggplot2', 'patchwork', 'ggraph', 'ggalluvial', 'circlize', \
    'svglite', 'ggrepel', 'cowplot', \
    'igraph', 'NMF', \
    'future', 'future.apply', \
    'hdf5r', 'curl', 'xml2', 'httr', 'jsonlite', \
    'nloptr', 'lme4', 'pbkrtest', 'car', 'rstatix', 'ggpubr' \
    ), repos='https://cloud.r-project.org', Ncpus=4)"

# Bioconductor 3.18 — matches R 4.3.x
RUN Rscript -e "BiocManager::install( \
    c('ComplexHeatmap', 'limma', 'MAST', 'BiocNeighbors'), \
    version='3.18', update=FALSE, ask=FALSE)"

# CellChat v2 from GitHub (version 2.2.0; jinworks fork)
RUN Rscript -e "remotes::install_github( \
    'jinworks/CellChat', \
    upgrade='never', \
    build_vignettes=FALSE)"

WORKDIR /work
