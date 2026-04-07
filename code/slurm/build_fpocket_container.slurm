#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=build-fpocket
#SBATCH --output=/arf/scratch/mozkurt/build-fpocket-%j.out
#SBATCH --error=/arf/scratch/mozkurt/build-fpocket-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G

# Build a new SIF that extends scanpy-env with fpocket, biopython,
# pdbfixer, and openmm.  Strategy: copy sandbox, install into it, rebuild SIF.

echo "Build node: $SLURMD_NODENAME  $(date)"

SANDBOX=$HOME/containers/fpocket-env
SIF=$HOME/containers/fpocket-env.sif
BASE_SIF=$HOME/containers/scanpy-env.sif

# Step 1: create sandbox from existing scanpy-env.sif
echo "=== Creating sandbox from $BASE_SIF ==="
apptainer build --sandbox --force "$SANDBOX" "$BASE_SIF"

# Step 1b: create /arf mount point inside sandbox
# TRUBA apptainer auto-binds /arf on compute nodes; without this directory
# inside the container every subsequent --writable exec fails with FATAL mount error.
echo "=== Creating /arf bind point inside sandbox ==="
mkdir -p "$SANDBOX/arf"

# Step 2a: install Python packages via pip (avoids conda solver deadlock on openmm)
echo "=== Installing biopython + pdbfixer + openmm via pip ==="
apptainer exec --writable "$SANDBOX" \
    /opt/conda/bin/pip install --no-cache-dir biopython pdbfixer openmm

# Step 2b: compile fpocket from source (skip bioconda — conda solver deadlocks in py3.12 env)
echo "=== Compiling fpocket from source ==="
apptainer exec --writable "$SANDBOX" bash -c "
    set -e
    cd /tmp
    git clone --depth=1 https://github.com/Discngine/fpocket.git
    cd fpocket
    make
    cp bin/fpocket bin/dpocket bin/mdpocket /opt/conda/bin/ 2>/dev/null || cp bin/fpocket /opt/conda/bin/
    echo 'fpocket compiled OK'
"

# Step 3: verify (use --no-mount bind to avoid /arf auto-bind interfering with verification)
echo "=== Verification ==="
apptainer exec --no-mount bind "$SANDBOX" \
    /opt/conda/bin/conda run -n base fpocket --version 2>&1 | head -2
apptainer exec --no-mount bind "$SANDBOX" \
    /opt/conda/bin/conda run -n base python -c "import Bio; print('biopython', Bio.__version__)"
apptainer exec --no-mount bind "$SANDBOX" \
    /opt/conda/bin/conda run -n base python -c "import pdbfixer; print('pdbfixer OK')"
apptainer exec --no-mount bind "$SANDBOX" \
    /opt/conda/bin/conda run -n base python -c "import openmm; print('openmm', openmm.__version__)"

# Step 4: build SIF
echo "=== Building SIF: $SIF ==="
apptainer build --force "$SIF" "$SANDBOX"
echo "Build exit code: $?"
ls -lh "$SIF"

echo "Done: $(date)"
