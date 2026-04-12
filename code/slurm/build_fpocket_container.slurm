#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=build-fpocket
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G

# Build a new SIF that extends scenic.sif with fpocket, biopython,
# pdbfixer, and openmm.  Strategy: copy sandbox, install into it, rebuild SIF.

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
mkdir -p "$PROJDIR/logs" "$PROJDIR/containers"

SANDBOX=$PROJDIR/containers/fpocket-env
SIF=$PROJDIR/containers/fpocket-env.sif
BASE_SIF=$PROJDIR/containers/scenic.sif

echo "Build node: $SLURMD_NODENAME  $(date)"

# Step 1: create sandbox from existing scenic.sif
echo "=== Creating sandbox from $BASE_SIF ==="
apptainer build --sandbox --force "$SANDBOX" "$BASE_SIF"

# Step 1b: create a bind-point directory inside the sandbox so that
# cluster auto-bind mounts (e.g. /arf on TRUBA) do not cause FATAL errors.
echo "=== Creating bind point inside sandbox ==="
mkdir -p "$SANDBOX/arf"

# Step 2a: install Python packages via pip (avoids conda solver deadlock on openmm)
echo "=== Installing biopython + pdbfixer + openmm via pip ==="
apptainer exec --writable "$SANDBOX" \
    pip install --no-cache-dir biopython pdbfixer openmm

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

# Step 3: verify (--no-mount bind prevents auto-bind mounts from interfering)
echo "=== Verification ==="
apptainer exec --no-mount bind "$SANDBOX" \
    fpocket --version 2>&1 | head -2
apptainer exec --no-mount bind "$SANDBOX" \
    python -c "import Bio; print('biopython', Bio.__version__)"
apptainer exec --no-mount bind "$SANDBOX" \
    python -c "import pdbfixer; print('pdbfixer OK')"
apptainer exec --no-mount bind "$SANDBOX" \
    python -c "import openmm; print('openmm', openmm.__version__)"

# Step 4: build SIF
echo "=== Building SIF: $SIF ==="
apptainer build --force "$SIF" "$SANDBOX"
echo "Build exit code: $?"
ls -lh "$SIF"

echo "Done: $(date)"
