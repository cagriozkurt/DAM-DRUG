#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=celloracle_setup
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/celloracle_setup-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/celloracle_setup-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=16G
# ─────────────────────────────────────────────────────────────────────────────
# OBSOLETE — CellOracle is now pre-installed in
#   docker://ghcr.io/cagriozkurt/dam-drug-scmultiomegrn:latest
# Do not submit this script. Kept for historical reference only.
# ─────────────────────────────────────────────────────────────────────────────
# Step 1 of 5.7 CellOracle:
#   - Install celloracle into scMultiomeGRN env
#   - Print version, available classes, and key method signatures
#   - Probe TFinfo and Oracle APIs
#   - Report for writing the inference script
#
# Results printed to .out log; paste to Claude to write the inference script.
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
mkdir -p $PROJDIR/logs

echo "=== CellOracle Setup Probe ==="
echo "Start: $(date)"

module load miniconda3

# ── Install celloracle ────────────────────────────────────────────────────────
echo ""
echo "--- Installing celloracle ---"

# Check if already installed
if conda run -n scMultiomeGRN python -c "import celloracle" 2>/dev/null; then
    echo "celloracle already installed."
else
    # Root cause: velocyto's setup.py imports Cython and numpy at metadata-generation
    # time, before pip's resolver has installed them. Must pre-install both, then
    # build velocyto with --no-build-isolation so setup.py can see the conda env.

    echo "Step 1: Pre-install Cython + numpy (velocyto build-time deps)..."
    conda run -n scMultiomeGRN pip install "cython" "numpy<2" 2>&1

    echo "Step 2: Install velocyto with --no-build-isolation..."
    conda run -n scMultiomeGRN pip install velocyto==0.17.17 --no-build-isolation 2>&1 || true

    echo "Step 3: Install celloracle (velocyto already satisfied)..."
    conda run -n scMultiomeGRN pip install celloracle==0.20.0 2>&1 || true

    # Final check
    conda run -n scMultiomeGRN python -c "import celloracle; print('celloracle', celloracle.__version__)" \
        || echo "ERROR: celloracle installation failed"
fi

# ── Probe API ─────────────────────────────────────────────────────────────────
echo ""
echo "--- CellOracle API probe ---"

cat > /tmp/celloracle_probe_${SLURM_JOB_ID}.py << 'PYEOF'
import celloracle as co
print(f"celloracle version: {co.__version__}")
print()

# Top-level public API
print("=== celloracle public attributes ===")
import inspect
members = [(n, type(v).__name__) for n, v in inspect.getmembers(co)
           if not n.startswith('_')]
for name, kind in members:
    print(f"  {name}: {kind}")

print()
print("=== Oracle class ===")
oracle = co.Oracle
print(inspect.getdoc(oracle))
print()
print("Oracle public methods:")
for name, method in inspect.getmembers(oracle, predicate=inspect.isfunction):
    if not name.startswith('_'):
        try:
            sig = str(inspect.signature(method))
        except Exception:
            sig = "(...)"
        print(f"  {name}{sig}")

print()
print("=== TFinfo class ===")
try:
    tfinfo = co.TFinfo
    print(inspect.getdoc(tfinfo))
    print()
    print("TFinfo public methods:")
    for name, method in inspect.getmembers(tfinfo, predicate=inspect.isfunction):
        if not name.startswith('_'):
            try:
                sig = str(inspect.signature(method))
            except Exception:
                sig = "(...)"
            print(f"  {name}{sig}")
except AttributeError:
    print("TFinfo not found at co.TFinfo — checking submodules...")
    import celloracle.motif_analysis as ma
    print("motif_analysis members:", dir(ma))

print()
print("=== celloracle.motif_analysis ===")
try:
    import celloracle.motif_analysis as ma
    print(dir(ma))
    if hasattr(ma, 'TFinfo'):
        tfinfo = ma.TFinfo
        print("TFinfo methods:")
        for name, method in inspect.getmembers(tfinfo, predicate=inspect.isfunction):
            if not name.startswith('_'):
                print(f"  {name}")
except Exception as e:
    print(f"  Error: {e}")

print()
print("=== Input format: what does Oracle.import_anndata_as_raw_count expect? ===")
try:
    sig = inspect.signature(co.Oracle.import_anndata_as_raw_count)
    print(f"  import_anndata_as_raw_count{sig}")
    doc = inspect.getdoc(co.Oracle.import_anndata_as_raw_count)
    if doc:
        print(doc[:500])
except Exception as e:
    print(f"  {e}")

print()
print("=== Perturbation simulation ===")
for method_name in ["simulate_shift", "simulate_count_GRN_KO",
                    "calculate_p_mass", "suggest_corrplot_n_grid"]:
    try:
        m = getattr(co.Oracle, method_name)
        sig = inspect.signature(m)
        print(f"  {method_name}{sig}")
    except Exception:
        print(f"  {method_name}: not found")

PYEOF

conda run -n scMultiomeGRN python /tmp/celloracle_probe_${SLURM_JOB_ID}.py

echo ""
echo "=== Setup complete: $(date) ==="
echo "Paste the above output to Claude to write the inference script."
