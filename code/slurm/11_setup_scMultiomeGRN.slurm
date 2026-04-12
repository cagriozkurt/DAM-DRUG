#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=scMultiomeGRN_setup
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# Step 1 of 2.8 scMultiomeGRN:
#   - Download from Zenodo 10.5281/zenodo.14848389
#   - Inspect package structure and API
#   - Probe input format requirements
#   - Report for writing the inference script
#
# Results printed to .out log; paste to Claude to write the inference script.
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
WORKDIR=$PROJDIR/tools/scMultiomeGRN
mkdir -p "$WORKDIR" "$PROJDIR/logs"

echo "=== scMultiomeGRN Setup ==="
echo "Start: $(date)"

# ── Download from Zenodo ───────────────────────────────────────────────────────
echo ""
echo "--- Downloading from Zenodo 10.5281/zenodo.14848389 ---"
cd "$WORKDIR"

# Zenodo API to get download links
ZENODO_RECORD="14848389"
ZENODO_API="https://zenodo.org/api/records/${ZENODO_RECORD}"

echo "Fetching record metadata..."
curl -s "$ZENODO_API" | python3 -c "
import json, sys
rec = json.load(sys.stdin)
print('Title:', rec.get('metadata', {}).get('title', 'N/A'))
print('DOI:', rec.get('doi', 'N/A'))
print('Files:')
for f in rec.get('files', []):
    print(f'  {f[\"key\"]}  ({f[\"size\"]:,} bytes)  {f[\"links\"][\"self\"]}')
"

# Download the archive (first file found)
echo ""
echo "Downloading archive..."
curl -s "$ZENODO_API" | python3 -c "
import json, sys, subprocess
rec = json.load(sys.stdin)
files = rec.get('files', [])
if not files:
    print('ERROR: no files found in Zenodo record')
    sys.exit(1)
# Download all files
for f in files:
    url = f['links']['self']
    name = f['key']
    print(f'Downloading {name} ...')
    subprocess.run(['wget', '-q', '--show-progress', '-O', name, url])
"

echo "Download done: $(date)"
ls -lh "$WORKDIR"

# ── Inspect and extract ────────────────────────────────────────────────────────
echo ""
echo "--- Inspecting package contents ---"
for f in "$WORKDIR"/*.zip "$WORKDIR"/*.tar.gz "$WORKDIR"/*.tar; do
    [ -f "$f" ] || continue
    echo "Archive: $f"
    case "$f" in
        *.zip)     python3 -c "import zipfile; z=zipfile.ZipFile('$f'); [print(' ',n) for n in z.namelist()[:50]]" ;;
        *.tar.gz)  tar -tzf "$f" | head -50 ;;
        *.tar)     tar -tf  "$f" | head -50 ;;
    esac
done

# ── Try pip install ────────────────────────────────────────────────────────────
echo ""
echo "--- Attempting pip install ---"

# Try direct pip install from PyPI first
pip install scMultiomeGRN 2>&1 | tail -5 || true

# Try installing from downloaded package
for f in "$WORKDIR"/*.zip "$WORKDIR"/*.tar.gz "$WORKDIR"/*.whl; do
    [ -f "$f" ] || continue
    echo "Trying: pip install $f"
    pip install "$f" 2>&1 | tail -10 || true
done

# Extract and look for setup.py / pyproject.toml
mkdir -p "$WORKDIR/extracted"
for f in "$WORKDIR"/*.zip; do
    [ -f "$f" ] || continue
    python3 -c "import zipfile; z=zipfile.ZipFile('$f'); z.extractall('$WORKDIR/extracted')" 2>/dev/null || true
done
for f in "$WORKDIR"/*.tar.gz; do
    [ -f "$f" ] || continue
    tar -xzf "$f" -C "$WORKDIR/extracted/" || true
done

echo ""
echo "--- Extracted contents ---"
find "$WORKDIR/extracted" -maxdepth 3 -name "*.py" | head -30
find "$WORKDIR/extracted" -maxdepth 2 -name "setup.py" -o -name "pyproject.toml" -o -name "requirements.txt" | head -10

# ── Print README and requirements ─────────────────────────────────────────────
echo ""
echo "--- README ---"
find "$WORKDIR/extracted" -maxdepth 3 -name "README*" | head -3 | xargs -I{} sh -c 'echo "=== {} ==="; head -80 "{}"'

echo ""
echo "--- requirements.txt ---"
find "$WORKDIR/extracted" -maxdepth 3 -name "requirements*.txt" | xargs -I{} sh -c 'echo "=== {} ==="; cat "{}"'

# ── Check if importable ────────────────────────────────────────────────────────
echo ""
echo "--- Import test ---"
python3 -c "
import sys
sys.path.insert(0, '$WORKDIR/extracted')

# Try various import names
for name in ['scMultiomeGRN', 'scmultiomegrn', 'scMultiome', 'grn']:
    try:
        m = __import__(name)
        print(f'SUCCESS: import {name}')
        print(dir(m))
        break
    except ImportError as e:
        print(f'  import {name} failed: {e}')
"

echo ""
echo "=== Setup complete: $(date) ==="
echo "Paste the above output to Claude to write the inference script."
