#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=tier2_prep
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 4 contingency — Tier 2 PDBQT generation (batch step)
#
# PREREQUISITE: run 26_fetch_tier2.py interactively on the login node first:
#   apptainer exec containers/scenic.sif python code/slurm/26_fetch_tier2.py
# That script fetches from ChEMBL (needs internet) and writes
#   data/compounds/tier2_approved.csv
# This SLURM script then generates 3D conformers + PDBQTs in parallel
# (no internet required — pure compute).
#
# Output:
#   data/compounds/tier2_pdbqt/*.pdbqt
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
OUT_DIR=$PROJDIR/data/compounds
PDBQT_DIR=$OUT_DIR/tier2_pdbqt
mkdir -p "$PDBQT_DIR" "$PROJDIR/logs"

SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

echo "Start: $(date)"
echo "Node:  $SLURMD_NODENAME"

# Verify prerequisite CSV exists (written by 26_fetch_tier2.py on login node)
TIER2_CSV=$OUT_DIR/tier2_approved.csv
[ -f "$TIER2_CSV" ] || { echo "ERROR: $TIER2_CSV not found — run 26_fetch_tier2.py on login node first"; exit 1; }
echo "Tier-2 CSV: $(wc -l < "$TIER2_CSV") rows"

echo ""
echo "=== 3D conformers + PDBQT ==="

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python - << PYEOF
import subprocess, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import os as _os
import pandas as pd

_PROJ = _os.environ.get("DAM_DRUG_DIR", str(Path.cwd()))
TIER2_CSV = Path(_PROJ) / 'data/compounds/tier2_approved.csv'
PDBQT_DIR = Path(_PROJ) / 'data/compounds/tier2_pdbqt'
NCPUS = 20

df = pd.read_csv(TIER2_CSV)
print(f"Generating PDBQTs for {len(df)} compounds on {NCPUS} workers...")

def make_pdbqt(row):
    cid   = row['molecule_chembl_id']
    smi   = row['smiles']
    out_p = PDBQT_DIR / f'{cid}.pdbqt'

    if out_p.exists() and out_p.stat().st_size > 100:
        return cid, 'skip'

    sdf_p = PDBQT_DIR / f'{cid}.sdf'

    # 3D conformer
    try:
        r1 = subprocess.run(
            ['obabel', f'-:{smi}', '--gen3d', '-O', str(sdf_p), '-p', '7.4', '--title', cid],
            capture_output=True, text=True, timeout=90
        )
    except subprocess.TimeoutExpired:
        sdf_p.unlink(missing_ok=True)
        return cid, 'obabel_timeout'
    except Exception as e:
        return cid, f'obabel_error: {e}'

    if not sdf_p.exists() or sdf_p.stat().st_size < 20:
        return cid, f'obabel_fail: {r1.stderr.strip()[:80]}'

    # PDBQT via meeko
    try:
        r2 = subprocess.run(
            ['mk_prepare_ligand.py', '-i', str(sdf_p), '-o', str(out_p)],
            capture_output=True, text=True, timeout=60
        )
    except subprocess.TimeoutExpired:
        sdf_p.unlink(missing_ok=True)
        return cid, 'meeko_timeout'
    except Exception as e:
        return cid, f'meeko_error: {e}'
    finally:
        sdf_p.unlink(missing_ok=True)

    if not out_p.exists() or out_p.stat().st_size < 50:
        return cid, f'meeko_fail: {r2.stderr.strip()[:80]}'

    return cid, 'ok'

results = {'ok': 0, 'skip': 0, 'fail': 0}
fail_log = []

with ProcessPoolExecutor(max_workers=NCPUS) as ex:
    futs = {ex.submit(make_pdbqt, row): row['molecule_chembl_id']
            for _, row in df.iterrows()}
    for i, fut in enumerate(as_completed(futs), 1):
        try:
            cid, status = fut.result()
        except Exception as e:
            cid = futs[fut]
            status = f'future_error: {e}'
        if status == 'ok':
            results['ok'] += 1
        elif status == 'skip':
            results['skip'] += 1
        else:
            results['fail'] += 1
            fail_log.append(f'{cid}: {status}')
        if i % 200 == 0:
            print(f"  {i}/{len(df)}  ok={results['ok']} skip={results['skip']} fail={results['fail']}")

print(f"\nFinal: ok={results['ok']}  skip={results['skip']}  fail={results['fail']}")
if fail_log:
    (Path(_PROJ) / 'data/compounds/tier2_failures.txt').write_text('\n'.join(fail_log[:50]))
    print(f"First failures written to tier2_failures.txt")

# Count valid PDBQTs
n_ready = len(list(PDBQT_DIR.glob('*.pdbqt')))
print(f"\nReady PDBQTs: {n_ready}")
PYEOF

echo ""
echo "Finished: $(date)"
echo "Tier-2 library at: $OUT_DIR/tier2_approved.csv + tier2_pdbqt/"
