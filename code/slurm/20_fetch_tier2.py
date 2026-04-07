"""
Run interactively on the TRUBA login node (needs internet):
    conda activate mmpbsa
    pip install requests chembl_webresource_client
    python3 20_fetch_tier2.py

Fetches all ChEMBL max_phase=4 compounds, applies Ro5 + CNS MPO filters
locally, writes data/compounds/tier2_approved.csv (~1500-2500 rows, ~5-15 min).
PDBQT generation is handled by 20_prep_tier2_library.slurm (batch job).
"""
import sys, time
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas not found")

_TRUBA = Path('/arf/scratch/mozkurt/DAM-DRUG/data/compounds')
_LOCAL = Path(__file__).resolve().parents[2] / 'data/compounds'
OUT_DIR   = _TRUBA if _TRUBA.exists() else _LOCAL
TIER1_CSV = OUT_DIR / 'tier1_cns_approved.csv'
OUT_CSV   = OUT_DIR / 'tier2_approved.csv'
OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"Output dir: {OUT_DIR}")

# ── Load Tier-1 IDs to exclude ───────────────────────────────────────────────
tier1_ids = set()
if TIER1_CSV.exists():
    t1 = pd.read_csv(TIER1_CSV)
    chembl_col = next((c for c in t1.columns if 'chembl' in c.lower()), None)
    if chembl_col:
        tier1_ids = set(t1[chembl_col].astype(str))
print(f"Tier-1 IDs to exclude: {len(tier1_ids)}")

# ── Fetch via chembl_webresource_client (handles pagination + retries) ────────
try:
    from chembl_webresource_client.new_client import new_client
    USE_CLIENT = True
    print("Using chembl_webresource_client")
except ImportError:
    USE_CLIENT = False
    print("chembl_webresource_client not found — falling back to requests")

def fetch_with_client():
    molecule = new_client.molecule
    # Simple query: max_phase=4 only; property filtering done locally
    # This avoids the 500 errors from complex combined server-side filters
    results = molecule.filter(max_phase=4).only([
        'molecule_chembl_id', 'pref_name',
        'molecule_structures', 'molecule_properties', 'max_phase',
    ])
    mols = []
    for i, mol in enumerate(results):
        mols.append(mol)
        if (i + 1) % 500 == 0:
            print(f"  Fetched {i+1} molecules...", flush=True)
    return mols

def fetch_with_requests():
    import requests
    BASE = 'https://www.ebi.ac.uk/chembl/api/data/molecule'
    # Minimal query — only max_phase filter to avoid 500 errors
    params = {'max_phase': 4, 'format': 'json', 'limit': 500, 'offset': 0}
    all_mols = []
    page = 0
    while True:
        for attempt in range(5):
            try:
                r = requests.get(BASE, params=params, timeout=90)
                r.raise_for_status()
                break
            except Exception as e:
                wait = 10 * (attempt + 1)
                print(f"  Page {page} attempt {attempt+1} failed: {e} — retrying in {wait}s")
                time.sleep(wait)
        else:
            print(f"  Page {page} failed after 5 attempts — stopping early")
            break
        data = r.json()
        mols = data.get('molecules', [])
        all_mols.extend(mols)
        print(f"  Page {page:3d}: +{len(mols):4d}  (total: {len(all_mols)})", flush=True)
        if not data.get('page_meta', {}).get('next'):
            break
        params['offset'] += 500
        page += 1
        time.sleep(0.5)  # be polite to the server
    return all_mols

print("\nFetching from ChEMBL (max_phase=4, all areas)...")
if USE_CLIENT:
    all_mols = fetch_with_client()
else:
    all_mols = fetch_with_requests()

print(f"\nTotal fetched: {len(all_mols)}")

# ── Local filtering (no RDKit needed — use ChEMBL-provided properties) ────────
# Ro5: MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10
# Simplified CNS MPO proxy (no pKa): MW ≤ 500, logP ≤ 5, TPSA ≤ 120, HBD ≤ 3
rows = []
n_excl = n_nosmi = n_ro5 = n_mpo = 0

for mol in all_mols:
    cid = (mol.get('molecule_chembl_id') or '')

    if cid in tier1_ids:
        n_excl += 1
        continue

    structs = mol.get('molecule_structures') or {}
    smi = (structs.get('canonical_smiles') or '').strip()
    if not smi:
        n_nosmi += 1
        continue

    props = mol.get('molecule_properties') or {}
    try:
        mw   = float(props.get('mw_freebase') or 999)
        logp = float(props.get('alogp') or 999)
        hbd  = int(float(props.get('hbd') or 99))
        hba  = int(float(props.get('hba') or 99))
        tpsa = float(props.get('psa') or 999)
        rtb  = int(float(props.get('rtb') or 99))
        ro5v = int(float(props.get('num_ro5_violations') or 0))
    except (TypeError, ValueError):
        n_ro5 += 1
        continue

    # Ro5
    if ro5v > 0 or mw > 500 or logp > 5 or hbd > 5 or hba > 10:
        n_ro5 += 1
        continue

    # Simplified CNS MPO proxy: needs to pass at least 3 of 5 checks
    cns = 0
    cns += 1 if mw   <= 360 else (0.5 if mw   <= 500 else 0)
    cns += 1 if logp <= 3.0 else (0.5 if logp <= 5.0 else 0)
    cns += 1 if 40 <= tpsa <= 90 else (0.5 if tpsa <= 120 else 0)
    cns += 1 if hbd  == 0  else (0.5 if hbd  <= 3  else 0)
    cns += 1 if rtb  <= 8  else (0.5 if rtb  <= 10 else 0)
    if cns < 3.0:
        n_mpo += 1
        continue

    rows.append({
        'molecule_chembl_id': cid,
        'smiles':             smi,
        'pref_name':          mol.get('pref_name') or '',
        'mw':    mw,
        'alogp': logp,
        'hbd':   hbd,
        'hba':   hba,
        'tpsa':  tpsa,
        'cns_mpo': cns,
    })

print(f"\nFilter summary:")
print(f"  Tier-1 excluded  : {n_excl}")
print(f"  No SMILES        : {n_nosmi}")
print(f"  Ro5 fail         : {n_ro5}")
print(f"  CNS MPO < 3.0    : {n_mpo}")
print(f"  Passed           : {len(rows)}")

df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)
print(f"\nWritten: {OUT_CSV}")
print(f"Next step: sbatch 20_prep_tier2_library.slurm")
