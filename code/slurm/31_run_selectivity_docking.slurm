#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=selectivity_dock
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 4 Selectivity docking — 6 MM-GBSA hits vs DRD2, HTR2A, hERG
#
# Compounds (validated by MM-GBSA, ΔG < −4 kcal/mol):
#   CHEMBL118      CELECOXIB     IKZF1_8RQC  −30.5 kcal/mol
#   CHEMBL973      TERIFLUNOMIDE IKZF1_8RQC  −15.2 kcal/mol
#   CHEMBL1096     AMLEXANOX     IKZF1_8RQC  −14.2 kcal/mol
#   CHEMBL3188267  CAPMATINIB    IRF8_AF2_DBD −14.0 kcal/mol
#   CHEMBL3039520  LASMIDITAN    IRF8_AF2_DBD −10.3 kcal/mol
#   CHEMBL4297528  RISDIPLAM     IKZF1_8RQC   −6.5 kcal/mol
#
# Off-targets: DRD2 (6CM4), HTR2A (7JP2), HERG (7CN1)
# Total docking tasks: 18  (6 compounds × 3 off-targets)
#
# Requires:
#   08_prep_selectivity_receptors.py run first (receptors + configs ready)
#   data/compounds/tier2_pdbqt/{CHEMBL_ID}.pdbqt   (already prepared)
#
# Selectivity index (SI):
#   SI = vina_ontarget / vina_offtarget  (both negative → SI > 1 = selective)
#   SI > 2.0 considered highly selective (target binds ≥2× more strongly)
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
OUT_DIR=$PROJDIR/results/phase4/vina_selectivity
mkdir -p "$OUT_DIR" "$PROJDIR/logs"

SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

# Override VINA by setting it in your environment before submitting, or place
# the binary at $PROJDIR/bin/vina.
export VINA=${VINA:-$PROJDIR/bin/vina}
[ -f "$VINA" ] || { echo "ERROR: Vina not found at $VINA — set VINA env var or place binary at \$PROJDIR/bin/vina"; exit 1; }
echo "Vina: $VINA"

echo "Start: $(date)"
echo "Node:  $SLURMD_NODENAME"

# ── Verify receptors exist ────────────────────────────────────────────────────
for name in DRD2 HTR2A HERG; do
    f="$PROJDIR/data/docking/receptors/${name}_prep.pdbqt"
    [ -f "$f" ] || { echo "ERROR: $f not found — run 08_prep_selectivity_receptors.py first"; exit 1; }
done

# ── Step 1: Vina docking ─────────────────────────────────────────────────────
echo ""
echo "=== Step 1: AutoDock Vina (6 hits × 3 off-targets, 18 workers) ==="

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python - << 'PYEOF'
import csv
import os
import shutil
import subprocess
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

PROJDIR   = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
PDBQT_DIR = PROJDIR / 'data/compounds/tier2_pdbqt'
OUT_DIR   = PROJDIR / 'results/phase4/vina_selectivity'

# 6 validated hits from MM-GBSA
HITS = [
    ("CHEMBL118",     "CELECOXIB",     "IKZF1_8RQC"),
    ("CHEMBL973",     "TERIFLUNOMIDE", "IKZF1_8RQC"),
    ("CHEMBL1096",    "AMLEXANOX",     "IKZF1_8RQC"),
    ("CHEMBL3188267", "CAPMATINIB",    "IRF8_AF2_DBD"),
    ("CHEMBL3039520", "LASMIDITAN",    "IRF8_AF2_DBD"),
    ("CHEMBL4297528", "RISDIPLAM",     "IKZF1_8RQC"),
]

OFF_TARGETS = {
    'DRD2':  PROJDIR / 'data/docking/configs/DRD2.conf',
    'HTR2A': PROJDIR / 'data/docking/configs/HTR2A.conf',
    'HERG':  PROJDIR / 'data/docking/configs/HERG.conf',
}

VINA = os.environ.get("VINA") or shutil.which('vina') or shutil.which('vina_1.2.3_linux_x86_64')

def parse_vina_score(text):
    for line in text.splitlines():
        s = line.strip()
        if s and s[0].isdigit():
            parts = s.split()
            if len(parts) >= 2:
                try: return float(parts[1])
                except: pass
        if 'VINA RESULT' in s:
            parts = s.split()
            if len(parts) >= 4:
                try: return float(parts[3])
                except: pass
    return None

def dock_one(args):
    off_name, conf_path, chembl_id = args
    lig_path = PDBQT_DIR / f'{chembl_id}.pdbqt'
    if not lig_path.exists():
        return chembl_id, off_name, None, 'missing_lig'

    out_p = OUT_DIR / f'{off_name}_{chembl_id}_out.pdbqt'
    if out_p.exists() and out_p.stat().st_size > 50:
        score = parse_vina_score(out_p.read_text())
        return chembl_id, off_name, score, 'skip'

    cmd = [VINA, '--config', str(conf_path),
           '--ligand', str(lig_path), '--out', str(out_p), '--cpu', '1']
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return chembl_id, off_name, None, 'timeout'

    if out_p.exists() and out_p.stat().st_size > 50:
        score = parse_vina_score(out_p.read_text())
    else:
        score = parse_vina_score(r.stdout) or parse_vina_score(r.stderr)
    return chembl_id, off_name, score, 'ok'

# Verify configs exist
for off_name, conf in OFF_TARGETS.items():
    if not conf.exists():
        raise FileNotFoundError(f"Config not found: {conf}")

# Build 18 tasks
tasks = []
for chembl_id, _, _ in HITS:
    for off_name, conf in OFF_TARGETS.items():
        tasks.append((off_name, conf, chembl_id))

print(f"Docking tasks: {len(tasks)}")

results = {}
with ProcessPoolExecutor(max_workers=20) as ex:
    futs = {ex.submit(dock_one, task): task for task in tasks}
    for fut in as_completed(futs):
        chembl_id, off_name, score, status = fut.result()
        results[(chembl_id, off_name)] = score
        print(f"  {off_name:6s} {chembl_id:14s}  Vina={score}  [{status}]")

# Write raw scores
with open(OUT_DIR / 'selectivity_vina_raw.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['chembl_id', 'off_target', 'vina_score'])
    for (cid, off), score in sorted(results.items()):
        w.writerow([cid, off, score])

print(f"\nRaw scores written: {OUT_DIR}/selectivity_vina_raw.csv")
PYEOF

# ── Step 2: Compute Selectivity Index ────────────────────────────────────────
echo ""
echo "=== Step 2: Selectivity Index ==="

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python - << 'PYEOF'
import csv
import os
from pathlib import Path

PROJDIR = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
RAW_CSV = PROJDIR / 'results/phase4/vina_selectivity/selectivity_vina_raw.csv'
OUT_CSV = PROJDIR / 'results/phase4/vina_selectivity/selectivity_table.csv'

# MM-GBSA-validated on-target Vina scores (from tier-2 docking)
ON_TARGET_VINA = {
    "CHEMBL118":     ("IKZF1_8RQC",   -6.286),
    "CHEMBL973":     ("IKZF1_8RQC",   -6.029),
    "CHEMBL1096":    ("IKZF1_8RQC",   -6.334),
    "CHEMBL3188267": ("IRF8_AF2_DBD", -7.230),
    "CHEMBL3039520": ("IRF8_AF2_DBD", -7.118),
    "CHEMBL4297528": ("IKZF1_8RQC",   -7.280),
}

PREF_NAME = {
    "CHEMBL118":     "CELECOXIB",
    "CHEMBL973":     "TERIFLUNOMIDE",
    "CHEMBL1096":    "AMLEXANOX",
    "CHEMBL3188267": "CAPMATINIB",
    "CHEMBL3039520": "LASMIDITAN",
    "CHEMBL4297528": "RISDIPLAM",
}

MMPBSA_DG = {
    "CHEMBL118":     -30.5,
    "CHEMBL973":     -15.24,
    "CHEMBL1096":    -14.18,
    "CHEMBL3188267": -14.03,
    "CHEMBL3039520": -10.33,
    "CHEMBL4297528":  -6.49,
}

OFF_TARGETS = ["DRD2", "HTR2A", "HERG"]

# Load off-target scores
off_scores = {}
with open(RAW_CSV) as f:
    for row in csv.DictReader(f):
        score = float(row['vina_score']) if row['vina_score'] else None
        off_scores[(row['chembl_id'], row['off_target'])] = score

# Compute SI and write table
rows = []
print(f"\n{'Compound':<16} {'On-target':<14} {'Vina':>6}  "
      f"{'DRD2':>7} {'SI_D2':>6}  "
      f"{'HTR2A':>7} {'SI_2A':>6}  "
      f"{'HERG':>7} {'SI_hERG':>7}  "
      f"{'Min_SI':>7}  dG(MM-GBSA)")
print("-" * 100)

for chembl_id, (on_target, on_vina) in ON_TARGET_VINA.items():
    name   = PREF_NAME[chembl_id]
    mmgbsa = MMPBSA_DG[chembl_id]

    si_vals = {}
    for off in OFF_TARGETS:
        off_v = off_scores.get((chembl_id, off))
        if off_v is not None and on_vina != 0:
            si = on_vina / off_v
        else:
            si = None
        si_vals[off] = (off_v, si)

    valid_si = [v[1] for v in si_vals.values() if v[1] is not None]
    min_si = min(valid_si) if valid_si else None
    flag = "PASS" if (min_si is not None and min_si > 2.0) else ("~" if (min_si and min_si > 1.0) else "FAIL")

    row = {
        'chembl_id': chembl_id,
        'pref_name': name,
        'on_target': on_target,
        'vina_ontarget': on_vina,
        'mmpbsa_dg': mmgbsa,
    }
    for off in OFF_TARGETS:
        off_v, si = si_vals[off]
        row[f'vina_{off.lower()}'] = off_v
        row[f'si_{off.lower()}'] = round(si, 2) if si else None
    row['min_si'] = round(min_si, 2) if min_si else None
    row['selective'] = flag
    rows.append(row)

    d2_v,  si_d2  = si_vals['DRD2']
    h2a_v, si_h2a = si_vals['HTR2A']
    herg_v,si_herg = si_vals['HERG']

    print(f"{name:<16} {on_target:<14} {on_vina:>6.3f}  "
          f"{d2_v or 'N/A':>7} {(f'{si_d2:.2f}' if si_d2 else 'N/A'):>6}  "
          f"{h2a_v or 'N/A':>7} {(f'{si_h2a:.2f}' if si_h2a else 'N/A'):>6}  "
          f"{herg_v or 'N/A':>7} {(f'{si_herg:.2f}' if si_herg else 'N/A'):>7}  "
          f"{(f'{min_si:.2f}' if min_si else 'N/A'):>7}  "
          f"{mmgbsa:>6.1f} kcal/mol  {flag}")

# Write CSV
fieldnames = ['chembl_id', 'pref_name', 'on_target', 'vina_ontarget', 'mmpbsa_dg',
              'vina_drd2', 'si_drd2', 'vina_htr2a', 'si_htr2a',
              'vina_herg', 'si_herg', 'min_si', 'selective']
with open(OUT_CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print(f"\nSelectivity table -> {OUT_CSV}")
print("\nNote: SI = vina_ontarget / vina_offtarget (both negative)")
print("      SI > 2.0 = target binds >=2x more strongly than off-target (PASS)")
print("      SI > 1.0 = some selectivity (~)")
print("      SI < 1.0 = off-target binds more strongly (FAIL)")
PYEOF

echo ""
echo "Finished: $(date)"
echo ""
echo "Results:"
echo "  Raw Vina scores:   $PROJDIR/results/phase4/vina_selectivity/selectivity_vina_raw.csv"
echo "  Selectivity table: $PROJDIR/results/phase4/vina_selectivity/selectivity_table.csv"
