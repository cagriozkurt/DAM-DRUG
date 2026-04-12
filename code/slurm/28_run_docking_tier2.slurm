#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=dock_tier2
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 4 Tier-2 docking — IRF8 + IKZF1 only
#
# Targets restricted to the two with convergent evidence:
#   IRF8_AF2_DBD    — DAM driver; padj=0.011 in GSE95587; core fold stable
#   IKZF1_8RQC_CRBN — strongest bulk signal (padj=0.004); molecular glue site
#
# Steps:
#   1. AutoDock Vina — all Tier-2 PDBQTs vs IRF8 + IKZF1 (20 workers, 300s/cpd)
#   2. Gnina CNN rescore — top 50 per target by Vina score
#   3. Consensus filter — Vina ≤ threshold AND CNN ≥ 0.70
#   4. Write tier2_consensus_hits.csv + tier2_shortlist.csv
#
# Requires:
#   data/compounds/tier2_pdbqt/*.pdbqt      (from 27_prep_tier2_library.slurm)
#   data/docking/receptors/IRF8_AF2_DBD_prep.pdbqt
#   data/docking/receptors/IKZF1_8RQC_CRBN_prep.pdbqt
#   data/docking/configs/IRF8_AF2_DBD.conf
#   data/docking/configs/IKZF1_8RQC.conf
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
PDBQT_DIR=$PROJDIR/data/compounds/tier2_pdbqt
OUT_DIR=$PROJDIR/results/phase4/vina_tier2
GNINA_OUT=$PROJDIR/results/phase4/gnina_tier2
mkdir -p "$OUT_DIR" "$GNINA_OUT" "$PROJDIR/logs"

SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

# Override VINA/GNINA by setting them in your environment before submitting,
# or place the binaries at $PROJDIR/bin/vina and $PROJDIR/bin/gnina.
export VINA=${VINA:-$PROJDIR/bin/vina}
export GNINA=${GNINA:-$PROJDIR/bin/gnina}

[ -f "$VINA" ]  || { echo "ERROR: Vina not found at $VINA — set VINA env var or place binary at \$PROJDIR/bin/vina"; exit 1; }
[ -f "$GNINA" ] || { echo "ERROR: gnina not found at $GNINA — set GNINA env var or place binary at \$PROJDIR/bin/gnina"; exit 1; }
echo "Vina:  $VINA"
echo "Gnina: $GNINA"

echo "Start: $(date)"
echo "Node:  $SLURMD_NODENAME"

N_PDBQT=$(ls "$PDBQT_DIR"/*.pdbqt 2>/dev/null | wc -l)
echo "Tier-2 PDBQTs: $N_PDBQT"
[ "$N_PDBQT" -eq 0 ] && { echo "ERROR: no PDBQTs in $PDBQT_DIR — run 27_prep_tier2_library.slurm first"; exit 1; }

# ── Step 1: Vina docking ─────────────────────────────────────────────────────
echo ""
echo "=== Step 1: AutoDock Vina (IRF8 + IKZF1, $N_PDBQT compounds, 20 workers) ==="

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python - << 'PYEOF'
import subprocess, os, time, shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

PROJDIR   = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
PDBQT_DIR = PROJDIR / 'data/compounds/tier2_pdbqt'
OUT_DIR   = PROJDIR / 'results/phase4/vina_tier2'
NCPUS     = 20

TARGETS = {
    'IRF8_AF2_DBD': {
        'receptor': PROJDIR / 'data/docking/receptors/IRF8_AF2_DBD_prep.pdbqt',
        'config':   PROJDIR / 'data/docking/configs/IRF8_AF2_DBD_prep.conf',
        'vina_thr': -6.5,
    },
    'IKZF1_8RQC': {
        'receptor': PROJDIR / 'data/docking/receptors/IKZF1_8RQC_CRBN_prep.pdbqt',
        'config':   PROJDIR / 'data/docking/configs/IKZF1_8RQC_CRBN.conf',
        'vina_thr': -6.0,
    },
}

# Verify receptor + config files exist
for tname, t in TARGETS.items():
    for key in ('receptor', 'config'):
        if not t[key].exists():
            raise FileNotFoundError(f"{tname} {key} not found: {t[key]}")

VINA = os.environ.get("VINA") or shutil.which('vina') or shutil.which('vina_1.2.3_linux_x86_64')

ligands = sorted(PDBQT_DIR.glob('*.pdbqt'))
print(f"Ligands: {len(ligands)}")

def parse_vina_score(text):
    """Parse best Vina score from stdout (v1.2.x) or PDBQT REMARK lines."""
    for line in text.splitlines():
        s = line.strip()
        # stdout table: '   1       -7.8  ...'
        if s and s[0].isdigit():
            parts = s.split()
            if len(parts) >= 2:
                try: return float(parts[1])
                except: pass
        # PDBQT REMARK: 'REMARK VINA RESULT:    -7.8 ...'
        if 'VINA RESULT' in s:
            parts = s.split()
            if len(parts) >= 4:
                try: return float(parts[3])
                except: pass
    return None

def dock_one(args):
    target_name, lig_path, receptor, config = args
    cid   = lig_path.stem
    out_p = OUT_DIR / f'{target_name}_{cid}_out.pdbqt'

    if out_p.exists() and out_p.stat().st_size > 50:
        score = parse_vina_score(out_p.read_text())
        return cid, target_name, score, 'skip'

    cmd = [
        VINA,
        '--config',  str(config),   # conf has receptor, box, exhaustiveness, num_modes
        '--ligand',  str(lig_path),
        '--out',     str(out_p),
        '--cpu', '1',
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return cid, target_name, None, 'timeout'

    # Parse score from output PDBQT first (reliable); fall back to stdout/stderr
    if out_p.exists() and out_p.stat().st_size > 50:
        score = parse_vina_score(out_p.read_text())
    else:
        score = parse_vina_score(r.stdout) or parse_vina_score(r.stderr)
    return cid, target_name, score, 'ok'

# Build task list
tasks = []
for tname, t in TARGETS.items():
    for lig in ligands:
        tasks.append((tname, lig, t['receptor'], t['config']))

print(f"Total docking tasks: {len(tasks)}")

all_results = []
t0 = time.time()
with ProcessPoolExecutor(max_workers=NCPUS) as ex:
    futs = {ex.submit(dock_one, task): task for task in tasks}
    for i, fut in enumerate(as_completed(futs), 1):
        cid, tname, score, status = fut.result()
        all_results.append({'compound': cid, 'target': tname,
                            'vina_score': score, 'status': status})
        if i % 500 == 0:
            elapsed = time.time() - t0
            rate = i / elapsed * 60
            print(f"  {i}/{len(tasks)} done  ({rate:.0f} tasks/min)  "
                  f"elapsed {elapsed/60:.1f} min")

df = pd.DataFrame(all_results)
df['vina_score'] = pd.to_numeric(df['vina_score'], errors='coerce')
df.to_csv(OUT_DIR / 'vina_tier2_all.csv', index=False)

# Summary per target
for tname in TARGETS:
    sub = df[df['target'] == tname].dropna(subset=['vina_score'])
    thr = TARGETS[tname]['vina_thr']
    n_pass = (sub['vina_score'] <= thr).sum()
    print(f"\n{tname}: {len(sub)} scored, {n_pass} with Vina <= {thr}")
    if len(sub) > 0:
        print(sub.nsmallest(10, 'vina_score')[['compound','vina_score']].to_string(index=False))

print(f"\nVina complete. Results: {OUT_DIR}/vina_tier2_all.csv")
PYEOF

# ── Step 2: Gnina CNN rescore (top 50 per target) ────────────────────────────
echo ""
echo "=== Step 2: Gnina CNN rescore (top 50 per target) ==="

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python - << 'PYEOF'
import subprocess, os, shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

PROJDIR   = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
PDBQT_DIR = PROJDIR / 'data/compounds/tier2_pdbqt'
VINA_DIR  = PROJDIR / 'results/phase4/vina_tier2'
GNINA_OUT = PROJDIR / 'results/phase4/gnina_tier2'
NCPUS     = 20

GNINA_BIN = os.environ.get("GNINA") or shutil.which('gnina')

TARGETS = {
    'IRF8_AF2_DBD': {
        'receptor': str(PROJDIR / 'data/docking/receptors/IRF8_AF2_DBD_prep.pdbqt'),
        'config':   PROJDIR / 'data/docking/configs/IRF8_AF2_DBD_prep.conf',
    },
    'IKZF1_8RQC': {
        'receptor': str(PROJDIR / 'data/docking/receptors/IKZF1_8RQC_CRBN_prep.pdbqt'),
        'config':   PROJDIR / 'data/docking/configs/IKZF1_8RQC_CRBN.conf',
    },
}

# Read box params from conf files
def read_conf(conf_path):
    box = {}
    for line in Path(conf_path).read_text().splitlines():
        line = line.strip()
        for k in ('center_x','center_y','center_z','size_x','size_y','size_z'):
            if line.startswith(k):
                box[k] = float(line.split('=')[1].split()[0])
    return box

df_vina = pd.read_csv(VINA_DIR / 'vina_tier2_all.csv')
df_vina['vina_score'] = pd.to_numeric(df_vina['vina_score'], errors='coerce')
df_vina = df_vina.dropna(subset=['vina_score'])

def _parse_gnina_cnn(pdbqt_path):
    if not Path(pdbqt_path).exists():
        return None
    for line in Path(pdbqt_path).read_text().splitlines():
        if 'CNNscore' in line:
            try: return float(line.split()[-1])
            except: pass
    return None

def rescore_one(args):
    tname, cid, receptor, box = args
    lig_vina = VINA_DIR / f'{tname}_{cid}_out.pdbqt'
    if not lig_vina.exists():
        lig_vina = PDBQT_DIR / f'{cid}.pdbqt'

    out_p = GNINA_OUT / f'{tname}_{cid}_gnina.pdbqt'
    if out_p.exists() and out_p.stat().st_size > 50:
        cnn = _parse_gnina_cnn(out_p)
        return cid, tname, cnn, 'skip'

    # Strip MODEL/ENDMDL wrapper from Vina output (gnina needs clean PDBQT)
    raw = lig_vina.read_text()
    if 'MODEL' in raw:
        lines = []
        in_model = False
        for line in raw.splitlines():
            if line.startswith('MODEL'): in_model = True; continue
            if line.startswith('ENDMDL'): break
            if in_model: lines.append(line)
        clean_p = GNINA_OUT / f'{tname}_{cid}_clean.pdbqt'
        clean_p.write_text('\n'.join(lines))
        lig_src = str(clean_p)
    else:
        lig_src = str(lig_vina)

    cmd = [
        GNINA_BIN,
        '--receptor', receptor,
        '--ligand',   lig_src,
        '--out',      str(out_p),
        '--center_x', str(box['center_x']),
        '--center_y', str(box['center_y']),
        '--center_z', str(box['center_z']),
        '--size_x',   str(box['size_x']),
        '--size_y',   str(box['size_y']),
        '--size_z',   str(box['size_z']),
        '--cnn', 'crossdock_default2018',
        '--no_gpu',
        '--num_modes', '1',
        '--cpu', '1',
    ]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return cid, tname, None, 'timeout'

    # Clean up temp file
    clean_tmp = GNINA_OUT / f'{tname}_{cid}_clean.pdbqt'
    if clean_tmp.exists():
        clean_tmp.unlink()

    cnn = _parse_gnina_cnn(out_p)
    return cid, tname, cnn, 'ok'

tasks = []
for tname, t in TARGETS.items():
    top50 = (df_vina[df_vina['target'] == tname]
             .nsmallest(50, 'vina_score')['compound'].tolist())
    box = read_conf(t['config'])
    for cid in top50:
        tasks.append((tname, cid, t['receptor'], box))

print(f"Gnina tasks: {len(tasks)}")

rows = []
with ProcessPoolExecutor(max_workers=NCPUS) as ex:
    futs = {ex.submit(rescore_one, task): task for task in tasks}
    for fut in as_completed(futs):
        cid, tname, cnn, status = fut.result()
        rows.append({'compound': cid, 'target': tname, 'cnn_score': cnn, 'status': status})

df_gnina = pd.DataFrame(rows)
df_gnina.to_csv(GNINA_OUT / 'gnina_tier2_top50.csv', index=False)
print(f"Written: gnina_tier2_top50.csv")

for tname in TARGETS:
    sub = df_gnina[(df_gnina['target'] == tname)].dropna(subset=['cnn_score'])
    print(f"\n{tname} top-10 by CNN:")
    print(sub.nlargest(10, 'cnn_score')[['compound','cnn_score']].to_string(index=False))
PYEOF

# ── Step 3: Consensus filter + shortlist ─────────────────────────────────────
echo ""
echo "=== Step 3: Consensus filter ==="

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python - << 'PYEOF'
import os
from pathlib import Path
import pandas as pd

PROJDIR  = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
VINA_CSV = PROJDIR / 'results/phase4/vina_tier2/vina_tier2_all.csv'
GNINA_CSV= PROJDIR / 'results/phase4/gnina_tier2/gnina_tier2_top50.csv'
OUT_CSV  = PROJDIR / 'results/phase4/tier2_consensus_hits.csv'
SHORT_CSV= PROJDIR / 'results/phase4/tier2_shortlist.csv'
TIER2_META = PROJDIR / 'data/compounds/tier2_approved.csv'

THRESHOLDS = {
    'IRF8_AF2_DBD': {'vina': -6.5, 'cnn': 0.70},
    'IKZF1_8RQC':   {'vina': -6.0, 'cnn': 0.70},
}

df_vina  = pd.read_csv(VINA_CSV).dropna(subset=['vina_score'])
df_gnina = pd.read_csv(GNINA_CSV).dropna(subset=['cnn_score'])
df_meta  = pd.read_csv(TIER2_META)[['molecule_chembl_id','pref_name','mw','alogp','cns_mpo']]
df_meta  = df_meta.rename(columns={'molecule_chembl_id':'compound'})

# Merge
df = df_vina.merge(df_gnina[['compound','target','cnn_score']], on=['compound','target'], how='left')
df = df.merge(df_meta, on='compound', how='left')

hits = []
for tname, thr in THRESHOLDS.items():
    sub = df[(df['target'] == tname) &
             (df['vina_score'] <= thr['vina']) &
             (df['cnn_score']  >= thr['cnn'])]
    hits.append(sub)
    print(f"{tname}: {len(sub)} consensus hits  "
          f"(Vina <= {thr['vina']}, CNN >= {thr['cnn']})")

df_hits = pd.concat(hits).sort_values(['target','vina_score'])
df_hits.to_csv(OUT_CSV, index=False)
print(f"\nAll consensus hits: {len(df_hits)}  -> {OUT_CSV}")

# Shortlist: top 10 per target
short = df_hits.groupby('target', group_keys=False).apply(
    lambda x: x.nsmallest(10, 'vina_score'))
short.to_csv(SHORT_CSV, index=False)

print(f"\n=== TIER-2 SHORTLIST (top 10 per target) ===")
print(short[['target','compound','pref_name','vina_score','cnn_score','mw','alogp','cns_mpo']].to_string(index=False))
PYEOF

echo ""
echo "Finished: $(date)"
echo ""
echo "Results:"
echo "  Vina scores:       $OUT_DIR/vina_tier2_all.csv"
echo "  Gnina scores:      $GNINA_OUT/gnina_tier2_top50.csv"
echo "  Consensus hits:    $PROJDIR/results/phase4/tier2_consensus_hits.csv"
echo "  Shortlist (top10): $PROJDIR/results/phase4/tier2_shortlist.csv"
