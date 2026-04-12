#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=core_rmsd
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --array=0-3
# ─────────────────────────────────────────────────────────────────────────────
# Core-RMSD reanalysis for 100 ns MD results
#
# Rationale: IRF8 is an AlphaFold2 prediction; residues 136-164 are a
# flexible C-terminus (per-residue RMSF up to 14.7 Å) that dominates the
# backbone RMSD calculation and inflates the ligand RMSD (via bad fitting).
# This script:
#   1. Reads the per-residue RMSF (already computed in 24_run_md100ns.slurm)
#   2. Identifies stable-core residues: RMSF < RMSF_CUTOFF_NM (default 0.4 nm)
#   3. Builds a Core_Backbone index group from backbone atoms of those residues
#   4. Recomputes backbone RMSD (fit=Core_Backbone, group=Core_Backbone)
#   5. Recomputes ligand RMSD (fit=Core_Backbone, group=LIG)
#   6. Re-reads existing hbond_num.xvg (mindist proxy, already computed)
#   7. Writes md100ns_core_result.txt with the corrected metrics
#
# Uses the PBC-corrected trajectory (md100ns_whole.xtc) already produced
# by 24_run_md100ns.slurm Step 11.
#
# Acceptance criteria (same as primary analysis):
#   Core backbone RMSD < 3.0 Å  (last-50 ns average)
#   Ligand RMSD       < 2.0 Å  (last-50 ns average, fit on Core_Backbone)
#   Contact occupancy > 30%    (mindist ≤ 3.5 Å, from existing hbond_num.xvg)
#
# Same manifest as 24_run_md100ns.slurm:
#   $PROJDIR/results/phase4/md/md100ns_manifest.txt
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
MANIFEST=$PROJDIR/results/phase4/md/md100ns_manifest.txt

mkdir -p "$PROJDIR/logs"

SRC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
if [ -z "$SRC" ]; then
    echo "ERROR: No entry at index $SLURM_ARRAY_TASK_ID in manifest"
    exit 1
fi

TARGET=$(basename "$(dirname "$SRC")")
CHEMBL=$(basename "$SRC")
WORKDIR=$PROJDIR/results/phase4/md/$TARGET/$CHEMBL

echo "Job:     $SLURM_JOB_ID  Array: $SLURM_ARRAY_JOB_ID[$SLURM_ARRAY_TASK_ID]"
echo "Node:    $SLURMD_NODENAME"
echo "Workdir: $WORKDIR"
echo "Start:   $(date)"

cd "$WORKDIR" || { echo "Cannot cd to $WORKDIR"; exit 1; }

# Skip if already done
if [ -f "MD100NS_CORE_DONE" ]; then
    echo "Core RMSD analysis already complete — skipping"
    exit 0
fi

# Verify prerequisites
for f in rmsf.xvg md100ns_whole.xtc md100ns.tpr em.gro index.ndx hbond_num.xvg; do
    [ -f "$f" ] || { echo "ERROR: missing $f — run 24_run_md100ns.slurm first"; exit 1; }
done

# ── Module / environment ────────────────────────────────────────────────────
source /etc/profile.d/modules.sh
module load apps/gromacs/2024.1-oneapi2024 || { echo "ERROR: Cannot load GROMACS"; exit 1; }
GMX=$(command -v gmx_mpi || command -v gmx || echo "")
[ -z "$GMX" ] && { echo "ERROR: GROMACS not found"; exit 1; }
echo "GROMACS: $GMX"

# ── Step 1: identify stable core + write Core_Backbone index group ──────────
echo "=== Step 1: stable core from RMSF ==="
python3 -c "
import re
from pathlib import Path

RMSF_CUTOFF_NM = 0.40   # 4.0 Å — residues above this are excluded from core
MIN_CORE_RESIDUES = 20  # fall back to full backbone if fewer than this

PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
               'HIS','HID','HIE','HIP','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER',
               'THR','TRP','TYR','VAL'}
BACKBONE_AT = {'N','CA','C','O'}

# ── read rmsf.xvg (residue index → RMSF in nm) ──────────────────────────
stable_res = set()
for line in Path('rmsf.xvg').read_text().splitlines():
    if line.startswith(('@','#')): continue
    parts = line.split()
    if len(parts) < 2: continue
    try:
        res_idx = int(float(parts[0]))
        rmsf    = float(parts[1])
        if rmsf < RMSF_CUTOFF_NM:
            stable_res.add(res_idx)
    except ValueError:
        pass

print(f'RMSF < {RMSF_CUTOFF_NM*10:.1f} A: {len(stable_res)} stable residues')

# ── map sequential GROMACS residue index → backbone atom numbers ─────────
# em.gro residues are in the same order as GROMACS internal numbering.
lines = Path('em.gro').read_text().splitlines()
n_atoms = int(lines[1].strip())

res_idx   = 0
prev_key  = None    # (resnum_str, resname) to detect residue boundaries
backbone  = []

for i, line in enumerate(lines[2:2+n_atoms], 1):
    if len(line) < 15: continue
    resname  = line[5:10].strip().upper()
    resnum   = line[0:5].strip()
    atname   = line[10:15].strip()
    if resname not in PROTEIN_RES: continue
    key = (resnum, resname)
    if key != prev_key:
        res_idx += 1
        prev_key = key
    if res_idx in stable_res and atname in BACKBONE_AT:
        backbone.append(i)

print(f'Core_Backbone: {len(backbone)} atoms from {len(stable_res)} residues')

if len(stable_res) < MIN_CORE_RESIDUES:
    print(f'WARNING: fewer than {MIN_CORE_RESIDUES} stable residues — falling back to full Backbone')
    # read existing Backbone group from index.ndx and rewrite Core_Backbone = Backbone
    idx_text = Path('index.ndx').read_text()
    in_bb = False
    backbone = []
    for line in idx_text.splitlines():
        if re.match(r'\[\s*Backbone\s*\]', line): in_bb = True; continue
        elif line.strip().startswith('['):         in_bb = False
        if in_bb:
            backbone.extend(int(x) for x in line.split() if x.isdigit())

# ── check Core_Backbone not already in index.ndx ─────────────────────────
idx_text = Path('index.ndx').read_text()
if '[ Core_Backbone ]' in idx_text:
    print('Core_Backbone already in index.ndx — skipping append')
else:
    with open('index.ndx', 'a') as f:
        f.write('[ Core_Backbone ]\n')
        for j, a in enumerate(backbone):
            f.write(f'{a:7d}')
            if (j+1) % 15 == 0: f.write('\n')
        f.write('\n\n')
    print('Core_Backbone appended to index.ndx')
"

# ── Step 2: backbone RMSD on core (fit=Core_Backbone, group=Core_Backbone) ─
echo "=== Step 2: core backbone RMSD ==="
printf "Core_Backbone\nCore_Backbone\n" | $GMX rms \
    -s md100ns.tpr -f md100ns_whole.xtc -n index.ndx \
    -o core_backbone_rmsd.xvg -tu ns \
    -nobackup 2>&1 | tail -3

# ── Step 3: ligand RMSD fitting on core backbone ────────────────────────────
echo "=== Step 3: ligand RMSD (fit on Core_Backbone) ==="
printf "Core_Backbone\nLIG\n" | $GMX rms \
    -s md100ns.tpr -f md100ns_whole.xtc -n index.ndx \
    -o core_ligand_rmsd.xvg -tu ns \
    -nobackup 2>&1 | tail -3

# ── Step 4: summary ──────────────────────────────────────────────────────────
echo "=== Step 4: summary ==="
python3 -c "
from pathlib import Path

def read_xvg(path):
    data = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(('@','#')): continue
        parts = line.split()
        if len(parts) >= 2:
            try: data.append((float(parts[0]), float(parts[1])))
            except ValueError: pass
    return data

def last_half_mean(data):
    half = data[len(data)//2:]
    return sum(v for _,v in half) / max(1, len(half))

def last_half_max(data):
    half = data[len(data)//2:]
    return max((v for _,v in half), default=0.0)

bb  = read_xvg('core_backbone_rmsd.xvg')
lig = read_xvg('core_ligand_rmsd.xvg')
hb  = read_xvg('hbond_num.xvg')   # reuse mindist proxy from primary analysis

# GROMACS rms outputs nm; × 10 → Å
bb_avg_A   = last_half_mean(bb)  * 10
lig_avg_A  = last_half_mean(lig) * 10
lig_max_A  = last_half_max(lig)  * 10

# Fraction of frames with ≥ 1 protein–ligand close contact (≤ 3.5 Å)
hb_occ = sum(1 for _,v in hb if v >= 1) / max(1, len(hb))

pass_bb  = bb_avg_A  < 3.0
pass_lig = lig_avg_A < 2.0
pass_hb  = hb_occ   > 0.30
verdict  = 'PASS' if (pass_bb and pass_lig and pass_hb) else 'FAIL'

meta   = Path('meta.txt').read_text()
target = [l.split('=')[1] for l in meta.splitlines() if l.startswith('target')][0]
cid    = [l.split('=')[1] for l in meta.splitlines() if l.startswith('chembl')][0]
vina   = [l.split('=')[1] for l in meta.splitlines() if l.startswith('vina')][0]
dg     = Path('$SRC/mmpbsa_result.txt').read_text().split(',')[3]

result = (f'{target},{cid},{vina},{dg},'
          f'{bb_avg_A:.2f},{lig_avg_A:.2f},{lig_max_A:.2f},{hb_occ:.2f},{verdict}')
print(f'RESULT: {result}')
Path('md100ns_core_result.txt').write_text(result + '\n')

print()
print(f'  Target / Compound   : {target} / {cid}')
print(f'  Vina score          : {vina} kcal/mol')
print(f'  MM-GBSA ΔG          : {dg} kcal/mol')
print(f'  Core BB RMSD (avg)  : {bb_avg_A:.2f} A  [<3.0: {pass_bb}]')
print(f'  Ligand RMSD (avg)   : {lig_avg_A:.2f} A  [<2.0: {pass_lig}]')
print(f'  Ligand RMSD (max)   : {lig_max_A:.2f} A')
print(f'  Contact occupancy   : {hb_occ:.0%}  [>30%: {pass_hb}]')
print(f'  Verdict             : {verdict}')
"

touch MD100NS_CORE_DONE
echo "Done: $(date)"
