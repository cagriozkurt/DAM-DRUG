#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=mmpbsa
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/mmpbsa-%A_%a.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/mmpbsa-%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --array=0-23   # 24 complexes prepared by 17_prep_mmpbsa.py

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
MANIFEST=$PROJDIR/results/phase4/mmpbsa/manifest.txt

# Read work directory for this array task
WORKDIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $MANIFEST)
if [ -z "$WORKDIR" ]; then
    echo "ERROR: No entry at index $SLURM_ARRAY_TASK_ID in manifest"
    exit 1
fi

echo "Job:     $SLURM_JOB_ID  Array: $SLURM_ARRAY_JOB_ID[$SLURM_ARRAY_TASK_ID]"
echo "Node:    $SLURMD_NODENAME"
echo "Workdir: $WORKDIR"
echo "Start:   $(date)"

cd $WORKDIR || { echo "Cannot cd to $WORKDIR"; exit 1; }

# ── Load modules ───────────────────────────────────────────────────────────
module load miniconda3
source activate mmpbsa   # provides tleap, acpype, gmx_MMPBSA

module load gromacs 2>/dev/null || \
    module load GROMACS 2>/dev/null || \
    echo "WARNING: gromacs module not found — trying PATH"

GMX=$(command -v gmx_mpi || command -v gmx || echo "")
if [ -z "$GMX" ]; then
    echo "ERROR: GROMACS not found"; exit 1
fi
echo "GROMACS: $GMX"

# gmx_MMPBSA and acpype should be pip-installed
MMPBSA=$(command -v gmx_MMPBSA || python3 -m gmx_MMPBSA 2>/dev/null && echo "python3 -m gmx_MMPBSA" || echo "")
if [ -z "$MMPBSA" ]; then
    echo "ERROR: gmx_MMPBSA not found. Install: pip3 install gmx_MMPBSA --user"
    exit 1
fi

# ── Check prerequisites ────────────────────────────────────────────────────
for f in receptor.pdb lig_GMX.itp lig_GMX.gro em.mdp nvt.mdp npt.mdp md.mdp mmpbsa.in; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing prerequisite: $f"
        exit 1
    fi
done

# Skip if already completed
if [ -f "MMPBSA_DONE" ]; then
    echo "Already complete — skipping"
    exit 0
fi

set -eo pipefail  # exit on any error, including inside pipes

# ── Step 1: Protein topology ───────────────────────────────────────────────
echo "=== Step 1: pdb2gmx ==="
$GMX pdb2gmx \
    -f receptor.pdb \
    -o protein.gro \
    -p protein.top \
    -ignh \
    -ff amber99sb-ildn \
    -water tip3p \
    -nobackup 2>&1 | tail -5

# ── Step 2: Combine protein + ligand ──────────────────────────────────────
echo "=== Step 2: combine protein + ligand ==="
# Merge coordinates: protein.gro + lig_GMX.gro → complex.gro
python3 -c "
from pathlib import Path

prot = Path('protein.gro').read_text().splitlines()
lig  = Path('lig_GMX.gro').read_text().splitlines()

# atom counts
n_prot = int(prot[1].strip())
n_lig  = int(lig[1].strip())
n_total = n_prot + n_lig

# box line = last line of protein.gro
box = prot[-1]

# atom lines: prot[2:-1]  lig[2:-1]
lines = [f'{prot[0].strip()}+LIG'] + [str(n_total)] + \
        prot[2:-1] + lig[2:-1] + [box]
Path('complex.gro').write_text('\n'.join(lines) + '\n')
print(f'Complex: {n_prot} protein + {n_lig} ligand = {n_total} atoms')
"

# Build topology: protein.top + ligand itp include
python3 -c "
from pathlib import Path

top = Path('protein.top').read_text()

# In GROMACS 2020+, [ atomtypes ] from ligand ITP must be defined BEFORE
# any [ moleculetype ] blocks.  Insert lig_GMX.itp immediately after the
# forcefield include (before the first [ moleculetype ]).
ff_line = '#include \"amber99sb-ildn.ff/forcefield.itp\"'
if ff_line in top:
    top = top.replace(ff_line,
                      ff_line + '\n#include \"lig_GMX.itp\"', 1)
else:
    # Fallback: prepend at very top
    top = '#include \"lig_GMX.itp\"\n' + top

# Add LIG molecule to [ molecules ] section
top = top.rstrip() + '\nLIG              1\n'

Path('topol.top').write_text(top)
print('topol.top written')
"

# ── Step 3: Define box ─────────────────────────────────────────────────────
echo "=== Step 3: editconf ==="
$GMX editconf \
    -f complex.gro \
    -o complex_box.gro \
    -c -d 1.0 -bt dodecahedron \
    -nobackup 2>&1 | tail -3

# ── Step 4: Solvate ────────────────────────────────────────────────────────
echo "=== Step 4: solvate ==="
$GMX solvate \
    -cp complex_box.gro \
    -cs spc216.gro \
    -o complex_solv.gro \
    -p topol.top \
    -nobackup 2>&1 | tail -3

# ── Step 5: Add ions (0.15 M NaCl) ────────────────────────────────────────
echo "=== Step 5: genion ==="
$GMX grompp \
    -f em.mdp \
    -c complex_solv.gro \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 2 \
    -nobackup 2>&1 | tail -5

echo "SOL" | $GMX genion \
    -s ions.tpr \
    -o complex_ions.gro \
    -p topol.top \
    -pname NA -nname CL \
    -neutral -conc 0.15 \
    -nobackup 2>&1 | tail -5

# ── Step 6: Energy minimisation ───────────────────────────────────────────
echo "=== Step 6: energy minimisation ==="
$GMX grompp \
    -f em.mdp \
    -c complex_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 2 \
    -nobackup 2>&1 | tail -3

$GMX mdrun \
    -v -deffnm em \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5

# ── Step 7: Make index & position restraints ──────────────────────────────
echo "=== Step 7: index ==="
# Write index.ndx directly from em.gro — avoids interactive make_ndx piping.
# Groups created (0-based):
#   0 System  1 Protein  2 LIG  3 Protein_LIG  4 SOL  5 ION  6 non-Protein  7 SOL_Ion
python3 -c "
from pathlib import Path

PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
               'HIS','HID','HIE','HIP','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER',
               'THR','TRP','TYR','VAL'}
ION_RES = {'NA','CL','MG','CA','ZN','K','NA+','CL-'}

lines = Path('em.gro').read_text().splitlines()
n = int(lines[1].strip())
prot, lig, sol, ion = [], [], [], []
lig_resnames = set()
for i, line in enumerate(lines[2:2+n], 1):
    if len(line) < 15: continue
    res = line[5:10].strip()
    res_up = res.upper()
    if res_up in PROTEIN_RES:        prot.append(i)
    elif res_up == 'SOL':            sol.append(i)
    elif res_up in ION_RES:          ion.append(i)
    else:
        lig.append(i)   # anything not protein/SOL/ion = ligand
        lig_resnames.add(res)
print(f'Ligand residue name(s) found: {lig_resnames}')

def grp(f, name, atoms):
    f.write(f'[ {name} ]\n')
    for j, a in enumerate(atoms):
        f.write(f'{a:6d}')
        if (j+1) % 15 == 0: f.write('\n')
    f.write('\n\n')

with open('index.ndx', 'w') as f:
    grp(f, 'System',      list(range(1, n+1)))
    grp(f, 'Protein',     prot)
    grp(f, 'LIG',         lig)
    grp(f, 'Protein_LIG', prot + lig)
    grp(f, 'SOL',         sol)
    grp(f, 'ION',         ion)
    grp(f, 'non-Protein', lig + sol + ion)
    grp(f, 'SOL_Ion',     sol + ion)

print(f'index.ndx: Protein={len(prot)}, LIG={len(lig)}, SOL={len(sol)}, ION={len(ion)}')
"

# Generate position restraints for protein+ligand heavy atoms
echo "Protein_LIG" | $GMX genrestr \
    -f em.gro \
    -n index.ndx \
    -o posres_complex.itp \
    -fc 1000 1000 1000 \
    -nobackup 2>&1 | tail -3

# ── Step 8: NVT equilibration (100 ps) ────────────────────────────────────
echo "=== Step 8: NVT ==="
$GMX grompp \
    -f nvt.mdp \
    -c em.gro \
    -r em.gro \
    -p topol.top \
    -n index.ndx \
    -o nvt.tpr \
    -maxwarn 2 \
    -nobackup 2>&1 | tail -3

$GMX mdrun \
    -deffnm nvt \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5

# ── Step 9: NPT equilibration (100 ps) ────────────────────────────────────
echo "=== Step 9: NPT ==="
$GMX grompp \
    -f npt.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -n index.ndx \
    -o npt.tpr \
    -maxwarn 2 \
    -nobackup 2>&1 | tail -3

$GMX mdrun \
    -deffnm npt \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5

# ── Step 10: Production MD (1 ns) ─────────────────────────────────────────
echo "=== Step 10: production MD ==="
$GMX grompp \
    -f md.mdp \
    -c npt.gro \
    -t npt.cpt \
    -p topol.top \
    -n index.ndx \
    -o md.tpr \
    -maxwarn 2 \
    -nobackup 2>&1 | tail -3

$GMX mdrun \
    -deffnm md \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5

# ── Step 11: gmx_MMPBSA ───────────────────────────────────────────────────
echo "=== Step 11: MM-GBSA ==="
# parmed (inside gmx_MMPBSA) resolves #include "amber99sb-ildn.ff/forcefield.itp"
# relative to the work directory — copy the ff from the GROMACS installation.
FF_SRC=$(find /arf/sw/apps/gromacs -name "amber99sb-ildn.ff" -type d 2>/dev/null | head -1)
[ -z "$FF_SRC" ] && FF_SRC="${GMXLIB}/amber99sb-ildn.ff"
cp -r "$FF_SRC" . 2>/dev/null || true
echo "Force field: $FF_SRC"

# Receptor = group 1 (Protein), Ligand = group 2 (LIG) — matches Python index above
gmx_MMPBSA \
    -O \
    -i mmpbsa.in \
    -cs md.tpr \
    -ci index.ndx \
    -cg 1 2 \
    -ct md.xtc \
    -cp topol.top \
    -nogui 2>&1

# ── Collect result ─────────────────────────────────────────────────────────
echo "=== Extracting ΔG ==="
python3 -c "
from pathlib import Path
import re

# gmx_MMPBSA writes FINAL_RESULTS_MMPBSA.dat
result_file = Path('FINAL_RESULTS_MMPBSA.dat')
if not result_file.exists():
    print('ERROR: FINAL_RESULTS_MMPBSA.dat not found')
    exit(1)

txt = result_file.read_text()
# Parse: 'DELTA TOTAL' line
m = re.search(r'DELTA TOTAL\s+([-\d.]+)\s+\+/-\s+([\d.]+)', txt)
if m:
    dg     = float(m.group(1))
    std    = float(m.group(2))
    meta   = Path('meta.txt').read_text()
    target = [l.split('=')[1] for l in meta.splitlines() if l.startswith('target')][0]
    cid    = [l.split('=')[1] for l in meta.splitlines() if l.startswith('chembl')][0]
    vina   = [l.split('=')[1] for l in meta.splitlines() if l.startswith('vina')][0]
    print(f'RESULT: {target},{cid},{vina},{dg},{std}')
    Path('mmpbsa_result.txt').write_text(f'{target},{cid},{vina},{dg},{std}\n')
else:
    print('ERROR: Could not parse DELTA TOTAL from results')
"

touch MMPBSA_DONE
echo "Done: $(date)"
