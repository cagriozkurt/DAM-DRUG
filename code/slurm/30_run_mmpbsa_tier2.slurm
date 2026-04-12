#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=mmpbsa_tier2
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --array=0-9   # 26 complexes in manifest

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
MANIFEST=$PROJDIR/results/phase4/mmpbsa_tier2/manifest.txt

mkdir -p "$PROJDIR/logs"

# Read work directory for this array task
WORKDIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
if [ -z "$WORKDIR" ]; then
    echo "ERROR: No entry at index $SLURM_ARRAY_TASK_ID in manifest"
    exit 1
fi

echo "Job:     $SLURM_JOB_ID  Array: $SLURM_ARRAY_JOB_ID[$SLURM_ARRAY_TASK_ID]"
echo "Node:    $SLURMD_NODENAME"
echo "Workdir: $WORKDIR"
echo "Start:   $(date)"

cd "$WORKDIR" || { echo "Cannot cd to $WORKDIR"; exit 1; }

# ── Initialize module environment ──────────────────────────────────────────
# TRUBA barbun batch jobs don't source /etc/profile.d automatically; we must
# initialize the module system explicitly before calling 'module load'.
source /etc/profile.d/modules.sh

module load apps/gromacs/2024.1-oneapi2024 || { echo "ERROR: Cannot load GROMACS module"; exit 1; }
module load miniconda3                      || { echo "ERROR: Cannot load miniconda3 module"; exit 1; }

# ── Activate conda environment (provides gmx_MMPBSA, acpype, AmberTools) ──
CONDA_BASE=$(conda info --base 2>/dev/null || true)
[ -z "$CONDA_BASE" ] && CONDA_BASE=$(dirname "$(dirname "$(which conda 2>/dev/null)")")
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate mmpbsa || { echo "ERROR: Cannot activate mmpbsa conda env"; exit 1; }

GMX=$(command -v gmx_mpi || command -v gmx || echo "")
[ -z "$GMX" ] && { echo "ERROR: GROMACS not found after module load"; exit 1; }
echo "GROMACS: $GMX"

command -v gmx_MMPBSA &>/dev/null || { echo "ERROR: gmx_MMPBSA not found — install via: pip install gmx_MMPBSA"; exit 1; }

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

# ── Step 0.5: Trim disordered termini from receptor ───────────────────────
# AF2 models (e.g. BHLHE41_AF2_bHLH) can have N/C-terminal tails that are
# 60-77 Å from the functional core.  These inflate the solvation box to
# >200k SOL molecules (solvation guard above would abort), and pdb2gmx/
# grompp then fail with "negative index" from integer overflow.
# Fix: remove any residue whose Cα is >50 Å from the Cα centroid.
# For compact crystal structures, no residues are removed (no-op).
echo "=== Step 0.5: receptor coordinate check/trim ==="
python3 -c "
import math
from pathlib import Path

pdb = Path('receptor.pdb').read_text().splitlines()

# Collect Cα positions per residue
ca_pos = {}
for line in pdb:
    if line.startswith('ATOM') and line[12:16].strip() == 'CA':
        try:
            resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ca_pos[resnum] = (x, y, z)
        except ValueError:
            pass

if len(ca_pos) < 10:
    print(f'Too few CA atoms ({len(ca_pos)}) — skipping trim')
    exit(0)

cx = sum(v[0] for v in ca_pos.values()) / len(ca_pos)
cy = sum(v[1] for v in ca_pos.values()) / len(ca_pos)
cz = sum(v[2] for v in ca_pos.values()) / len(ca_pos)
print(f'CA centroid: ({cx:.1f}, {cy:.1f}, {cz:.1f}) A')

CUTOFF = 50.0
keep = {r for r, (x,y,z) in ca_pos.items()
        if math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2) <= CUTOFF}

n_orig = len(ca_pos)
n_keep = len(keep)
removed = sorted(set(ca_pos.keys()) - keep)
print(f'Residues within {CUTOFF:.0f}A: {n_keep}/{n_orig}')

if not removed:
    print('Receptor is compact — no trimming needed')
    exit(0)

suffix = '...' if len(removed) > 10 else ''
print(f'Removing {len(removed)} distant residues: {removed[:10]}{suffix}')

out = []
for line in pdb:
    if line.startswith(('ATOM','HETATM','ANISOU')):
        try:
            resnum = int(line[22:26].strip())
            if resnum not in keep:
                continue
        except ValueError:
            pass
    elif line.startswith('TER'):
        # Only keep TER lines that follow a kept residue
        out.append(line)
        continue
    out.append(line)

Path('receptor.pdb').write_text('\n'.join(out) + '\n')
print(f'Trimmed receptor.pdb: {n_keep} residues retained')
"

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

# GRO coordinate fields: fixed-width columns 20-28, 28-36, 36-44 (nm)
def parse_xyz(line):
    if len(line) < 44: return 0.0, 0.0, 0.0
    return float(line[20:28]), float(line[28:36]), float(line[36:44])

def translate_line(line, dx, dy, dz):
    if len(line) < 44: return line
    x, y, z = parse_xyz(line)
    return line[:20] + f'{x-dx:8.3f}{y-dy:8.3f}{z-dz:8.3f}' + line[44:]

# Translate both protein and ligand so protein centroid is at origin.
# AF2 models can have large absolute-coordinate offsets that inflate the
# solvation box by orders of magnitude; centering before editconf fixes this.
prot_atoms = prot[2:-1]
lig_atoms  = lig[2:-1]
n = len(prot_atoms)
cx = sum(parse_xyz(l)[0] for l in prot_atoms) / n
cy = sum(parse_xyz(l)[1] for l in prot_atoms) / n
cz = sum(parse_xyz(l)[2] for l in prot_atoms) / n
prot_t = [translate_line(l, cx, cy, cz) for l in prot_atoms]
lig_t  = [translate_line(l, cx, cy, cz) for l in lig_atoms]
print(f'Centroid shift: ({cx:.2f}, {cy:.2f}, {cz:.2f}) nm')

lines = [f'{prot[0].strip()}+LIG'] + [str(n_total)] + prot_t + lig_t + [box]
Path('complex.gro').write_text('\n'.join(lines) + '\n')
print(f'Complex: {n_prot} protein + {n_lig} ligand = {n_total} atoms')
"

# ── Copy local forcefield & inject ligand atomtypes ───────────────────────
# GROMACS 2024 only allows [ atomtypes ] inside the forcefield ITP include
# chain.  Any standalone [ atomtypes ] block (even before the forcefield) is
# rejected.  Fix: copy amber99sb-ildn.ff locally, inject GAFF2 types into the
# local ffnonbonded.itp (where they're inside the forcefield chain), then
# include lig_mol.itp (without [ atomtypes ]) after the local forcefield.
FF_SRC=$(find "${GMXDATA:-}" -name "amber99sb-ildn.ff" -type d 2>/dev/null | head -1)
[ -z "$FF_SRC" ] && FF_SRC="${GMXLIB}/amber99sb-ildn.ff"
# Always re-copy a clean forcefield — never reuse a directory from a prior run,
# which may have atomtypes injected in the wrong position (after [ pairtypes ]).
rm -rf amber99sb-ildn.ff
cp -r "$FF_SRC" .
echo "Forcefield: $FF_SRC"

python3 -c "
import re
from pathlib import Path

itp = Path('lig_GMX.itp').read_text()

# ── Collect existing atomtype names from local ffnonbonded.itp ────────────
ffnb_path = Path('amber99sb-ildn.ff/ffnonbonded.itp')
ffnb_text = ffnb_path.read_text()
existing = set()
in_at = False
for line in ffnb_text.splitlines():
    s = line.strip()
    if re.match(r'\[\s*atomtypes\s*\]', s):
        in_at = True; continue
    elif s.startswith('['):
        in_at = False
    if in_at and s and not s.startswith(';'):
        existing.add(s.split()[0])

# ── Extract new GAFF2 types from lig_GMX.itp (skip duplicates) ───────────
new_at_lines = []
in_at = False
for line in itp.splitlines(keepends=True):
    if re.match(r'\[\s*atomtypes\s*\]', line):
        in_at = True; continue
    elif re.match(r'\[', line):
        in_at = False
    if in_at:
        s = line.strip()
        if s and not s.startswith(';'):
            name = s.split()[0]
            if name not in existing:
                new_at_lines.append(line)
                existing.add(name)

if new_at_lines:
    # Insert INSIDE the existing [ atomtypes ] section of ffnonbonded.itp.
    # GROMACS 2024 requires atomtypes to precede pairtypes/bondtypes/etc.
    # amber99sb-ildn.ff/ffnonbonded.itp has exactly one section ([ atomtypes ])
    # with no following [ ... ] header, so insert_idx will never be set by the
    # loop.  We handle three cases:
    #   1. [ atomtypes ] followed by another section  → insert before that section
    #   2. [ atomtypes ] is the last/only section     → append at end of file
    #   3. No [ atomtypes ] at all                    → append a new section (rare)
    ffnb_lines = ffnb_text.splitlines(keepends=True)
    insert_idx = None
    found_at = False
    in_at2 = False
    for idx, ln in enumerate(ffnb_lines):
        s2 = ln.strip()
        if re.match(r'\[\s*atomtypes\s*\]', s2):
            in_at2 = True
            found_at = True
        elif in_at2 and s2.startswith('['):
            insert_idx = idx   # insert before next section
            break
    if insert_idx is not None:
        # Case 1: insert before the next section header
        ffnb_lines[insert_idx:insert_idx] = new_at_lines
        ffnb_path.write_text(''.join(ffnb_lines))
        print(f'Injected {len(new_at_lines)} atomtypes inside [ atomtypes ] (before next section)')
    elif found_at:
        # Case 2: [ atomtypes ] is the last section — append directly (no new header)
        ffnb_path.write_text(ffnb_text.rstrip('\n') + '\n' + ''.join(new_at_lines))
        print(f'Injected {len(new_at_lines)} atomtypes at end of [ atomtypes ] section')
    else:
        # Case 3: no [ atomtypes ] at all — create one
        with open(ffnb_path, 'a') as f:
            f.write('\n[ atomtypes ]\n')
            f.writelines(new_at_lines)
        print(f'Injected {len(new_at_lines)} atomtypes (new [ atomtypes ] section created)')
else:
    print('No new atomtypes (all already in forcefield)')

# ── Write lig_mol.itp — lig_GMX.itp without [ atomtypes ] ────────────────
mol_lines = []
in_at = False
for line in itp.splitlines(keepends=True):
    if re.match(r'\[\s*atomtypes\s*\]', line):
        in_at = True
    elif re.match(r'\[', line):
        in_at = False
    if not in_at:
        mol_lines.append(line)
Path('lig_mol.itp').write_text(''.join(mol_lines))
print(f'lig_mol.itp: {len(mol_lines)} lines (no [ atomtypes ])')
"

# ── Build topology ─────────────────────────────────────────────────────────
python3 -c "
import re
from pathlib import Path

top = Path('protein.top').read_text()

# Insert lig_mol.itp immediately AFTER the forcefield include.
# pdb2gmx writes './amber99sb-ildn.ff/forcefield.itp' (with ./) when the ff
# is found in the current directory, so match both forms with a regex.
top, n = re.subn(
    r'(#include\s+\"\.?/?amber99sb-ildn\.ff/forcefield\.itp\")',
    r'\1\n#include \"lig_mol.itp\"',
    top, count=1)
if n == 0:
    print('WARNING: forcefield include not found — inserting lig_mol.itp at top (may fail)')
    top = '#include \"lig_mol.itp\"\n' + top
else:
    print('Inserted #include \"lig_mol.itp\" after forcefield include')

# Read the EXACT moleculetype name from lig_mol.itp (case-sensitive).
# parmed is case-sensitive: if [ molecules ] has 'LIG' but [ moleculetype ]
# defines 'lig', parmed raises 'no template defined'.  Use the actual name.
lig_itp = Path('lig_mol.itp').read_text()
m = re.search(r'\[\s*moleculetype\s*\][^\[]*?\n\s*([A-Za-z0-9_]+)', lig_itp)
mol_name = m.group(1) if m else 'LIG'
print(f'Ligand moleculetype name: {mol_name}')

top = top.rstrip() + f'\n{mol_name}              1\n'
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

# Sanity-check solvation: abort early if system is abnormally large
# (>200k SOL molecules almost always means the receptor had extreme coordinates)
SOL_COUNT=$(python3 -c "
import re
txt = open('topol.top').read()
m = re.search(r'^SOL\s+(\d+)', txt, re.MULTILINE)
print(m.group(1) if m else '0')
" 2>/dev/null || echo "0")
echo "Solvated system: $SOL_COUNT SOL molecules"
if [ "$SOL_COUNT" -gt 200000 ]; then
    echo "ERROR: Abnormally large system ($SOL_COUNT SOL). Check receptor.pdb for extreme coordinates."
    exit 1
fi

# ── Step 5: Add ions (0.15 M NaCl) ────────────────────────────────────────
echo "=== Step 5: genion ==="
$GMX grompp \
    -f em.mdp \
    -c complex_solv.gro \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 5 \
    -nobackup > grompp_ions.log 2>&1 || { echo "ERROR: grompp for ions failed:"; tail -40 grompp_ions.log; exit 1; }

echo "SOL" | $GMX genion \
    -s ions.tpr \
    -o complex_ions.gro \
    -p topol.top \
    -pname NA -nname CL \
    -neutral -conc 0.15 \
    -nobackup > genion.log 2>&1 || { echo "ERROR: genion failed:"; tail -20 genion.log; exit 1; }

# ── Invalidate stale MD checkpoints (system rebuilt from scratch above) ────
# Steps 1-5 always rebuild pdb2gmx→genion, so old nvt/npt/md files from a
# prior run may have a different atom count or ion configuration.
echo "Clearing stale MD/MMPBSA files from prior runs"
rm -f nvt.cpt nvt.gro nvt.tpr nvt.edr nvt.log nvt.xtc
rm -f npt.cpt npt.gro npt.tpr npt.edr npt.log npt.xtc
rm -f md.xtc md.cpt md.tpr md.edr md.log md_whole.xtc md_run.mdp
rm -f FINAL_RESULTS_MMPBSA.dat mmpbsa_result.txt MMPBSA_DONE
rm -f COM.prmtop REC.prmtop LIG.prmtop

# ── Step 6: Energy minimisation ───────────────────────────────────────────
echo "=== Step 6: energy minimisation ==="
$GMX grompp \
    -f em.mdp \
    -c complex_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 5 \
    -nobackup > grompp_em.log 2>&1 || { echo "ERROR: grompp EM failed:"; tail -40 grompp_em.log; exit 1; }

$GMX mdrun \
    -v -deffnm em \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5
[ -f "em.gro" ] || { echo "ERROR: EM mdrun failed — em.gro missing"; exit 1; }

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
        # Use 7-char field so 6-digit atom numbers (>99999) still get a
        # leading space and are whitespace-separated in the index file.
        # With 6d, atoms >=100000 run together (e.g. '100001100002') and
        # GROMACS's free-format parser reads them as one token that
        # overflows int32 -> negative index -2 / Invalid atom number -1.
        f.write(f'{a:7d}')
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
    -maxwarn 5 \
    -nobackup > grompp_nvt.log 2>&1 || { echo "ERROR: grompp NVT failed:"; tail -40 grompp_nvt.log; exit 1; }

$GMX mdrun \
    -deffnm nvt \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5
[ -f "nvt.cpt" ] || { echo "ERROR: NVT mdrun failed — nvt.cpt missing"; exit 1; }

# ── Step 9: NPT equilibration (100 ps) ────────────────────────────────────
echo "=== Step 9: NPT ==="
# Always rewrite npt.mdp to guarantee Berendsen barostat (stable during equil)
cat > npt.mdp << 'NPTMDP'
; NPT equilibration (100 ps, Berendsen barostat, NO position restraints)
; Position restraints are deliberately OFF here: -DPOSRES would cause the
; restraint forces to dominate the pressure tensor (producing millions of bar)
; which leads to catastrophic Berendsen box-rescaling and SETTLE failures.
; Without restraints, the initial box pressure is a few hundred bar at most,
; well within Berendsen's stable range.
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 10
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
tcoupl      = V-rescale
tc-grps     = Protein_LIG  SOL_Ion
tau_t       = 0.1          0.1
ref_t       = 300          300
pcoupl      = Berendsen
pcoupltype  = isotropic
tau_p       = 1.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
gen_vel     = no
NPTMDP

$GMX grompp \
    -f npt.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -n index.ndx \
    -o npt.tpr \
    -maxwarn 5 \
    -nobackup > grompp_npt.log 2>&1 || { echo "ERROR: grompp NPT failed:"; tail -40 grompp_npt.log; exit 1; }

$GMX mdrun \
    -deffnm npt \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5
[ -f "npt.gro" ] || { echo "ERROR: NPT mdrun failed — npt.gro missing"; exit 1; }

# ── Step 10: Production MD (1 ns) ─────────────────────────────────────────
echo "=== Step 10: production MD ==="
# mmpbsa.in requires startframe=251 through endframe=500 (500 frames total).
# md.mdp ships with nstxout-compressed=5000 (every 10 ps → only 100 frames).
# Patch to nstxout-compressed=1000 (every 2 ps → 500 frames) without touching
# the user's original md.mdp.
sed 's/nstxout-compressed\s*=\s*[0-9]*/nstxout-compressed = 1000/' md.mdp > md_run.mdp
echo "md_run.mdp: nstxout-compressed patched to 1000 (500 frames over 1 ns)"

# Always run grompp — md.tpr is required by gmx_MMPBSA even when xtc exists
$GMX grompp \
    -f md_run.mdp \
    -c npt.gro \
    -t npt.cpt \
    -p topol.top \
    -n index.ndx \
    -o md.tpr \
    -maxwarn 5 \
    -nobackup > grompp_md.log 2>&1 || { echo "ERROR: grompp production MD failed:"; tail -40 grompp_md.log; exit 1; }

if [ -f "md.xtc" ]; then echo "md.xtc exists — skipping production mdrun"; else
$GMX mdrun \
    -deffnm md \
    -ntomp $SLURM_CPUS_PER_TASK \
    -nobackup 2>&1 | tail -5
[ -f "md.xtc" ] || { echo "ERROR: production mdrun failed — md.xtc missing"; exit 1; }
fi  # end production mdrun skip

# ── Step 10.5: Make trajectory whole (remove PBC split-molecule artefacts) ─
# sander GB calculations are non-periodic.  If the protein crosses a box
# boundary during MD, sander sees a ~90 Å "bond" instead of 1.5 Å, producing
# BOND energies of ~86 million kcal/mol (even though force constants and
# equilibrium distances in the prmtop are correct).  Fix: image the complex
# as a single whole molecule before passing to gmx_MMPBSA.
echo "=== Step 10.5: image trajectory (pbc mol) ==="
# Select Protein_LIG as output (strips solvent/ions — gmx_MMPBSA only needs complex)
# -pbc mol makes each molecule whole across PBC without needing a centering group
echo "Protein_LIG" | $GMX trjconv \
    -f md.xtc \
    -s md.tpr \
    -n index.ndx \
    -o md_whole.xtc \
    -pbc mol \
    -nobackup 2>&1 | tail -3
[ -f "md_whole.xtc" ] || { echo "ERROR: trjconv imaging failed"; exit 1; }

# ── Step 11: gmx_MMPBSA ───────────────────────────────────────────────────
echo "=== Step 11: MM-GBSA ==="
echo "Force field: amber99sb-ildn.ff (local, with injected GAFF2 types)"

# Receptor = group 1 (Protein), Ligand = group 2 (LIG) — matches Python index above
gmx_MMPBSA \
    -O \
    -i mmpbsa.in \
    -cs md.tpr \
    -ci index.ndx \
    -cg 1 2 \
    -ct md_whole.xtc \
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
# gmx_MMPBSA writes 'ΔTOTAL  avg  SD(Prop)  SD  SEM(Prop)  SEM'
# Columns: Average, SD(Prop.), SD, SEM(Prop.), SEM
m = re.search(r'ΔTOTAL\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', txt)
if m:
    dg  = float(m.group(1))   # Average ΔG
    sem = float(m.group(4))   # SEM(Prop.) — propagated SEM; more informative than
                               # column-5 SEM which rounds to 0.00 for near-zero ΔG
    meta   = Path('meta.txt').read_text()
    target = [l.split('=')[1] for l in meta.splitlines() if l.startswith('target')][0]
    cid    = [l.split('=')[1] for l in meta.splitlines() if l.startswith('chembl')][0]
    vina   = [l.split('=')[1] for l in meta.splitlines() if l.startswith('vina')][0]
    print(f'RESULT: {target},{cid},{vina},{dg},{sem}')
    Path('mmpbsa_result.txt').write_text(f'{target},{cid},{vina},{dg},{sem}\n')
else:
    print('ERROR: Could not parse ΔTOTAL from results')
    print(txt[-500:])  # print end of file for debugging
"

touch MMPBSA_DONE
echo "Done: $(date)"
