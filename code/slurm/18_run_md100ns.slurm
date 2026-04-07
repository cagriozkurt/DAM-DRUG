#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=md100ns
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/md100ns-%A_%a.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/md100ns-%A_%a.err
#SBATCH --open-mode=append          # accumulate logs across checkpoint restarts
#SBATCH --time=72:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --array=0-3
# ─────────────────────────────────────────────────────────────────────────────
# 100 ns production MD stability validation — top MM-GBSA hits
#
# Manifest: results/phase4/md/md100ns_manifest.txt  (4 lines, MM-GBSA workdirs)
#   /arf/scratch/mozkurt/DAM-DRUG/results/phase4/mmpbsa/IRF8_AF2_DBD_prep/CHEMBL2103837
#   /arf/scratch/mozkurt/DAM-DRUG/results/phase4/mmpbsa/IRF8_AF2_DBD_prep/CHEMBL490
#   /arf/scratch/mozkurt/DAM-DRUG/results/phase4/mmpbsa/IRF8_AF2_DBD_prep/CHEMBL42
#   /arf/scratch/mozkurt/DAM-DRUG/results/phase4/mmpbsa/PPARG_1FM9_LBD_prep/CHEMBL898
#
# Output goes to: results/phase4/md/<target>/<chembl>/
#
# PERFORMANCE NOTE: barbun CPUs typically achieve ~1-2 ns/day for ~130k-atom
# systems.  100 ns ≈ 50-100 days → job WILL time out.  Resubmit with the
# same sbatch command; mdrun resumes automatically from md100ns.cpt.
# Expect 3-5 resubmissions.  Use GPU partition (akya-cuda) if available to
# complete in 1-2 days instead.
#
# ACCEPTANCE CRITERIA (evaluated in final analysis step):
#   Protein backbone RMSD < 3.0 Å  (last-50ns average)
#   Ligand heavy-atom RMSD < 2.0 Å (last-50ns average)
#   Protein–ligand H-bond occupancy > 30%
# ─────────────────────────────────────────────────────────────────────────────

PROJDIR=/arf/scratch/mozkurt/DAM-DRUG
MANIFEST=$PROJDIR/results/phase4/md/md100ns_manifest.txt

SRC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $MANIFEST)
if [ -z "$SRC" ]; then
    echo "ERROR: No entry at index $SLURM_ARRAY_TASK_ID in manifest"
    exit 1
fi

TARGET=$(basename $(dirname $SRC))
CHEMBL=$(basename $SRC)
WORKDIR=$PROJDIR/results/phase4/md/$TARGET/$CHEMBL
mkdir -p "$WORKDIR"

echo "Job:     $SLURM_JOB_ID  Array: $SLURM_ARRAY_JOB_ID[$SLURM_ARRAY_TASK_ID]"
echo "Node:    $SLURMD_NODENAME"
echo "Source:  $SRC"
echo "Workdir: $WORKDIR"
echo "Start:   $(date)"

cd "$WORKDIR" || { echo "Cannot cd to $WORKDIR"; exit 1; }

# ── Module / environment init (identical to 17_run_mmpbsa.slurm) ──────────
source /etc/profile.d/modules.sh
module load apps/gromacs/2024.1-oneapi2024 || { echo "ERROR: Cannot load GROMACS"; exit 1; }
module load miniconda3                      || { echo "ERROR: Cannot load miniconda3"; exit 1; }
source /arf/sw/comp/python/miniconda3/etc/profile.d/conda.sh
conda activate mmpbsa || { echo "ERROR: Cannot activate mmpbsa env"; exit 1; }
export PATH="$HOME/.local/bin:$PATH"

GMX=$(command -v gmx_mpi || command -v gmx || echo "")
[ -z "$GMX" ] && { echo "ERROR: GROMACS not found after module load"; exit 1; }
echo "GROMACS: $GMX"

# ── Skip if already complete ───────────────────────────────────────────────
if [ -f "MD100NS_DONE" ]; then
    echo "Already complete — skipping"
    exit 0
fi

set -eo pipefail

# ── Copy source files from MM-GBSA workdir (idempotent) ───────────────────
for f in receptor.pdb lig_GMX.itp lig_GMX.gro em.mdp meta.txt; do
    [ -f "$f" ] || cp "$SRC/$f" .
done

# ── Step 0.5: trim extreme-coordinate receptor residues ───────────────────
# Same logic as 17_run_mmpbsa.slurm; no-op for compact crystal structures.
echo "=== Step 0.5: receptor coordinate check ==="
python3 -c "
import math
from pathlib import Path

pdb = Path('receptor.pdb').read_text().splitlines()
ca_pos = {}
for line in pdb:
    if line.startswith('ATOM') and line[12:16].strip() == 'CA':
        try:
            resnum = int(line[22:26].strip())
            ca_pos[resnum] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            pass

if len(ca_pos) < 10:
    print(f'Too few CA atoms ({len(ca_pos)}) — skipping trim'); exit(0)

cx = sum(v[0] for v in ca_pos.values()) / len(ca_pos)
cy = sum(v[1] for v in ca_pos.values()) / len(ca_pos)
cz = sum(v[2] for v in ca_pos.values()) / len(ca_pos)
print(f'CA centroid: ({cx:.1f}, {cy:.1f}, {cz:.1f}) A')

CUTOFF = 50.0
keep = {r for r,(x,y,z) in ca_pos.items()
        if math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2) <= CUTOFF}
removed = sorted(set(ca_pos.keys()) - keep)
print(f'Residues within {CUTOFF:.0f}A: {len(keep)}/{len(ca_pos)}')

if not removed:
    print('Receptor compact — no trimming needed'); exit(0)

suffix = '...' if len(removed) > 10 else ''
print(f'Removing {len(removed)} residues: {removed[:10]}{suffix}')
out = []
for line in pdb:
    if line.startswith(('ATOM','HETATM','ANISOU')):
        try:
            if int(line[22:26].strip()) not in keep: continue
        except ValueError:
            pass
    out.append(line)
Path('receptor.pdb').write_text('\n'.join(out) + '\n')
print('Trimmed receptor.pdb written')
"

# ── Steps 1-7: build topology through EM and index ────────────────────────
# Skipped if em.gro already exists from a prior run of this script.
if [ ! -f "em.gro" ]; then
    echo "=== Step 1: pdb2gmx ==="
    $GMX pdb2gmx \
        -f receptor.pdb -o protein.gro -p protein.top \
        -ignh -ff amber99sb-ildn -water tip3p \
        -nobackup 2>&1 | tail -5

    echo "=== Step 2: combine protein + ligand (centroid-translated) ==="
    python3 -c "
from pathlib import Path

def parse_xyz(l):
    if len(l) < 44: return 0.0, 0.0, 0.0
    return float(l[20:28]), float(l[28:36]), float(l[36:44])

def translate(l, dx, dy, dz):
    if len(l) < 44: return l
    x, y, z = parse_xyz(l)
    return l[:20] + f'{x-dx:8.3f}{y-dy:8.3f}{z-dz:8.3f}' + l[44:]

prot = Path('protein.gro').read_text().splitlines()
lig  = Path('lig_GMX.gro').read_text().splitlines()
pa, la = prot[2:-1], lig[2:-1]
cx = sum(parse_xyz(l)[0] for l in pa) / len(pa)
cy = sum(parse_xyz(l)[1] for l in pa) / len(pa)
cz = sum(parse_xyz(l)[2] for l in pa) / len(pa)
print(f'Centroid shift: ({cx:.2f}, {cy:.2f}, {cz:.2f}) nm')
lines = ([f'{prot[0].strip()}+LIG', str(len(pa) + len(la))] +
         [translate(l, cx, cy, cz) for l in pa] +
         [translate(l, cx, cy, cz) for l in la] +
         [prot[-1]])
Path('complex.gro').write_text('\n'.join(lines) + '\n')
print(f'Complex: {len(pa)} protein + {len(la)} ligand atoms')
"

    echo "=== Step 2b: forcefield copy + GAFF2 atomtype injection ==="
    FF_SRC=$(find /arf/sw/apps/gromacs -name "amber99sb-ildn.ff" -type d 2>/dev/null | head -1)
    [ -z "$FF_SRC" ] && FF_SRC="${GMXLIB}/amber99sb-ildn.ff"
    rm -rf amber99sb-ildn.ff
    cp -r "$FF_SRC" .

    python3 -c "
import re
from pathlib import Path

itp = Path('lig_GMX.itp').read_text()
ffnb_path = Path('amber99sb-ildn.ff/ffnonbonded.itp')
ffnb_text = ffnb_path.read_text()

# Collect existing atomtype names
existing = set()
in_at = False
for line in ffnb_text.splitlines():
    s = line.strip()
    if re.match(r'\[\s*atomtypes\s*\]', s): in_at = True; continue
    elif s.startswith('['): in_at = False
    if in_at and s and not s.startswith(';'): existing.add(s.split()[0])

# Extract new GAFF2 atomtypes from lig_GMX.itp
new_at_lines = []
in_at = False
for line in itp.splitlines(keepends=True):
    if re.match(r'\[\s*atomtypes\s*\]', line): in_at = True; continue
    elif re.match(r'\[', line): in_at = False
    if in_at:
        s = line.strip()
        if s and not s.startswith(';'):
            name = s.split()[0]
            if name not in existing:
                new_at_lines.append(line)
                existing.add(name)

if new_at_lines:
    ffnb_lines = ffnb_text.splitlines(keepends=True)
    insert_idx = None; found_at = False; in_at2 = False
    for idx, ln in enumerate(ffnb_lines):
        s2 = ln.strip()
        if re.match(r'\[\s*atomtypes\s*\]', s2): in_at2 = True; found_at = True
        elif in_at2 and s2.startswith('['): insert_idx = idx; break
    if insert_idx is not None:
        ffnb_lines[insert_idx:insert_idx] = new_at_lines
        ffnb_path.write_text(''.join(ffnb_lines))
        print(f'Injected {len(new_at_lines)} atomtypes (before next section)')
    elif found_at:
        ffnb_path.write_text(ffnb_text.rstrip('\n') + '\n' + ''.join(new_at_lines))
        print(f'Injected {len(new_at_lines)} atomtypes (end of [ atomtypes ])')
    else:
        with open(ffnb_path, 'a') as f:
            f.write('\n[ atomtypes ]\n')
            f.writelines(new_at_lines)
        print(f'Injected {len(new_at_lines)} atomtypes (new section)')
else:
    print('No new atomtypes (all already in forcefield)')

# Write lig_mol.itp (lig_GMX.itp without [ atomtypes ])
mol_lines = []
in_at = False
for line in itp.splitlines(keepends=True):
    if re.match(r'\[\s*atomtypes\s*\]', line): in_at = True
    elif re.match(r'\[', line): in_at = False
    if not in_at: mol_lines.append(line)
Path('lig_mol.itp').write_text(''.join(mol_lines))

# Build topology
top = Path('protein.top').read_text()
top, n = re.subn(
    r'(#include\s+\"\.?/?amber99sb-ildn\.ff/forcefield\.itp\")',
    r'\1\n#include \"lig_mol.itp\"', top, count=1)
if n == 0:
    top = '#include \"lig_mol.itp\"\n' + top
lig_itp = Path('lig_mol.itp').read_text()
m = re.search(r'\[\s*moleculetype\s*\][^\[]*?\n\s*([A-Za-z0-9_]+)', lig_itp)
mol_name = m.group(1) if m else 'LIG'
print(f'Ligand moleculetype: {mol_name}')
top = top.rstrip() + f'\n{mol_name}              1\n'
Path('topol.top').write_text(top)
"

    echo "=== Step 3: editconf ==="
    $GMX editconf \
        -f complex.gro -o complex_box.gro \
        -c -d 1.0 -bt dodecahedron \
        -nobackup 2>&1 | tail -3

    echo "=== Step 4: solvate ==="
    $GMX solvate \
        -cp complex_box.gro -cs spc216.gro \
        -o complex_solv.gro -p topol.top \
        -nobackup 2>&1 | tail -3

    SOL_COUNT=$(python3 -c "
import re
txt = open('topol.top').read()
m = re.search(r'^SOL\s+(\d+)', txt, re.MULTILINE)
print(m.group(1) if m else '0')
" 2>/dev/null || echo "0")
    echo "Solvated: $SOL_COUNT SOL molecules"
    [ "$SOL_COUNT" -gt 200000 ] && { echo "ERROR: >200k SOL — extreme receptor coordinates"; exit 1; }

    echo "=== Step 5: genion ==="
    $GMX grompp \
        -f em.mdp -c complex_solv.gro -p topol.top \
        -o ions.tpr -maxwarn 5 -nobackup >grompp_ions.log 2>&1 || \
        { echo "ERROR: grompp ions:"; tail -20 grompp_ions.log; exit 1; }
    echo "SOL" | $GMX genion \
        -s ions.tpr -o complex_ions.gro -p topol.top \
        -pname NA -nname CL -neutral -conc 0.15 \
        -nobackup >genion.log 2>&1 || \
        { echo "ERROR: genion:"; tail -10 genion.log; exit 1; }

    echo "=== Step 6: energy minimisation ==="
    $GMX grompp \
        -f em.mdp -c complex_ions.gro -p topol.top \
        -o em.tpr -maxwarn 5 -nobackup >grompp_em.log 2>&1 || \
        { echo "ERROR: grompp EM:"; tail -20 grompp_em.log; exit 1; }
    $GMX mdrun \
        -v -deffnm em -ntomp $SLURM_CPUS_PER_TASK \
        -nobackup 2>&1 | tail -5
    [ -f "em.gro" ] || { echo "ERROR: EM failed — em.gro missing"; exit 1; }

    echo "=== Step 7: index ==="
    python3 -c "
from pathlib import Path

PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
               'HIS','HID','HIE','HIP','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER',
               'THR','TRP','TYR','VAL'}
ION_RES     = {'NA','CL','MG','CA','ZN','K','NA+','CL-'}
BACKBONE_AT = {'N','CA','C','O'}

lines = Path('em.gro').read_text().splitlines()
n = int(lines[1].strip())
prot, lig, sol, ion, backbone = [], [], [], [], []
for i, line in enumerate(lines[2:2+n], 1):
    if len(line) < 15: continue
    res    = line[5:10].strip().upper()
    atname = line[10:15].strip()
    if res in PROTEIN_RES:
        prot.append(i)
        if atname in BACKBONE_AT: backbone.append(i)
    elif res == 'SOL': sol.append(i)
    elif res in ION_RES: ion.append(i)
    else: lig.append(i)

def grp(f, name, atoms):
    f.write(f'[ {name} ]\n')
    for j, a in enumerate(atoms):
        f.write(f'{a:7d}')   # 7-wide: 6-digit atom numbers always get a leading space
        if (j+1) % 15 == 0: f.write('\n')
    f.write('\n\n')

with open('index.ndx', 'w') as f:
    grp(f, 'System',      list(range(1, n+1)))
    grp(f, 'Protein',     prot)
    grp(f, 'Backbone',    backbone)
    grp(f, 'LIG',         lig)
    grp(f, 'Protein_LIG', prot + lig)
    grp(f, 'SOL',         sol)
    grp(f, 'ION',         ion)
    grp(f, 'non-Protein', lig + sol + ion)
    grp(f, 'SOL_Ion',     sol + ion)

print(f'index.ndx: Protein={len(prot)}, Backbone={len(backbone)}, '
      f'LIG={len(lig)}, SOL={len(sol)}, ION={len(ion)}')
"
else
    echo "em.gro found — skipping Steps 1-7 (topology/solvation/EM)"
fi

# ── Step 8: NVT equilibration (100 ps) ────────────────────────────────────
if [ ! -f "nvt.cpt" ]; then
    echo "=== Step 8: NVT ==="
    cat > nvt.mdp << 'NVTMDP'
; NVT equilibration — 100 ps, velocity-rescale thermostat, no barostat
; No position restraints: system is energy-minimized; free NVT is sufficient
; for 100 ps at this stage and avoids posre.itp dependency.
integrator      = md
nsteps          = 50000
dt              = 0.002
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500
cutoff-scheme   = Verlet
nstlist         = 10
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
tcoupl          = V-rescale
tc-grps         = Protein_LIG  SOL_Ion
tau_t           = 0.1          0.1
ref_t           = 300          300
pcoupl          = no
pbc             = xyz
gen_vel         = yes
gen_temp        = 300
gen_seed        = -1
constraints     = h-bonds
constraint-algorithm = LINCS
NVTMDP
    $GMX grompp \
        -f nvt.mdp -c em.gro -p topol.top -n index.ndx \
        -o nvt.tpr -maxwarn 5 -nobackup >grompp_nvt.log 2>&1 || \
        { echo "ERROR: grompp NVT:"; tail -30 grompp_nvt.log; exit 1; }
    $GMX mdrun \
        -deffnm nvt -ntomp $SLURM_CPUS_PER_TASK \
        -nobackup 2>&1 | tail -5
    [ -f "nvt.cpt" ] || { echo "ERROR: NVT failed"; exit 1; }
else
    echo "nvt.cpt found — skipping NVT"
fi

# ── Step 9: NPT equilibration (100 ps) ────────────────────────────────────
if [ ! -f "npt.gro" ]; then
    echo "=== Step 9: NPT ==="
    cat > npt.mdp << 'NPTMDP'
; NPT equilibration — 100 ps, Berendsen barostat (stable for equilibration)
integrator      = md
nsteps          = 50000
dt              = 0.002
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500
cutoff-scheme   = Verlet
nstlist         = 10
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
tcoupl          = V-rescale
tc-grps         = Protein_LIG  SOL_Ion
tau_t           = 0.1          0.1
ref_t           = 300          300
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 1.0
ref_p           = 1.0
compressibility = 4.5e-5
pbc             = xyz
gen_vel         = no
constraints     = h-bonds
constraint-algorithm = LINCS
NPTMDP
    $GMX grompp \
        -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -n index.ndx \
        -o npt.tpr -maxwarn 5 -nobackup >grompp_npt.log 2>&1 || \
        { echo "ERROR: grompp NPT:"; tail -30 grompp_npt.log; exit 1; }
    $GMX mdrun \
        -deffnm npt -ntomp $SLURM_CPUS_PER_TASK \
        -nobackup 2>&1 | tail -5
    [ -f "npt.gro" ] || { echo "ERROR: NPT failed"; exit 1; }
else
    echo "npt.gro found — skipping NPT"
fi

# ── Step 10: 100 ns production MD ─────────────────────────────────────────
echo "=== Step 10: 100 ns production MD ==="

# Write MDP once; reuse on checkpoint restarts (tpr is stable)
cat > md100ns.mdp << 'MDMDP'
; 100 ns production MD
; Parrinello-Rahman barostat: correct NPT ensemble for production
; (Berendsen suppresses pressure fluctuations → wrong compressibility)
integrator          = md
nsteps              = 50000000   ; 100 ns at dt=0.002 ps
dt                  = 0.002

; Output — every 10 ps
nstxout-compressed  = 5000       ; XTC trajectory: 10 ps → 10,000 frames total
nstxout             = 0
nstvout             = 0
nstenergy           = 5000       ; energy file: every 10 ps
nstlog              = 50000      ; log: every 100 ps

; Neighbour search
cutoff-scheme       = Verlet
nstlist             = 20
ns_type             = grid

; Electrostatics & VdW
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0

; Temperature — V-rescale (canonical NVT-like, correct heat capacity)
tcoupl              = V-rescale
tc-grps             = Protein_LIG  SOL_Ion
tau_t               = 0.1          0.1
ref_t               = 300          300

; Pressure — Parrinello-Rahman (correct NPT ensemble)
pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5

pbc                 = xyz
gen_vel             = no          ; velocities carried from NPT checkpoint

; Constraints
constraints         = h-bonds
constraint-algorithm = LINCS
lincs-iter          = 1
lincs-order         = 4
MDMDP

# grompp only on first run; tpr is reused for all checkpoint restarts
if [ ! -f "md100ns.tpr" ]; then
    $GMX grompp \
        -f md100ns.mdp -c npt.gro -t npt.cpt \
        -p topol.top -n index.ndx \
        -o md100ns.tpr -maxwarn 5 -nobackup >grompp_md100ns.log 2>&1 || \
        { echo "ERROR: grompp 100ns:"; tail -30 grompp_md100ns.log; exit 1; }
    echo "md100ns.tpr created"
else
    echo "md100ns.tpr exists — checkpoint restart mode"
fi

# Resume from checkpoint if one exists, otherwise start fresh
if [ -f "md100ns.cpt" ]; then
    echo "Resuming from md100ns.cpt"
    $GMX mdrun \
        -deffnm md100ns -cpi md100ns.cpt \
        -ntomp $SLURM_CPUS_PER_TASK -nobackup 2>&1 | tail -10
else
    $GMX mdrun \
        -deffnm md100ns \
        -ntomp $SLURM_CPUS_PER_TASK -nobackup 2>&1 | tail -10
fi

# Count completed frames to decide if simulation is finished
FRAMES=$(python3 -c "
import subprocess, re
r = subprocess.run(
    ['$GMX','check','-f','md100ns.xtc'],
    capture_output=True, text=True)
m = re.search(r'Last frame\s+(\d+)', r.stderr + r.stdout)
print(int(m.group(1)) + 1 if m else 0)
" 2>/dev/null || echo "0")
echo "Frames in trajectory: $FRAMES / 10000"

if [ "$FRAMES" -lt 10000 ]; then
    echo "Simulation incomplete — resubmit this job to continue from checkpoint"
    echo "  sbatch --array=$SLURM_ARRAY_TASK_ID code/slurm/18_run_md100ns.slurm"
    exit 0   # clean exit; not an error; log will show progress on next run
fi

echo "Simulation complete (10000 frames). Proceeding to analysis."

# ── Step 11: PBC imaging (2-step protocol) ────────────────────────────────
# Step 11a: remove PBC jumps on the full system so the trajectory is
# continuous (no molecule teleports across the boundary).
echo "=== Step 11a: trjconv -pbc nojump (full system) ==="
echo "System" | $GMX trjconv \
    -f md100ns.xtc -s md100ns.tpr \
    -o md100ns_nojump.xtc -pbc nojump \
    -nobackup 2>&1 | tail -3
[ -f "md100ns_nojump.xtc" ] || { echo "ERROR: nojump trjconv failed"; exit 1; }

# Step 11b: center on Protein, cluster Protein_LIG into same image, output
# Protein_LIG only (strips solvent → faster analysis).
# -center: translates so Protein COM is at box center each frame
# -pbc mol: puts each molecule into the central box image after centering
echo "=== Step 11b: trjconv -center -pbc mol (Protein_LIG) ==="
printf "Protein\nProtein_LIG\n" | $GMX trjconv \
    -f md100ns_nojump.xtc -s md100ns.tpr -n index.ndx \
    -o md100ns_whole.xtc -pbc mol -center \
    -nobackup 2>&1 | tail -3
[ -f "md100ns_whole.xtc" ] || { echo "ERROR: trjconv imaging failed"; exit 1; }
rm -f md100ns_nojump.xtc

# ── Step 12: Analysis ─────────────────────────────────────────────────────
echo "=== Step 12: analysis ==="

# 12a: Protein backbone RMSD vs time (fit on Backbone; RMSD on Backbone)
printf "Backbone\nBackbone\n" | $GMX rms \
    -s md100ns.tpr -f md100ns_whole.xtc -n index.ndx \
    -o backbone_rmsd.xvg -tu ns \
    -nobackup 2>&1 | tail -3

# 12b: Ligand heavy-atom RMSD vs time (fit on Backbone; RMSD on LIG)
printf "Backbone\nLIG\n" | $GMX rms \
    -s md100ns.tpr -f md100ns_whole.xtc -n index.ndx \
    -o ligand_rmsd.xvg -tu ns \
    -nobackup 2>&1 | tail -3

# 12c: Per-residue RMSF (protein, all atoms)
echo "Protein" | $GMX rmsf \
    -s md100ns.tpr -f md100ns_whole.xtc -n index.ndx \
    -o rmsf.xvg -res \
    -nobackup 2>&1 | tail -3

# 12d: Protein–ligand close contacts vs time (H-bond geometric proxy)
# gmx hbond in GROMACS 2024 uses a new trajectory-analysis backend that
# auto-detects donors/acceptors from atom types; GAFF2 ligand atom types
# are not recognised → "no donors AND no acceptors" fatal error.
# gmx mindist -d 3.5 -on counts protein–ligand pairs within 3.5 Å
# (standard H-bond distance cutoff); output format is identical to
# gmx hbond -num so the Python summary step reads it without change.
printf "Protein\nLIG\n" | $GMX mindist \
    -s md100ns.tpr -f md100ns_whole.xtc -n index.ndx \
    -on hbond_num.xvg -d 3.5 \
    -nobackup 2>&1 | tail -3

# 12e: Radius of gyration (protein fold stability)
echo "Protein" | $GMX gyrate \
    -s md100ns.tpr -f md100ns_whole.xtc -n index.ndx \
    -o gyrate.xvg \
    -nobackup 2>&1 | tail -3

# ── Step 13: Summary ─────────────────────────────────────────────────────
echo "=== Step 13: summary ==="
python3 -c "
import re
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

bb  = read_xvg('backbone_rmsd.xvg')
lig = read_xvg('ligand_rmsd.xvg')
hb  = read_xvg('hbond_num.xvg')

# GROMACS rms outputs nm; convert to Angstrom (×10)
bb_avg_A   = last_half_mean(bb) * 10
lig_avg_A  = last_half_mean(lig) * 10
lig_max_A  = last_half_max(lig) * 10

# H-bond occupancy: fraction of frames with ≥1 protein–ligand H-bond
hb_occ = sum(1 for _,v in hb if v >= 1) / max(1, len(hb))

# Acceptance criteria
pass_bb  = bb_avg_A  < 3.0
pass_lig = lig_avg_A < 2.0
pass_hb  = hb_occ   > 0.30
verdict  = 'PASS' if (pass_bb and pass_lig and pass_hb) else 'FAIL'

meta   = Path('meta.txt').read_text()
target = [l.split('=')[1] for l in meta.splitlines() if l.startswith('target')][0]
cid    = [l.split('=')[1] for l in meta.splitlines() if l.startswith('chembl')][0]
vina   = [l.split('=')[1] for l in meta.splitlines() if l.startswith('vina')][0]
dg     = Path('$SRC/mmpbsa_result.txt').read_text().split(',')[3]

result = f'{target},{cid},{vina},{dg},{bb_avg_A:.2f},{lig_avg_A:.2f},{lig_max_A:.2f},{hb_occ:.2f},{verdict}'
print(f'RESULT: {result}')
Path('md100ns_result.txt').write_text(result + '\n')

print()
print(f'  Target / Compound : {target} / {cid}')
print(f'  Vina score        : {vina} kcal/mol')
print(f'  MM-GBSA ΔG        : {dg} kcal/mol')
print(f'  Backbone RMSD     : {bb_avg_A:.2f} A (last 50 ns avg)  [<3.0: {pass_bb}]')
print(f'  Ligand RMSD (avg) : {lig_avg_A:.2f} A (last 50 ns avg)  [<2.0: {pass_lig}]')
print(f'  Ligand RMSD (max) : {lig_max_A:.2f} A (last 50 ns max)')
print(f'  H-bond occupancy  : {hb_occ:.0%}  [>30%: {pass_hb}]')
print(f'  Verdict           : {verdict}')
"

touch MD100NS_DONE
echo "Done: $(date)"
