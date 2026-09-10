#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_gmx_prep
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --array=0-2
# ─────────────────────────────────────────────────────────────────────────────
# WP3.2 — GROMACS system prep for one ternary complex (array index = glue rank).
# Amber14SB (protein) + GAFF2 (ligand, from acpype) + TIP3P, dodecahedron box,
# 150 mM NaCl neutralising. EM -> NVT -> NPT. Zn kept as ion with position
# restraints during equilibration.
#
# GROMACS: TRUBA system module apps/gromacs/2024.1-oneapi2024 (as in the
# accepted paper's MD).
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
module load apps/gromacs/2024.1-oneapi2024
GMX="gmx_mpi"

TOPCSV=$PROJDIR/results/phase7/glue_design/glue_candidates_top.csv
GID=$(awk -F, -v n=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==n{print $2}' "$TOPCSV")   # col 2 = glue_id, +1 header
WD=$PROJDIR/results/phase7/glue_md/$GID
cd "$WD"
echo "[$GID] prep — start $(date)"

ACP=$WD/${GID}.acpype
LIG_ITP=$ACP/${GID}_GMX.itp
LIG_GRO=$ACP/${GID}_GMX.gro
[ -s "$LIG_ITP" ] || { echo "ERROR: run s4_01_build_ternary.py first"; exit 1; }

# ── 1. protein topology (Amber14SB, TIP3P) ─────────────────────────────────
mkdir -p mdp
$GMX pdb2gmx -f complex.pdb -o proc.gro -p topol.top -i posre.itp \
    -ff amber14sb -water tip3p -ignh <<< $'1\n1'   # adjust chain termini as needed

# ── 2. splice ligand into topol.top ───────────────────────────────────────
python - << PY
from pathlib import Path
top = Path("topol.top"); lines = top.read_text().splitlines()
gid = "$GID"
# insert #include for the ligand itp after the forcefield include
out = []
for ln in lines:
    out.append(ln)
    if ln.strip().startswith('#include') and 'forcefield.itp' in ln:
        out.append(f'#include "{gid}.acpype/{gid}_GMX.itp"')
# add ligand to [ molecules ]
txt = "\n".join(out)
if "\n[ molecules ]" in txt:
    txt = txt.rstrip() + f"\n{gid}                 1\n"
top.write_text(txt)
print("ligand spliced into topol.top")
PY

# combine ligand coords into the processed structure
python - << PY
# concatenate proc.gro (protein) + ligand gro into complex_box-ready gro
from pathlib import Path
p = Path("proc.gro").read_text().splitlines()
l = Path("$LIG_GRO").read_text().splitlines()
np_ = int(p[1]); nl = int(l[1])
box = p[-1]
atoms = p[2:2+np_] + l[2:2+nl]
Path("complex_full.gro").write_text(
    p[0] + "\n" + str(np_+nl) + "\n" + "\n".join(atoms) + "\n" + box + "\n")
print("complex_full.gro:", np_+nl, "atoms")
PY

# ── 3. box + solvate + ions ───────────────────────────────────────────────
$GMX editconf -f complex_full.gro -o box.gro -bt dodecahedron -d 1.2 -c
$GMX solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

cat > mdp/ions.mdp <<'EOF'
integrator = steep
nsteps     = 5000
emtol      = 1000
EOF
$GMX grompp -f mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | $GMX genion -s ions.tpr -o ionised.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15

# ── 3b. Zn2+ harmonic distance restraints ────────────────────────────────
# Added AFTER all [ molecules ] edits (solvate/genion) so the
# [ intermolecular_interactions ] block stays last in topol.top. Global
# 1-based indices come from complex_full.gro (solute only); solvent/ions are
# appended after the solute, so those indices are unchanged in the final
# system. k = 10000 kJ/mol/nm^2; r0 0.230 nm Zn-S(Cys), 0.200 nm Zn-N(His).
python "$PROJDIR/code/phase7_robustness/s4_md/s4_zn_restraints.py" \
    complex_full.gro topol.top
grep -q 'intermolecular_interactions' topol.top || \
    { echo "ERROR: Zn restraints not written to topol.top"; exit 1; }
# EM/NVT/NPT/production below all read the updated topol.top; the
# [ intermolecular_interactions ] harmonic bonds are always active (not gated
# by -DPOSRES), so Zn stays caged through equilibration and T-REMD.

# ── 4. EM ─────────────────────────────────────────────────────────────────
cat > mdp/em.mdp <<'EOF'
integrator = steep
nsteps     = 50000
emtol      = 500
emstep     = 0.01
nstlist    = 10
cutoff-scheme = Verlet
rlist      = 1.0
rcoulomb   = 1.0
rvdw       = 1.0
coulombtype = PME
EOF
$GMX grompp -f mdp/em.mdp -c ionised.gro -p topol.top -o em.tpr -maxwarn 2
$GMX mdrun -deffnm em -ntomp $SLURM_CPUS_PER_TASK

# ── 5. NVT (300 K, restrain protein heavy atoms + Zn) ──────────────────────
cat > mdp/nvt.mdp <<'EOF'
integrator  = md
nsteps      = 100000          ; 200 ps
dt          = 0.002
tcoupl      = V-rescale
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
pcoupl      = no
gen-vel     = yes
gen-temp    = 300
gen-seed    = -1
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
define      = -DPOSRES
EOF
$GMX grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
$GMX mdrun -deffnm nvt -ntomp $SLURM_CPUS_PER_TASK

# ── 6. NPT (300 K, 1 bar) ────────────────────────────────────────────────
cat > mdp/npt.mdp <<'EOF'
integrator  = md
nsteps      = 250000          ; 500 ps
dt          = 0.002
tcoupl      = V-rescale
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
define      = -DPOSRES
EOF
$GMX grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
$GMX mdrun -deffnm npt -ntomp $SLURM_CPUS_PER_TASK

echo "[$GID] equilibrated -> npt.gro / npt.cpt  ($(date))"
echo "Next: sbatch code/slurm/paper2/s4_03_run_tremd.slurm"
