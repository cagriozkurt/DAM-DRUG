#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_tremd
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --open-mode=append
#SBATCH --time=72:00:00
#SBATCH --partition=akya-cuda
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=48G
#SBATCH --array=0-2
# ─────────────────────────────────────────────────────────────────────────────
# WP3.3 — Temperature Replica-Exchange MD (T-REMD) for one ternary complex.
# 8 replicas, 300-320 K (geometric ladder), GROMACS -multidir -replex.
# 100 ns per replica; exchange attempt every 1 ps.
# Advancement gate (checked in s4_04): mean ligand core-RMSD < 3.5 A over the
# final 20 ns across ALL replicas.
#
# barbun/akya CPUs ~1-2 ns/day; on a V100/A100 GPU ~20-50 ns/day for ~130k atoms
# -> 100 ns in ~2-5 days. Resubmit the same sbatch; mdrun resumes from .cpt.
#
# Two bugs found running the RunPod smoke test (10ps dry run, 2026-09-11)
# that apply equally to this production script -- fixed here before first
# launch (job 6346140 cancelled pre-start, never actually ran):
#  1. `sed "s/__REF_T__/.../"` (no /g) only replaced the FIRST of the two
#     __REF_T__ placeholders on the `ref-t = __REF_T__ __REF_T__` line ->
#     grompp ERROR "Invalid value for mdp option ref-t" (real number pair
#     required, one side left as literal text).
#  2. gmx's generic "Ion" index family includes Zn2+ alongside Na/Cl, so
#     "Water" | "Ion" pulls Zn into the water/ion thermostat group WHILE
#     "Protein" | "Other" already contains it too (Zn sits contiguously
#     between the protein and ligand blocks in the coordinate file) ->
#     grompp FATAL "Atom in multiple T-Coupling groups". Fixed by merging
#     from "NA"/"CL" (not "Ion") for the solvent group, explicitly folding
#     "ZN" into the solute group, and using gmx's own auto-generated merge
#     names (Protein_Other_ZN / Water_NA_CL) directly in tc-grps instead of
#     hardcoded `name 20/21` (GROMACS 2024+ injects extra default index
#     groups at 20/21, silently renaming the WRONG pre-existing groups).
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
source /etc/profile.d/modules.sh
module load apps/gromacs/2023.2-cuda || { echo "ERROR: cannot load GROMACS"; exit 1; }  # GPU build for T-REMD
GMX="gmx_mpi"

TOPCSV=$PROJDIR/results/phase7/glue_design/glue_candidates_top.csv
GID=$(awk -F, -v n=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==n{print $2}' "$TOPCSV")
WD=$PROJDIR/results/phase7/glue_md/$GID
cd "$WD"
echo "[$GID] T-REMD — start $(date)"

# 8-replica geometric temperature ladder 300 -> 320 K
TEMPS=(300.00 302.79 305.60 308.44 311.30 314.19 317.10 320.04)
NREP=${#TEMPS[@]}

# ── build per-replica dirs + tpr on first run ─────────────────────────────
if [ ! -f rep0/prod.tpr ]; then
  for i in $(seq 0 $((NREP-1))); do
    mkdir -p "rep$i"
    sed "s/__REF_T__/${TEMPS[$i]}/g; s/__GENSEED__/$((1000+i))/g" \
        > "rep$i/prod.mdp" <<'EOF'
integrator  = md
nsteps      = 50000000        ; 100 ns @ 2 fs
dt          = 0.002
nstxout-compressed = 5000     ; 10 ps
nstenergy   = 5000
nstlog      = 5000
continuation = yes
constraints = h-bonds
constraint-algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
DispCorr    = EnerPres
tcoupl      = V-rescale
tc-grps     = Protein_Other_ZN Water_NA_CL
tau-t       = 0.1 0.1
ref-t       = __REF_T__ __REF_T__
pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
gen-vel     = no
gen-seed    = __GENSEED__
EOF
    $GMX make_ndx -f npt.gro -o "rep$i/index.ndx" <<< $'"Protein" | "Other" | "ZN"\n"Water" | "NA" | "CL"\nq'
    $GMX grompp -f "rep$i/prod.mdp" -c npt.gro -t npt.cpt -p topol.top \
        -n "rep$i/index.ndx" -o "rep$i/prod.tpr" -maxwarn 3
  done
fi

# ── run with replica exchange (attempt every 1000 steps = 2 ps) ───────────
mpirun -np $NREP $GMX mdrun -v -deffnm prod \
    -multidir rep0 rep1 rep2 rep3 rep4 rep5 rep6 rep7 \
    -replex 1000 -nb gpu -bonded gpu -pme gpu -cpi prod.cpt

echo "[$GID] T-REMD segment done $(date). Resubmit until 100 ns reached."
