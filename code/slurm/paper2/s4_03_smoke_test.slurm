#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_tremd_smoke
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=akya-cuda
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=10   # akya-cuda enforces exactly 10 cores/GPU (40
#                               # cores, 4 GPUs per node); plain --ntasks or
#                               # --cpus-per-task is rejected, only
#                               # --ntasks-per-node satisfies the GRES check
#                               # (verified 2026-09-11). mpirun -np 8 below
#                               # uses 8 of the 10 allocated slots.
#SBATCH --mem=48G
#SBATCH --array=0
# ─────────────────────────────────────────────────────────────────────────────
# WP3.3-smoke — 10 ps T-REMD dry run for ONE complex (GLUE0231, array=0).
# Verifies, before committing GPU allocation to the real 100 ns x 8-replica
# production run (s4_03_run_tremd.slurm):
#   (a) the harmonic Zn2+ restraints survive grompp with no new warnings
#       (already confirmed separately: s4_02 ran EM/NVT/500 ps NPT clean with
#       them active) and don't blow up once replica exchange starts perturbing
#       coordinates across 8 independent replicas;
#   (b) the -replex MPI exchange handshake completes multiple attempts across
#       8 GPU-offloaded ranks without an interconnect deadlock.
# Writes to smoke_rep0..7/ (separate from the real run's rep0..7/) so this
# never collides with or gets mistaken for production output.
#
# Requires: s4_02_gromacs_prep.slurm already completed for this complex
# (reads npt.gro / npt.cpt / topol.top from the array-0 working directory).
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail
PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
source /etc/profile.d/modules.sh
module load apps/gromacs/2023.2-cuda || { echo "ERROR: cannot load GROMACS"; exit 1; }
GMX="gmx_mpi"

TOPCSV=$PROJDIR/results/phase7/glue_design/glue_candidates_top.csv
GID=$(awk -F, -v n=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==n{print $2}' "$TOPCSV")
WD=$PROJDIR/results/phase7/glue_md/$GID
cd "$WD"
[ -s npt.gro ] || { echo "ERROR: npt.gro missing — run s4_02_gromacs_prep.slurm first"; exit 1; }
echo "[$GID] T-REMD SMOKE TEST — start $(date)"

TEMPS=(300.00 302.79 305.60 308.44 311.30 314.19 317.10 320.04)
NREP=${#TEMPS[@]}

for i in $(seq 0 $((NREP-1))); do
    mkdir -p "smoke_rep$i"
    sed "s/__REF_T__/${TEMPS[$i]}/; s/__GENSEED__/$((1000+i))/" \
        > "smoke_rep$i/prod.mdp" <<'EOF'
integrator  = md
nsteps      = 5000             ; 10 ps SMOKE TEST @ 2 fs -- NOT production
dt          = 0.002
nstxout-compressed = 100
nstenergy   = 100
nstlog      = 100
continuation = yes
constraints = h-bonds
constraint-algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
DispCorr    = EnerPres
tcoupl      = V-rescale
tc-grps     = Protein_LIG Water_and_ions
tau-t       = 0.1 0.1
ref-t       = __REF_T__ __REF_T__
pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
refcoord-scaling = com
gen-vel     = no
gen-seed    = __GENSEED__
EOF
    [ -s smoke_index.ndx ] || $GMX make_ndx -f npt.gro -o smoke_index.ndx <<< $'"Protein" | "Other"\nname 20 Protein_LIG\n"Water" | "Ion"\nname 21 Water_and_ions\nq'
    $GMX grompp -f "smoke_rep$i/prod.mdp" -c npt.gro -t npt.cpt -p topol.top \
        -n smoke_index.ndx -o "smoke_rep$i/prod.tpr" -maxwarn 3
done

echo "[$GID] grompp x8 OK — starting replica exchange (attempt every 100 steps = 200 fs)"
mpirun -np $NREP $GMX mdrun -v -deffnm prod \
    -multidir smoke_rep0 smoke_rep1 smoke_rep2 smoke_rep3 smoke_rep4 smoke_rep5 smoke_rep6 smoke_rep7 \
    -replex 100 -nb gpu -bonded gpu -pme gpu

echo "[$GID] SMOKE TEST done $(date)"
echo "Check smoke_rep0/prod.log for 'Repl' exchange-attempt lines and no LINCS/segfault."
