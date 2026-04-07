"""
DAM-DRUG Phase 4 — MM-GBSA preparation
========================================
For each compound in candidates_shortlist.csv:
  1. Extract model-1 Vina pose (PDBQT → PDB)
  2. Fetch SMILES from ChEMBL → generate GAFF2 topology via acpype
  3. Copy receptor prepared PDB
  4. Write GROMACS MDP files and topology template
  5. Append to manifest.txt consumed by 17_run_mmpbsa.slurm

Run on TRUBA login node (no GPU needed — prep only).

Usage:
    python 17_prep_mmpbsa.py

Dependencies (must be available):
    obabel (~/.local/bin/obabel)
    acpype   (pip install acpype)
    rdkit    (already installed)
    requests (stdlib)
"""

import csv
import json
import logging
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT      = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
VINA_OUT     = PROJECT / "results/phase4/vina"
MMPBSA_DIR   = PROJECT / "results/phase4/mmpbsa"
RECEPTOR_DIR = PROJECT / "data/structures/prepared"
SHORTLIST    = PROJECT / "results/phase4/candidates_shortlist.csv"
MANIFEST     = MMPBSA_DIR / "manifest.txt"
OBABEL       = os.environ.get("OBABEL", str(Path.home() / ".local/bin/obabel"))
ACPYPE       = os.environ.get("ACPYPE", str(Path.home() / "miniconda3/envs/mmpbsa/bin/acpype"))
SDF_DIR      = PROJECT / "data/compounds/sdf_3d"

# receptor stem → prepared PDB filename
RECEPTOR_PDB = {
    "PPARG_1FM9_LBD_prep":   "PPARG_1FM9_LBD_prep.pdb",
    "IRF8_AF2_DBD_prep":     "IRF8_AF2_DBD_prep.pdb",
    "MAF_4EOT_bZIP_prep":    "MAF_4EOT_bZIP_prep.pdb",
    "RUNX1_1LJM_Runt_prep":  "RUNX1_1LJM_Runt_prep.pdb",
    "BHLHE41_AF2_bHLH_prep": "BHLHE41_AF2_bHLH_prep.pdb",
}


# ── GROMACS MDP file templates ─────────────────────────────────────────────

EM_MDP = """\
; Energy minimisation
integrator  = steep
emtol       = 500.0
emstep      = 0.01
nsteps      = 5000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""

NVT_MDP = """\
; NVT equilibration (100 ps, position restraints)
define      = -DPOSRES
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
pcoupl      = no
pbc         = xyz
gen_vel     = yes
gen_temp    = 300
gen_seed    = 42
"""

NPT_MDP = """\
; NPT equilibration (100 ps, position restraints)
define      = -DPOSRES
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
"""

MD_MDP = """\
; Production MD (1 ns)
integrator  = md
nsteps      = 500000
dt          = 0.002
nstxout     = 0
nstvout     = 0
nstfout     = 0
nstxout-compressed = 5000
compressed-x-grps  = Protein LIG
nstenergy   = 5000
nstlog      = 5000
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
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
gen_vel     = no
"""

MMPBSA_IN = """\
&general
   startframe=251, endframe=500, interval=1,
   verbose=2, keep_files=0,
/
&gb
   igb=2, saltcon=0.15,
/
"""

TOPOL_TEMPLATE = """\
; GROMACS topology for {target} + {chembl_id}
#include "amber99sb-ildn.ff/forcefield.itp"
#include "lig_GMX.itp"

[ system ]
{target}_{chembl_id}

[ molecules ]
Protein_chain_A  1
LIG              1
"""


def fetch_smiles(chembl_id: str) -> str:
    """Fetch canonical SMILES from ChEMBL REST API."""
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            data = json.loads(r.read())
        smiles = data["molecule_structures"]["canonical_smiles"]
        return smiles
    except Exception as e:
        log.warning(f"  ChEMBL API failed for {chembl_id}: {e}")
        return ""


def get_formal_charge(smiles: str) -> int:
    """Compute formal charge from SMILES using RDKit."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        return sum(a.GetFormalCharge() for a in mol.GetAtoms())
    except Exception:
        return 0


def extract_model1_pdb(pdbqt_path: Path, out_pdb: Path) -> bool:
    """Extract model 1 from Vina PDBQT and convert to PDB via obabel."""
    # Strip MODEL/ENDMDL/REMARK to get model 1
    lines = []
    in_model = False
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("MODEL 1"):
            in_model = True
            continue
        if line.startswith("ENDMDL") and in_model:
            break
        if in_model and not line.startswith("REMARK"):
            lines.append(line)
    if not lines:
        return False
    tmp_pdbqt = out_pdb.parent / "_tmp_model1.pdbqt"
    tmp_pdbqt.write_text("\n".join(lines) + "\n")
    result = subprocess.run(
        [OBABEL, str(tmp_pdbqt), "-O", str(out_pdb), "-h"],
        capture_output=True, text=True
    )
    tmp_pdbqt.unlink(missing_ok=True)
    return out_pdb.exists() and out_pdb.stat().st_size > 0


def run_acpype(sdf_path: Path, work_dir: Path, charge: int) -> bool:
    """Run acpype to generate GAFF2 GROMACS topology for the ligand.

    Uses AM1-BCC charges via tleap from conda env 'mmpbsa'
    (module load miniconda3 && conda activate mmpbsa).
    """
    result = subprocess.run(
        [ACPYPE, "-i", str(sdf_path), "-c", "bcc", "-a", "gaff2",
         "-n", str(charge), "-b", "lig", "-o", "gmx"],
        cwd=str(work_dir),
        capture_output=True, text=True, timeout=300
    )
    # acpype creates <basename>.acpype/ directory
    acpype_dir = work_dir / "lig.acpype"
    itp = acpype_dir / "lig_GMX.itp"
    gro = acpype_dir / "lig_GMX.gro"
    if itp.exists() and gro.exists():
        shutil.copy(itp, work_dir / "lig_GMX.itp")
        shutil.copy(gro, work_dir / "lig_GMX.gro")
        # Also copy posre file if exists (acpype names it posre_lig.itp)
        for posre_name in ("posre_lig.itp", "lig_GMX_posre.itp"):
            posre = acpype_dir / posre_name
            if posre.exists():
                shutil.copy(posre, work_dir / "lig_GMX_posre.itp")
                break
        return True
    log.warning(f"  acpype output missing.\n  stdout: {result.stdout[-300:]}\n  stderr: {result.stderr[-300:]}")
    return False


def prep_complex(target: str, chembl_id: str, vina_score: float) -> bool:
    work_dir = MMPBSA_DIR / target / chembl_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Extract Vina pose ────────────────────────────────────────────────
    vina_pdbqt = VINA_OUT / target / f"{chembl_id}.pdbqt"
    ligand_pdb = work_dir / "ligand_vina.pdb"
    if not vina_pdbqt.exists():
        log.error(f"  Vina PDBQT not found: {vina_pdbqt}")
        return False
    if not extract_model1_pdb(vina_pdbqt, ligand_pdb):
        log.error(f"  Failed to extract model 1 from {vina_pdbqt}")
        return False
    log.info(f"  Pose extracted → {ligand_pdb.name}")

    # ── 2. Get SMILES & generate SDF ────────────────────────────────────────
    # Try cached SDF from ligand prep stage first
    cached_sdf = SDF_DIR / f"{chembl_id}.sdf"
    ligand_sdf = work_dir / "ligand.sdf"
    if cached_sdf.exists():
        shutil.copy(cached_sdf, ligand_sdf)
        log.info(f"  SDF from cache: {cached_sdf.name}")
    else:
        smiles = fetch_smiles(chembl_id)
        if not smiles:
            log.error(f"  No SMILES for {chembl_id}")
            return False
        (work_dir / "ligand.smi").write_text(f"{smiles}\t{chembl_id}\n")
        result = subprocess.run(
            [OBABEL, "-ismi", str(work_dir / "ligand.smi"),
             "-O", str(ligand_sdf), "--gen3d", "-p", "7.4"],
            capture_output=True, text=True
        )
        if not ligand_sdf.exists():
            log.error(f"  obabel SDF generation failed: {result.stderr[:100]}")
            return False
        log.info(f"  SDF generated from SMILES")
        time.sleep(0.3)  # be polite to ChEMBL API

    # Get formal charge for acpype
    smiles_for_charge = ""
    if (work_dir / "ligand.smi").exists():
        smiles_for_charge = (work_dir / "ligand.smi").read_text().split()[0]
    charge = get_formal_charge(smiles_for_charge)

    # ── 3. Run acpype (GAFF2 topology) ──────────────────────────────────────
    if not (work_dir / "lig_GMX.itp").exists():
        log.info(f"  Running acpype (charge={charge})...")
        if not run_acpype(ligand_sdf, work_dir, charge):
            log.error(f"  acpype failed for {chembl_id}")
            return False
        log.info(f"  GAFF2 topology ready")
    else:
        log.info(f"  GAFF2 topology already exists, skipping acpype")

    # ── 4. Copy receptor PDB ─────────────────────────────────────────────────
    rec_pdb_name = RECEPTOR_PDB.get(target)
    if not rec_pdb_name:
        log.error(f"  Unknown target: {target}")
        return False
    rec_src = RECEPTOR_DIR / rec_pdb_name
    rec_dst = work_dir / "receptor.pdb"
    if not rec_src.exists():
        log.error(f"  Receptor PDB not found: {rec_src}")
        return False
    shutil.copy(rec_src, rec_dst)

    # ── 5. Write MDP files ───────────────────────────────────────────────────
    (work_dir / "em.mdp").write_text(EM_MDP)
    (work_dir / "nvt.mdp").write_text(NVT_MDP)
    (work_dir / "npt.mdp").write_text(NPT_MDP)
    (work_dir / "md.mdp").write_text(MD_MDP)
    (work_dir / "mmpbsa.in").write_text(MMPBSA_IN)

    # ── 6. Write topology template ───────────────────────────────────────────
    (work_dir / "topol_template.top").write_text(
        TOPOL_TEMPLATE.format(target=target, chembl_id=chembl_id)
    )

    # ── 7. Write metadata ────────────────────────────────────────────────────
    (work_dir / "meta.txt").write_text(
        f"target={target}\nchembl_id={chembl_id}\nvina_score={vina_score}\n"
    )

    return True


def main():
    MMPBSA_DIR.mkdir(parents=True, exist_ok=True)

    shortlist = list(csv.DictReader(open(SHORTLIST)))
    log.info(f"Preparing {len(shortlist)} complexes for MM-GBSA")

    manifest_lines = []
    n_ok = n_fail = 0

    for row in shortlist:
        target    = row["target"]
        chembl_id = row["chembl_id"]
        vina      = float(row["vina_score"])
        log.info(f"\n── {target} / {chembl_id}  (Vina={vina:.2f}) ──")

        ok = prep_complex(target, chembl_id, vina)
        if ok:
            work_dir = MMPBSA_DIR / target / chembl_id
            manifest_lines.append(str(work_dir))
            n_ok += 1
            log.info(f"  ✓ Ready")
        else:
            n_fail += 1
            log.warning(f"  ✗ Failed — will be skipped in SLURM array")

    MANIFEST.write_text("\n".join(manifest_lines) + "\n")

    log.info(f"\n{'='*60}")
    log.info(f"Prepared: {n_ok}   Failed: {n_fail}")
    log.info(f"Manifest: {MANIFEST}  ({n_ok} entries)")
    log.info(f"Next: sbatch code/slurm/17_run_mmpbsa.slurm  (array 0-{n_ok-1})")


if __name__ == "__main__":
    main()
