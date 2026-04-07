"""
DAM-DRUG Phase 4 — MM-GBSA preparation for Tier-2 hits
========================================================
Prepares work directories for the 10 selected Tier-2 consensus hits
(6 vs IKZF1_8RQC, 4 vs IRF8_AF2_DBD) for gmx_MMPBSA via 17_run_mmpbsa.slurm.

Run on TRUBA login node (no GPU needed — prep only):
    conda activate mmpbsa
    python code/phase4_docking/22_prep_mmpbsa_tier2.py

Requires:
    obabel (~/.local/bin/obabel)
    acpype (conda env mmpbsa)
    rdkit  (conda env mmpbsa)
"""

import json
import logging
import os
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT      = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
VINA_TIER2   = PROJECT / "results/phase4/vina_tier2"
MMPBSA_DIR   = PROJECT / "results/phase4/mmpbsa_tier2"
RECEPTOR_DIR = PROJECT / "data/structures/prepared"
MANIFEST     = MMPBSA_DIR / "manifest.txt"
OBABEL       = os.environ.get("OBABEL", str(Path.home() / ".local/bin/obabel"))
ACPYPE       = os.environ.get("ACPYPE", "acpype")   # in conda env PATH

# Tier-2 target name → receptor PDB filename
RECEPTOR_PDB = {
    "IRF8_AF2_DBD": "IRF8_AF2_DBD_prep.pdb",
    "IKZF1_8RQC":   "8RQC_CRBN_ZF2_prep.pdb",
}

# Selected Tier-2 compounds for MM-GBSA
# (target, chembl_id, vina_score, rationale)
SELECTED = [
    # ── IKZF1_8RQC ──────────────────────────────────────────────────────────
    ("IKZF1_8RQC", "CHEMBL4297528", -7.28,  "RISDIPLAM: SMN2 splicing/RNA processing"),
    ("IKZF1_8RQC", "CHEMBL3137320", -6.931, "TALAZOPARIB: PARP inhibitor, neuroprotection"),
    ("IKZF1_8RQC", "CHEMBL118",     -6.286, "CELECOXIB: COX-2/neuroinflammation"),
    ("IKZF1_8RQC", "CHEMBL1096",    -6.334, "AMLEXANOX: TBK1/IKKe inhibitor, innate immune"),
    ("IKZF1_8RQC", "CHEMBL1229211", -6.496, "DOLUTEGRAVIR: CNS-penetrant, highest CNN=0.998"),
    ("IKZF1_8RQC", "CHEMBL973",     -6.029, "TERIFLUNOMIDE: MS drug, CNS immune modulation"),
    # ── IRF8_AF2_DBD ─────────────────────────────────────────────────────────
    ("IRF8_AF2_DBD", "CHEMBL3261331", -7.276, "RESMETIROM: THR-beta agonist, neuroinflammation"),
    ("IRF8_AF2_DBD", "CHEMBL3188267", -7.23,  "CAPMATINIB: MET inhibitor, microglia signaling"),
    ("IRF8_AF2_DBD", "CHEMBL3039520", -7.118, "LASMIDITAN: 5-HT1F agonist, CNS-selective"),
    ("IRF8_AF2_DBD", "CHEMBL550348",  -7.098, "DEFERASIROX: iron chelator, AD-relevant"),
]


# ── GROMACS MDP file templates (same as 17_prep_mmpbsa.py) ─────────────────

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
; NPT equilibration (100 ps)
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


def fetch_smiles(chembl_id: str) -> str:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            data = json.loads(r.read())
        return data["molecule_structures"]["canonical_smiles"]
    except Exception as e:
        log.warning(f"  ChEMBL API failed for {chembl_id}: {e}")
        return ""


def get_formal_charge(smiles: str) -> int:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        return sum(a.GetFormalCharge() for a in mol.GetAtoms())
    except Exception:
        return 0


def extract_model1_pdb(pdbqt_path: Path, out_pdb: Path) -> bool:
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
    tmp = out_pdb.parent / "_tmp_model1.pdbqt"
    tmp.write_text("\n".join(lines) + "\n")
    result = subprocess.run(
        [OBABEL, str(tmp), "-O", str(out_pdb), "-h"],
        capture_output=True, text=True
    )
    tmp.unlink(missing_ok=True)
    return out_pdb.exists() and out_pdb.stat().st_size > 0


def run_acpype(sdf_path: Path, work_dir: Path, charge: int) -> bool:
    """Try AM1-BCC charges first; fall back to Gasteiger if antechamber/sqm fails."""
    for charge_method in ("bcc", "gas"):
        # Clean up failed attempt directory before retry
        acpype_dir = work_dir / "lig.acpype"
        if acpype_dir.exists():
            shutil.rmtree(acpype_dir)

        log.info(f"  acpype -c {charge_method} ...")
        result = subprocess.run(
            [ACPYPE, "-i", str(sdf_path), "-c", charge_method, "-a", "gaff2",
             "-n", str(charge), "-b", "lig", "-o", "gmx"],
            cwd=str(work_dir),
            capture_output=True, text=True, timeout=300
        )
        itp = acpype_dir / "lig_GMX.itp"
        gro = acpype_dir / "lig_GMX.gro"
        if itp.exists() and gro.exists():
            shutil.copy(itp, work_dir / "lig_GMX.itp")
            shutil.copy(gro, work_dir / "lig_GMX.gro")
            for posre_name in ("posre_lig.itp", "lig_GMX_posre.itp"):
                posre = acpype_dir / posre_name
                if posre.exists():
                    shutil.copy(posre, work_dir / "lig_GMX_posre.itp")
                    break
            log.info(f"  GAFF2 topology ready (charge method: {charge_method})")
            return True
        log.warning(f"  acpype -{charge_method} failed")
        if charge_method == "bcc":
            log.info(f"  Retrying with Gasteiger charges (-c gas)...")

    log.error(f"  acpype failed with both bcc and gas charge methods")
    return False


def prep_complex(target: str, chembl_id: str, vina_score: float) -> bool:
    work_dir = MMPBSA_DIR / target / chembl_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Extract Vina pose from tier2 output PDBQT
    vina_pdbqt = VINA_TIER2 / f"{target}_{chembl_id}_out.pdbqt"
    ligand_pdb = work_dir / "ligand_vina.pdb"
    if not vina_pdbqt.exists():
        log.error(f"  Vina PDBQT not found: {vina_pdbqt}")
        return False
    if not extract_model1_pdb(vina_pdbqt, ligand_pdb):
        log.error(f"  Failed to extract model 1 from {vina_pdbqt}")
        return False
    log.info(f"  Pose extracted → {ligand_pdb.name}")

    # 2. Get SMILES → SDF
    ligand_sdf = work_dir / "ligand.sdf"
    if not ligand_sdf.exists():
        smiles = fetch_smiles(chembl_id)
        if not smiles:
            log.error(f"  No SMILES for {chembl_id}")
            return False
        smi_file = work_dir / "ligand.smi"
        smi_file.write_text(f"{smiles}\t{chembl_id}\n")
        result = subprocess.run(
            [OBABEL, "-ismi", str(smi_file),
             "-O", str(ligand_sdf), "--gen3d", "-p", "7.4"],
            capture_output=True, text=True
        )
        if not ligand_sdf.exists():
            log.error(f"  obabel SDF generation failed: {result.stderr[:200]}")
            return False
        log.info(f"  SDF generated from SMILES")
        time.sleep(0.3)
    else:
        smiles = ""
        smi_file = work_dir / "ligand.smi"
        if smi_file.exists():
            smiles = smi_file.read_text().split()[0]
        log.info(f"  SDF already exists, skipping obabel")

    # 3. Run acpype
    if not (work_dir / "lig_GMX.itp").exists():
        smi_file = work_dir / "ligand.smi"
        smiles = smi_file.read_text().split()[0] if smi_file.exists() else ""
        charge = get_formal_charge(smiles)
        log.info(f"  Running acpype (charge={charge})...")
        if not run_acpype(ligand_sdf, work_dir, charge):
            log.error(f"  acpype failed for {chembl_id}")
            return False
    else:
        log.info(f"  GAFF2 topology already exists, skipping acpype")

    # 4. Copy receptor PDB
    rec_pdb_name = RECEPTOR_PDB.get(target)
    if not rec_pdb_name:
        log.error(f"  Unknown target: {target}")
        return False
    rec_src = RECEPTOR_DIR / rec_pdb_name
    if not rec_src.exists():
        log.error(f"  Receptor PDB not found: {rec_src}")
        return False
    shutil.copy(rec_src, work_dir / "receptor.pdb")

    # 5. Write MDP + mmpbsa.in
    (work_dir / "em.mdp").write_text(EM_MDP)
    (work_dir / "nvt.mdp").write_text(NVT_MDP)
    (work_dir / "npt.mdp").write_text(NPT_MDP)
    (work_dir / "md.mdp").write_text(MD_MDP)
    (work_dir / "mmpbsa.in").write_text(MMPBSA_IN)

    # 6. Write metadata
    (work_dir / "meta.txt").write_text(
        f"target={target}\nchembl={chembl_id}\nvina={vina_score}\n"
    )

    return True


def main():
    MMPBSA_DIR.mkdir(parents=True, exist_ok=True)
    log.info(f"Preparing {len(SELECTED)} Tier-2 complexes for MM-GBSA")

    manifest_lines = []
    n_ok = n_fail = 0

    for target, chembl_id, vina, rationale in SELECTED:
        log.info(f"\n── {target} / {chembl_id}  Vina={vina:.3f}  ({rationale}) ──")
        ok = prep_complex(target, chembl_id, vina)
        if ok:
            manifest_lines.append(str(MMPBSA_DIR / target / chembl_id))
            n_ok += 1
            log.info(f"  OK")
        else:
            n_fail += 1
            log.warning(f"  FAILED — will be skipped in SLURM array")

    MANIFEST.write_text("\n".join(manifest_lines) + "\n")
    log.info(f"\n{'='*60}")
    log.info(f"Prepared: {n_ok}   Failed: {n_fail}")
    log.info(f"Manifest: {MANIFEST}  ({n_ok} entries)")
    log.info(f"Next: edit 17_run_mmpbsa.slurm  #SBATCH --array=0-{n_ok-1}")
    log.info(f"      update --output/--error paths to mmpbsa_tier2-%A_%a")
    log.info(f"      sbatch code/slurm/17_run_mmpbsa.slurm")


if __name__ == "__main__":
    main()
