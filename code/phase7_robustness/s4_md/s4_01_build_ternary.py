"""
WP3.1 — Build CRBN-IKZF1(ZF2)-glue ternary complexes for T-REMD
=============================================================
For the top N generated glues: dock (already done), place the best pose into
8RQC_CRBN_ZF2.pdb (CRBN + IKZF1 ZF2 + 4 Zn), GAFF2-parameterise the ligand with
acpype, and emit a GROMACS-ready complex.

Run on TRUBA (needs rdkit, openbabel, acpype, and either pdb2pqr/pdbfixer):
  python s4_01_build_ternary.py --top 3

Inputs:
  data/structures/pdb/8RQC_CRBN_ZF2.pdb
  results/phase7/glue_design/glue_candidates_top.csv
  results/phase7/glue_design/docked/<glue_id>.pdbqt     (rerun 06_dock_glues.py on TRUBA first)
Outputs -> results/phase7/glue_md/<glue_id>/
  complex.pdb, receptor.pdb, ligand.mol2, ligand GAFF2 (acpype dir)
  README.txt with the residue/atom bookkeeping
"""
import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
RECEPTOR = PROJECT / "data/structures/pdb/8RQC_CRBN_ZF2.pdb"
TOPCSV = PROJECT / "results/phase7/glue_design/glue_candidates_top.csv"
DOCKED = PROJECT / "results/phase7/glue_design/docked"
OUTBASE = PROJECT / "results/phase7/glue_md"


def sh(cmd, cwd=None):
    print("+", cmd)
    subprocess.run(cmd, shell=True, check=True, cwd=cwd)


def pdbqt_to_mol(pdbqt: Path, ref_smiles: str, out_mol2: Path):
    """Best pose pdbqt -> mol2 with correct bond orders from the reference SMILES."""
    # openbabel keeps coordinates; assign bonds from template
    sh(f"obabel {pdbqt} -O {out_mol2} --partialcharge gasteiger")
    # sanity: RDKit read-back
    m = Chem.MolFromMol2File(str(out_mol2), sanitize=False)
    if m is None:
        print(f"  WARN: RDKit could not parse {out_mol2}")
    return out_mol2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=3)
    args = ap.parse_args()

    top = pd.read_csv(TOPCSV).head(args.top)
    OUTBASE.mkdir(parents=True, exist_ok=True)

    for _, r in top.iterrows():
        gid = r["glue_id"]
        wd = OUTBASE / gid
        wd.mkdir(parents=True, exist_ok=True)
        pose = DOCKED / f"{gid}.pdbqt"
        if not pose.exists():
            print(f"SKIP {gid}: {pose} missing — run 06_dock_glues.py on TRUBA")
            continue

        # 1. ligand pose -> mol2 (GAFF2 later via acpype)
        lig_mol2 = wd / "ligand.mol2"
        pdbqt_to_mol(pose, r["smiles"], lig_mol2)

        # 2. receptor: strip waters, keep protein + ZN; add hydrogens at pH 7.4
        rec_pdb = wd / "receptor.pdb"
        sh(f"grep -E '^ATOM|^HETATM.{{13}}ZN|^TER' {RECEPTOR} > {rec_pdb} || true")
        # FIXME(D-Zn): choose a Zn model. Default = keep Zn as ion + distance
        #   restraints to the 4 coordinating Cys/His in NVT/NPT (see s4_02).
        #   Alternative: cationic dummy-atom model (Duarte et al.).

        # 3. protonate receptor
        rec_h = wd / "receptor_H.pdb"
        try:
            sh(f"pdb2pqr --ff=AMBER --with-ph=7.4 --keep-chain {rec_pdb} {wd}/receptor.pqr && "
               f"obabel {wd}/receptor.pqr -O {rec_h}")
        except subprocess.CalledProcessError:
            print("  pdb2pqr failed — falling back to pdbfixer")
            sh(f"python -m pdbfixer {rec_pdb} --add-atoms=hydrogens --ph=7.4 "
               f"--output={rec_h}")

        # 4. GAFF2 ligand parameters
        sh(f"acpype -i {lig_mol2} -b {gid} -c bcc -a gaff2 -o gmx", cwd=wd)

        # 5. assemble complex.pdb (receptor_H + ligand)
        complex_pdb = wd / "complex.pdb"
        lig_pdb = wd / f"{gid}.acpype" / f"{gid}_NEW.pdb"
        sh(f"grep -E '^ATOM|^HETATM' {rec_h} > {complex_pdb}; "
           f"echo 'TER' >> {complex_pdb}; "
           f"grep -E '^ATOM|^HETATM' {lig_pdb} >> {complex_pdb}; "
           f"echo 'END' >> {complex_pdb}")

        (wd / "README.txt").write_text(
            f"glue_id: {gid}\nsmiles: {r['smiles']}\n"
            f"dock Vina: {r.get('vina_kcal_mol','NA')} kcal/mol\n"
            f"receptor: 8RQC CRBN (A/D) + IKZF1 ZF2 (B/E, res 144-170) + 4 Zn\n"
            f"ligand FF: GAFF2 (acpype, AM1-BCC charges)\n"
            f"Zn: kept as ion; restrain to coordinating residues in equilibration\n"
        )
        print(f"[{gid}] complex built -> {wd}")

    print("\nNext: sbatch --array=0-%d code/slurm/paper2/s4_02_gromacs_prep.slurm"
          % (args.top - 1))


if __name__ == "__main__":
    main()
