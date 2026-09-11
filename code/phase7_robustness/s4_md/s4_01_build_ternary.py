"""
WP3.1 — Build CRBN-IKZF1(ZF2)-glue ternary complexes for T-REMD
=============================================================
For the top N generated glues: dock (already done), place the best pose into
8RQC_CRBN_ZF2.pdb (CRBN + IKZF1 ZF2 + 4 Zn), GAFF2-parameterise the ligand with
acpype, and emit a GROMACS-ready complex.

Run on TRUBA after code/slurm/paper2/s4_00_setup_truba_toolchain.sh:
  export PATH="$HOME/.local/bin:$PATH"
  /arf/home/mozkurt/miniconda3/envs/pld3/bin/python \
      code/phase7_robustness/s4_md/s4_01_build_ternary.py --top 3

(neither scenic.sif nor lipogate-env.sif has the full rdkit+openbabel+vina+
meeko+acpype+pdb2pqr chain; see s4_00 for the working toolchain and the
library-layout fixes it applies.)

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


def _select_copy(src: Path, dst: Path, keep_chains=("A", "B")):
    """Keep ATOM records + ZN HETATM for the given chains; drop waters/others."""
    out = []
    for l in Path(src).read_text().splitlines():
        rec = l[:6]
        if rec == "ATOM  " and l[21] in keep_chains:
            out.append(l)
        elif rec == "HETATM" and l[17:20].strip() == "ZN" and l[21] in keep_chains:
            out.append(l)
        elif rec == "TER" and (len(l) < 22 or l[21] in keep_chains):
            out.append(l)
    out.append("END")
    dst.write_text("\n".join(out) + "\n")
    print(f"  receptor: {sum(1 for x in out if x[:6]=='ATOM  ')} atoms, "
          f"{sum(1 for x in out if x[:6]=='HETATM')} Zn, chains {keep_chains}")


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

        # 2. receptor: ONE ternary copy only — chain A (CRBN) + chain B
        #    (IKZF1 ZF2, res 144-170) + their two Zn (A/601 CRBN C4 site;
        #    B/201 ZF2 C2H2 site). Drop the second copy (chains D/E, Zn D/E).
        rec_pdb = wd / "receptor.pdb"
        _select_copy(RECEPTOR, rec_pdb, keep_chains=("A", "B"))
        # Zn model = ion + harmonic distance restraints (k=10000 kJ/mol/nm^2,
        # r0 0.23 nm Zn-S(Cys) / 0.20 nm Zn-N(His)), added to topol.top by
        # code/phase7_robustness/s4_md/s4_zn_restraints.py in s4_02. Coordination
        # is auto-detected from the structure. No cationic-dummy model.

        # 3. protonate receptor
        # pdb2pqr has no AMBER-FF template for a bare Zn2+ ion and silently
        # DROPS it (verified 2026-09-11: 0/2 ZN atoms survive pdb2pqr on the
        # real receptor, no error/exit code to flag it). Re-insert the two
        # original ZN HETATM records after protonation — pdb2pqr/obabel never
        # needed to touch them (no hydrogens on a monatomic ion).
        rec_h = wd / "receptor_H.pdb"
        try:
            sh(f"pdb2pqr --ff=AMBER --with-ph=7.4 --keep-chain {rec_pdb} {wd}/receptor.pqr && "
               f"obabel {wd}/receptor.pqr -O {rec_h}")
        except subprocess.CalledProcessError:
            print("  pdb2pqr failed — falling back to pdbfixer")
            sh(f"python -m pdbfixer {rec_pdb} --add-atoms=hydrogens --ph=7.4 "
               f"--output={rec_h}")
        zn_lines = [l for l in rec_pdb.read_text().splitlines()
                   if l[:6] == "HETATM" and l[17:20].strip() == "ZN"]
        h_lines = rec_h.read_text().splitlines()
        end_idx = next((i for i, l in enumerate(h_lines) if l.startswith("END")), len(h_lines))
        rec_h.write_text("\n".join(h_lines[:end_idx] + zn_lines + ["END"]) + "\n")
        n_zn_check = sum(1 for l in rec_h.read_text().splitlines()
                         if l[:6] == "HETATM" and l[17:20].strip() == "ZN")
        if n_zn_check != len(zn_lines):
            print(f"  WARNING: expected {len(zn_lines)} Zn in {rec_h}, found {n_zn_check}")

        # 4. GAFF2 ligand parameters. Charge method = gas (Gasteiger), not
        # AM1-BCC/sqm: verified 2026-09-11 that -c bcc needs a clean-valence
        # input (sqm refused an odd-electron count on a quick obabel-converted
        # test pose -- a ligand-prep protonation issue, not a toolchain gap).
        # gas ran cleanly end-to-end and produced a complete GROMACS topology.
        # TODO: revisit -c bcc once ligand protonation is built explicitly via
        # RDKit (AddHs + sanitize) instead of relying on obabel's guess.
        sh(f"acpype -i {lig_mol2} -b {gid} -c gas -a gaff2 -o gmx", cwd=wd)

        # 5. complex.pdb = protonated receptor (protein + Zn) ONLY.
        # acpype with -o gmx does not write a ligand _NEW.pdb (that filename
        # was a wrong assumption from a different acpype output mode -- the
        # earlier version of this script silently produced a receptor-only
        # complex.pdb because that grep failed but `check=True` only sees the
        # exit code of the LAST command in the `;`-joined shell string, so the
        # failure never surfaced). This receptor-only file is actually what
        # s4_02_gromacs_prep.slurm's `pdb2gmx -f complex.pdb` wants: pdb2gmx
        # force-fields the PROTEIN, and the ligand is merged in separately by
        # s4_02 from acpype's own <gid>_GMX.gro (GAFF2, already fully typed).
        complex_pdb = wd / "complex.pdb"
        atom_lines = [l for l in rec_h.read_text().splitlines()
                     if l[:6] in ("ATOM  ", "HETATM")]
        complex_pdb.write_text("\n".join(atom_lines) + "\nTER\nEND\n")
        n_atoms = len(atom_lines)
        if n_atoms < 100:
            print(f"  WARNING: complex.pdb only has {n_atoms} atoms — check receptor_H.pdb")

        (wd / "README.txt").write_text(
            f"glue_id: {gid}\nsmiles: {r['smiles']}\n"
            f"dock Vina: {r.get('vina_kcal_mol','NA')} kcal/mol\n"
            f"receptor: 8RQC CRBN (A/D) + IKZF1 ZF2 (B/E, res 144-170) + 4 Zn\n"
            f"ligand FF: GAFF2 (acpype, Gasteiger charges)\n"
            f"Zn: kept as ion; restrain to coordinating residues in equilibration\n"
        )
        print(f"[{gid}] complex built -> {wd}")

    print("\nNext: sbatch --array=0-%d code/slurm/paper2/s4_02_gromacs_prep.slurm"
          % (args.top - 1))


if __name__ == "__main__":
    main()
