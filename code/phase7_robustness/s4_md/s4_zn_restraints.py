"""
WP3 — Zn2+ harmonic distance restraints for the ternary-complex topology
======================================================================
IKZF1 ZF2 and CRBN each hold a structural Zn2+. Plain Lennard-Jones parameters
let the ion drift or leave at 320 K T-REMD, collapsing the degron hairpin whose
stability is the readout. This appends harmonic distance restraints (GROMACS
bonded function type 6 — harmonic, no exclusions) between each Zn2+ and its
coordinating atoms, into an [ intermolecular_interactions ] block at the end of
topol.top.

Coordination is auto-detected from the structure (any Cys-SG / His-NE2 / His-ND1
/ Asp-Glu carboxylate O within CUTOFF of a ZN), so it is robust to residue
renumbering between the 8RQC construct and pdb2gmx output.

Restraint parameters (per user spec / metalloprotein MD practice):
  Zn--S(Cys)      r0 = 0.230 nm
  Zn--N(His)      r0 = 0.200 nm
  Zn--O(Asp/Glu)  r0 = 0.200 nm
  k = 10000 kJ mol^-1 nm^-2   (all)

Usage:
  python s4_zn_restraints.py <complex_full.gro> <topol.top>
  # run AFTER pdb2gmx + ligand splice + concatenation into complex_full.gro,
  # BEFORE editconf/solvate. Global 1-based indices from the solute-only .gro
  # stay valid once solvent/ions are appended to [ molecules ].
"""
import sys
from pathlib import Path

import numpy as np

CUTOFF = 0.29                 # nm
K = 10000.0                   # kJ/mol/nm^2
R0 = {"S": 0.230, "N": 0.200, "O": 0.200}
COORD_ATOMS = {"SG": "S", "NE2": "N", "ND1": "N",
               "OD1": "O", "OD2": "O", "OE1": "O", "OE2": "O"}


def read_gro(path):
    lines = Path(path).read_text().splitlines()
    natoms = int(lines[1])
    atoms = []
    for i in range(2, 2 + natoms):
        ln = lines[i]
        resnum = int(ln[0:5]); resname = ln[5:10].strip()
        aname = ln[10:15].strip(); anum = int(ln[15:20])
        x, y, z = float(ln[20:28]), float(ln[28:36]), float(ln[36:44])
        atoms.append({"idx": i - 1, "resnum": resnum, "resname": resname,
                      "aname": aname, "pos": np.array([x, y, z])})
    return atoms


def main():
    gro, top = sys.argv[1], sys.argv[2]
    atoms = read_gro(gro)
    zns = [a for a in atoms if a["resname"].upper() == "ZN" or a["aname"].upper() == "ZN"]
    if not zns:
        print("no ZN found in", gro, "— nothing to restrain")
        return

    coord = [a for a in atoms if a["aname"] in COORD_ATOMS]
    cpos = np.array([a["pos"] for a in coord])

    bonds = []   # (ai, aj, r0, comment)
    for zn in zns:
        d = np.linalg.norm(cpos - zn["pos"], axis=1)
        for j in np.where(d < CUTOFF)[0]:
            a = coord[j]
            elem = COORD_ATOMS[a["aname"]]
            bonds.append((zn["idx"], a["idx"], R0[elem],
                          f"; ZN{zn['resnum']} -- {a['resname']}{a['resnum']}.{a['aname']} "
                          f"(obs {d[j]*10:.2f} A)"))
    if not bonds:
        print("WARNING: ZN present but no coordinating atom within "
              f"{CUTOFF*10:.1f} A — check the input structure")
        return

    block = ["", "; ── Zn2+ harmonic distance restraints (s4_zn_restraints.py) ──",
             "[ intermolecular_interactions ]", "[ bonds ]",
             ";  ai    aj  funct       b0(nm)      kb(kJ/mol/nm^2)"]
    for ai, aj, r0, cm in bonds:
        block.append(f"{ai:6d} {aj:6d}      6   {r0:10.4f}   {K:12.1f}   {cm}")
    block.append("")

    tp = Path(top)
    txt = tp.read_text().rstrip()
    if "[ intermolecular_interactions ]" in txt:
        print("topol.top already has [ intermolecular_interactions ] — not touching it")
        return
    tp.write_text(txt + "\n" + "\n".join(block) + "\n")
    print(f"added {len(bonds)} Zn restraints to {top}:")
    for ai, aj, r0, cm in bonds:
        print(f"  {ai:5d}-{aj:5d}  r0={r0:.3f} nm  {cm}")


if __name__ == "__main__":
    main()
