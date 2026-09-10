"""
Section 4B.6 — Positional-plausibility docking of top glue candidates
====================================================================
Docks the top-ranked generated glues into the 8RQC CRBN-IKZF1(ZF2) ternary
interface, reusing the exact receptor and grid box from the published screen.
This is a docking sanity check on where the molecule can sit — NOT MD, NOT an
affinity claim.

Inputs:
  results/phase7/glue_design/glue_candidates_top.csv
  data/docking/receptors/IKZF1_8RQC_CRBN_prep.pdbqt
  data/docking/configs/IKZF1_8RQC_CRBN.conf   (center 0.566,-2.235,4.760; size 20^3)
Outputs:
  results/phase7/glue_design/glue_docking.csv
  results/phase7/glue_design/docked/<glue_id>.pdbqt   (best pose)

Reference: lenalidomide anchor docked the same way (baseline).
Note: exhaustiveness lowered to 16 (published screen used 32) for a local run.
"""
import os
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
GD = PROJECT / "results/phase7/glue_design"
RECEPTOR = PROJECT / "data/docking/receptors/IKZF1_8RQC_CRBN_prep.pdbqt"
DOCKED = GD / "docked"
DOCKED.mkdir(parents=True, exist_ok=True)

CENTER = (0.566, -2.235, 4.760)
BOX = (20.0, 20.0, 20.0)
EXHAUST = 16
N_POSES = 10
N_TOP = 25
ANCHOR_SMILES = "O=C1CCC(N2Cc3cccc(N)c3C2=O)C(=O)N1"


def prep_ligand_pdbqt(smiles, seed=42):
    from meeko import MoleculePreparation
    try:
        from meeko import PDBQTWriterLegacy
        legacy = True
    except ImportError:
        legacy = False
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=seed) != 0:
        if AllChem.EmbedMolecule(m, randomSeed=seed, useRandomCoords=True) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    prep = MoleculePreparation()
    setups = prep.prepare(m)
    setup = setups[0] if isinstance(setups, (list, tuple)) else setups
    if legacy:
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
        return pdbqt if ok else None
    return setup.write_pdbqt_string()


def dock_one(v, smiles, out_pdbqt):
    lig = prep_ligand_pdbqt(smiles)
    if lig is None:
        return None
    v.set_ligand_from_string(lig)
    v.dock(exhaustiveness=EXHAUST, n_poses=N_POSES)
    e = v.energies(n_poses=1)
    v.write_poses(str(out_pdbqt), n_poses=1, overwrite=True)
    return float(e[0][0])   # total kcal/mol of best pose


def main():
    from vina import Vina
    if not RECEPTOR.exists():
        sys.exit(f"receptor not found: {RECEPTOR}")
    keep = pd.read_csv(GD / "glue_candidates_top.csv").head(N_TOP)

    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(rigid_pdbqt_filename=str(RECEPTOR))
    v.compute_vina_maps(center=list(CENTER), box_size=list(BOX))

    rows = []
    # baseline: lenalidomide anchor
    a = dock_one(v, ANCHOR_SMILES, DOCKED / "anchor_lenalidomide.pdbqt")
    print(f"anchor (lenalidomide): {a}")
    rows.append({"glue_id": "anchor_lenalidomide", "smiles": ANCHOR_SMILES,
                 "vina_kcal_mol": a, "linker": "", "cns_mpo_proxy": ""})

    for _, r in keep.iterrows():
        try:
            score = dock_one(v, r["smiles"], DOCKED / f"{r['glue_id']}.pdbqt")
        except Exception as e:
            print(f"  {r['glue_id']} failed: {e}")
            score = None
        rows.append({"glue_id": r["glue_id"], "smiles": r["smiles"],
                     "vina_kcal_mol": score, "linker": r["linker"],
                     "cns_mpo_proxy": r["cns_mpo_proxy"]})
        print(f"  {r['glue_id']:10s} {score}")

    out = pd.DataFrame(rows)
    out.to_csv(GD / "glue_docking.csv", index=False)
    scored = out[out["vina_kcal_mol"].notna() & (out["glue_id"] != "anchor_lenalidomide")]
    print(f"\n{len(scored)}/{len(keep)} docked")
    if len(scored):
        print(f"Vina range: {scored['vina_kcal_mol'].min():.2f} .. "
              f"{scored['vina_kcal_mol'].max():.2f} kcal/mol  (anchor {a:.2f})")
        print(scored.nsmallest(10, "vina_kcal_mol")[
            ["glue_id", "linker", "vina_kcal_mol"]].to_string(index=False))
    print(f"\nWrote {GD}/glue_docking.csv and docked/*.pdbqt")


if __name__ == "__main__":
    main()
