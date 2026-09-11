"""
Section 4 tail — off-target counter-docking of the top-25 CRBN glue candidates
================================================================================
Docks the same top-25 glue_candidates_top.csv set already docked into the
8RQC on-target interface (results/phase7/glue_design/docked/) against three
CNS off-target receptors prepped in Phase 1 (code/phase4_docking/
08_prep_selectivity_receptors.py): DRD2 (6CM4), HTR2A (6A94), HERG (7CN1).

Reuses the exact receptor .pdbqt + grid box (center/size) from Phase 1's
data/docking/configs/*.conf so off-target scores are directly comparable to
the accepted paper's own selectivity_table.csv methodology.

Output: results/phase7/selectivity/glue_offtarget_docking.csv
  (glue_id, smiles, vina_drd2, vina_htr2a, vina_herg)
"""
import os
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
GD = PROJECT / "results/phase7/glue_design"
SEL = PROJECT / "results/phase7/selectivity"
RECV = PROJECT / "data/docking/receptors"
N_TOP = 25

OFFTARGETS = {
    "DRD2":  {"receptor": RECV / "DRD2_prep.pdbqt",
              "center": (9.925, 5.846, -9.582), "box": (22, 22, 22)},
    "HTR2A": {"receptor": RECV / "HTR2A_prep.pdbqt",
              "center": (16.106, -0.864, 57.074), "box": (22, 22, 22)},
    "HERG":  {"receptor": RECV / "HERG_prep.pdbqt",
              "center": (141.441, 141.443, 173.423), "box": (26, 26, 26)},
}
EXHAUST = 16  # matches 06_dock_glues.py's local-run reduction (published screen used 32)
N_POSES = 10


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


def dock_one(v, smiles):
    lig = prep_ligand_pdbqt(smiles)
    if lig is None:
        return None
    v.set_ligand_from_string(lig)
    v.dock(exhaustiveness=EXHAUST, n_poses=N_POSES)
    e = v.energies(n_poses=1)
    return float(e[0][0])


def main():
    from vina import Vina
    keep = pd.read_csv(GD / "glue_candidates_top.csv").head(N_TOP)
    print(f"counter-docking {len(keep)} top candidates against "
          f"{list(OFFTARGETS.keys())}")

    results = {name: {} for name in OFFTARGETS}
    for name, cfg in OFFTARGETS.items():
        recv = cfg["receptor"]
        if not recv.exists():
            print(f"WARNING: receptor missing, skipping {name}: {recv}")
            continue
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(rigid_pdbqt_filename=str(recv))
        v.compute_vina_maps(center=list(cfg["center"]), box_size=list(cfg["box"]))
        for _, r in keep.iterrows():
            try:
                score = dock_one(v, r["smiles"])
            except Exception as e:
                print(f"  {name} {r['glue_id']}: FAILED ({e})")
                score = None
            results[name][r["glue_id"]] = score
            print(f"  {name} {r['glue_id']}: {score}")

    rows = []
    for _, r in keep.iterrows():
        row = {"glue_id": r["glue_id"], "smiles": r["smiles"],
               "vina_ontarget_8rqc": r["vina_kcal_mol"]}
        for name in OFFTARGETS:
            row[f"vina_{name.lower()}"] = results[name].get(r["glue_id"])
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(SEL / "glue_offtarget_docking.csv", index=False)
    print(f"-> {SEL / 'glue_offtarget_docking.csv'}")


if __name__ == "__main__":
    main()
