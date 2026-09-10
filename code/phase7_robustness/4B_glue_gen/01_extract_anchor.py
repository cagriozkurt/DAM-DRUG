"""
Section 4B.1 — Define the CRBN glutarimide anchor and growth vector
==================================================================
- Extract the CELMoD ligand QFC from PDB 8RQC (reference only).
- Define the fixed CRBN warhead: lenalidomide isoindolinone-glutarimide.
- Mark the growth exit vector at the 4-amino position (CELMoD linker attach point).

Output: results/phase7/glue_design/anchor.json
"""
import os
import json
from pathlib import Path
import urllib.request

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
PDB_8RQC = PROJECT / "data/structures/pdb/8RQC.pdb"
OUT = PROJECT / "results/phase7/glue_design"
OUT.mkdir(parents=True, exist_ok=True)

# lenalidomide: 3-(4-amino-1-oxo-1,3-dihydro-2H-isoindol-2-yl)piperidine-2,6-dione
ANCHOR_SMILES = "O=C1CCC(N2Cc3cccc(N)c3C2=O)C(=O)N1"
# same, with the 4-amino nitrogen flagged as the grow handle (primary aromatic amine)
GROW_HANDLE_SMARTS = "[c:1][NH2:2]"          # react here
# glutarimide + isoindolinone substructure that every product MUST retain
ANCHOR_CORE_SMARTS = "O=C1CCC(N2Cc3ccccc3C2=O)C(=O)N1"
# CRBN tri-tryptophan cage (8RQC numbering, chain A) — for the docking step
CRBN_CAGE_RES = ["HIS378", "TRP380", "TRP386", "TRP400"]


def qfc_from_ccd():
    for url in (
        "https://files.rcsb.org/ligands/download/QFC_ideal.sdf",
        "https://files.rcsb.org/ligands/download/QFC.sdf",
    ):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                data = r.read().decode()
            m = Chem.MolFromMolBlock(data)
            if m:
                return Chem.MolToSmiles(m), url
        except Exception as e:
            print(f"  CCD fetch failed ({url}): {e}")
    return None, None


def qfc_from_pdb():
    if not PDB_8RQC.exists():
        return None
    lines = [l for l in PDB_8RQC.read_text().splitlines()
             if l.startswith("HETATM") and l[17:20].strip() == "QFC"]
    # first copy only (one chain)
    if not lines:
        return None
    chain = lines[0][21]
    block = [l for l in lines if l[21] == chain]
    pdb_block = "\n".join(block) + "\nEND\n"
    m = Chem.MolFromPDBBlock(pdb_block, sanitize=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
        return Chem.MolToSmiles(m)
    except Exception:
        return Chem.MolToSmiles(m, canonical=False)


def main():
    anchor = Chem.MolFromSmiles(ANCHOR_SMILES)
    assert anchor is not None
    core = Chem.MolFromSmarts(ANCHOR_CORE_SMARTS)
    assert anchor.HasSubstructMatch(core), "anchor core SMARTS does not match anchor"

    qfc_smiles, qfc_src = qfc_from_ccd()
    if qfc_smiles is None:
        qfc_smiles = qfc_from_pdb()
        qfc_src = "8RQC.pdb HETATM (bond perception, approximate)"

    rec = {
        "warhead": "lenalidomide isoindolinone-glutarimide",
        "anchor_smiles": Chem.MolToSmiles(anchor),
        "anchor_mw": round(Descriptors.MolWt(anchor), 2),
        "anchor_core_smarts": ANCHOR_CORE_SMARTS,
        "grow_handle_smarts": GROW_HANDLE_SMARTS,
        "grow_handle_desc": "4-amino aromatic nitrogen of the isoindolinone (CELMoD linker attachment point)",
        "qfc_reference_smiles": qfc_smiles,
        "qfc_reference_source": qfc_src,
        "crbn_cage_residues_8rqc": CRBN_CAGE_RES,
        "target_surface": "IKZF1 ZF2 beta-hairpin degron (8RQC chains B/E, res 144-170; G-loop Gln146)",
        "docking_box_from": "data/docking/configs/IKZF1_8RQC_CRBN.conf (center 0.566,-2.235,4.760; size 20^3)",
    }
    (OUT / "anchor.json").write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec, indent=2))

    img = Draw.MolToImage(anchor, size=(400, 300))
    img.save(str(OUT / "anchor.png"))
    print(f"\nWrote {OUT}/anchor.json and anchor.png")


if __name__ == "__main__":
    main()
