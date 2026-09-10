"""
Section 4B.3 — Grow molecular-glue candidates off the CRBN anchor
===============================================================
Anchor = lenalidomide isoindolinone-glutarimide, grow handle = 4-amino N.
For each pool fragment, attach it to the anchor N through a small linker set
(direct bond, amide, reverse-amide, sulfonamide, acetamide spacer, urea).
Every product retains the intact glutarimide + isoindolinone by construction.

Inputs:
  results/phase7/glue_design/anchor.json
  results/phase7/glue_design/fragment_pool.csv
Outputs:
  results/phase7/glue_design/generated_raw.csv   (chembl-free ids, smiles, linker, source_frag)
  results/phase7/glue_design/generated_raw.sdf   (1 ETKDG conformer each)
"""
import os
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")
random.seed(42)
np.random.seed(42)

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
GD = PROJECT / "results/phase7/glue_design"

MAX_PRODUCT_MW = 550          # hard ceiling pre-filter (gate is <450 later)
MAX_PRODUCTS = 20000

# anchor with the 4-amino N carrying a dummy (map 1) — a secondary-amine handle
ANCHOR_FRAG = "O=C1CCC(N2Cc3cccc(N[*:1])c3C2=O)C(=O)N1"

# linkers: two dummies, map 1 -> bonds to anchor N, map 2 -> bonds to fragment
LINKERS = {
    "direct":       "[*:1][*:2]",
    "amide":        "[*:1]C(=O)[*:2]",          # N-C(=O)-frag
    "reverse_amide":"[*:1]C(=O)[*:2]".replace("C(=O)", "NC(=O)"),  # N-NC(=O) invalid; see below
    "sulfonamide":  "[*:1]S(=O)(=O)[*:2]",
    "acetamide":    "[*:1]C(=O)C[*:2]",         # N-C(=O)-CH2-frag
    "urea":         "[*:1]C(=O)N[*:2]",         # N-C(=O)-NH-frag
    "carbamate":    "[*:1]C(=O)O[*:2]",
}
LINKERS.pop("reverse_amide")   # drop the malformed entry


def set_frag_map(smiles, mapnum):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    dummies = [a for a in m.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) != 1:
        return None
    dummies[0].SetAtomMapNum(mapnum)
    return m


def assemble(anchor_m, linker_m, frag_m):
    """molzip anchor(map1) + linker(map1,map2) + fragment(map2)."""
    try:
        combo = Chem.CombineMols(Chem.CombineMols(anchor_m, linker_m), frag_m)
        prod = Chem.molzip(combo)
        Chem.SanitizeMol(prod)
        for a in prod.GetAtoms():
            a.SetAtomMapNum(0)
        return Chem.MolFromSmiles(Chem.MolToSmiles(prod))
    except Exception:
        return None


def main():
    anchor_json = json.loads((GD / "anchor.json").read_text())
    core = Chem.MolFromSmarts(anchor_json["anchor_core_smarts"])
    pool = pd.read_csv(GD / "fragment_pool.csv")
    print(f"anchor + {len(pool)} fragments x {len(LINKERS)} linkers")

    anchor_m = set_frag_map(ANCHOR_FRAG, 1)
    assert anchor_m is not None

    linker_mols = {}
    for name, s in LINKERS.items():
        lm = Chem.MolFromSmiles(s)
        if lm is None:
            print(f"  bad linker {name}: {s}")
            continue
        linker_mols[name] = lm

    seen = set()
    rows = []
    for _, fr in pool.iterrows():
        fm = set_frag_map(fr["smiles"], 2)
        if fm is None:
            continue
        for lname, lm in linker_mols.items():
            prod = assemble(anchor_m, lm, fm)
            if prod is None:
                continue
            if not prod.HasSubstructMatch(core):
                continue
            mw = Descriptors.MolWt(prod)
            if mw > MAX_PRODUCT_MW:
                continue
            ik = Chem.MolToInchiKey(prod)
            if ik in seen:
                continue
            seen.add(ik)
            rows.append({
                "glue_id": f"GLUE{len(rows):04d}",
                "smiles": Chem.MolToSmiles(prod),
                "inchikey": ik,
                "linker": lname,
                "source_fragment": fr["smiles"],
                "mw": round(mw, 2),
            })
            if len(rows) >= MAX_PRODUCTS:
                break
        if len(rows) >= MAX_PRODUCTS:
            break

    df = pd.DataFrame(rows)
    df.to_csv(GD / "generated_raw.csv", index=False)
    print(f"{len(df)} unique valid products (anchor retained)")
    print(df["linker"].value_counts().to_string())

    # 1 conformer each
    w = Chem.SDWriter(str(GD / "generated_raw.sdf"))
    n_embed = 0
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            continue
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42, useRandomCoords=False) != 0:
            if AllChem.EmbedMolecule(m, randomSeed=42, useRandomCoords=True) != 0:
                continue
        try:
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            pass
        m.SetProp("_Name", r["glue_id"])
        w.write(m)
        n_embed += 1
    w.close()
    print(f"embedded {n_embed} conformers -> generated_raw.sdf")


if __name__ == "__main__":
    main()
