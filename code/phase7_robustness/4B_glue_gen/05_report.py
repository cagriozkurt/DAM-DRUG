"""
Section 4B.5 — Glue-generation report
=====================================
Input:  results/phase7/glue_design/generated_library.csv
Outputs:
  results/phase7/glue_design/glue_candidates_top.csv   (ranked passers)
  results/phase7/glue_design/glue_top.sdf              (<=50, 3D conformer)
  results/phase7/glue_design/glue_grid.png             (top 24, 2D)
"""
import os
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
GD = PROJECT / "results/phase7/glue_design"
TOP_N_SDF = 50


def main():
    df = pd.read_csv(GD / "generated_library.csv")
    df = df.drop_duplicates("inchikey")
    tier = "passes_all_gates" if df["passes_all_gates"].sum() >= 10 else "passes_fallback"
    keep = df[df[tier]].sort_values(["cns_mpo_proxy", "qed"], ascending=False).reset_index(drop=True)
    keep.insert(0, "rank", keep.index + 1)
    keep.to_csv(GD / "glue_candidates_top.csv", index=False)
    print(f"tier used: {tier}; {len(keep)} unique ranked candidates")

    # 3D SDF for the top N
    w = Chem.SDWriter(str(GD / "glue_top.sdf"))
    for _, r in keep.head(TOP_N_SDF).iterrows():
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            continue
        mh = Chem.AddHs(m)
        if AllChem.EmbedMolecule(mh, randomSeed=42) != 0:
            AllChem.EmbedMolecule(mh, randomSeed=42, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(mh)
        except Exception:
            pass
        mh.SetProp("_Name", r["glue_id"])
        for c in ["linker", "mw", "clogp", "tpsa", "cns_mpo_proxy", "qed"]:
            mh.SetProp(c, str(r[c]))
        w.write(mh)
    w.close()

    # 2D grid
    mols, legs = [], []
    for _, r in keep.head(24).iterrows():
        m = Chem.MolFromSmiles(r["smiles"])
        if m:
            mols.append(m)
            legs.append(f"{r['glue_id']} MPO{r['cns_mpo_proxy']} "
                        f"cLogP{r['clogp']:.1f} QED{r['qed']:.2f}")
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 230), legends=legs)
    img.save(str(GD / "glue_grid.png"))
    print(f"Wrote glue_candidates_top.csv, glue_top.sdf, glue_grid.png")
    print(keep.head(10)[["rank", "glue_id", "linker", "mw", "clogp", "tpsa",
                         "cns_mpo_proxy", "qed"]].to_string(index=False))


if __name__ == "__main__":
    main()
