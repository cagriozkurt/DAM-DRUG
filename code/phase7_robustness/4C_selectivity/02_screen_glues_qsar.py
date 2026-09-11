"""
Section 4 tail — QSAR triage of CNS-fallback glue candidates for hERG liability
================================================================================
Cheap first-pass filter before the expensive counter-docking step: applies the
trained RandomForest hERG model to all 184 candidates that pass the CNS
fallback gates (results/phase7/glue_design/generated_library.csv,
passes_fallback == True).

Output: results/phase7/selectivity/glue_herg_qsar.csv
  (glue_id, smiles, herg_liability_prob, herg_predicted_active)
"""
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
GD = PROJECT / "results/phase7/glue_design"
SEL = PROJECT / "results/phase7/selectivity"


def morgan_fp(smiles, bits, radius):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=bits)
    arr = np.zeros((bits,), dtype=np.int8)
    for b in fp.GetOnBits():
        arr[b] = 1
    return arr


def main():
    with open(SEL / "herg_qsar_model.pkl", "rb") as f:
        bundle = pickle.load(f)
    model, bits, radius = bundle["model"], bundle["fp_bits"], bundle["fp_radius"]

    lib = pd.read_csv(GD / "generated_library.csv")
    cand = lib[lib["passes_fallback"] == True].copy()  # noqa: E712
    print(f"screening {len(cand)} CNS-fallback candidates for hERG liability")

    fps, keep = [], []
    for i, smi in enumerate(cand["smiles"]):
        fp = morgan_fp(smi, bits, radius)
        if fp is not None:
            fps.append(fp)
            keep.append(i)
    cand = cand.iloc[keep].reset_index(drop=True)
    X = np.vstack(fps)
    proba = model.predict_proba(X)[:, 1]
    cand["herg_liability_prob"] = proba
    cand["herg_predicted_active"] = (proba >= 0.5).astype(int)

    out = cand[["glue_id", "smiles", "mw", "clogp", "tpsa", "cns_mpo_proxy",
                "herg_liability_prob", "herg_predicted_active"]].sort_values(
        "herg_liability_prob")
    out.to_csv(SEL / "glue_herg_qsar.csv", index=False)
    n_flag = int(out["herg_predicted_active"].sum())
    print(f"{len(out)} scored -> {SEL / 'glue_herg_qsar.csv'}")
    print(f"{n_flag}/{len(out)} flagged as predicted hERG-active "
          f"(liability_prob >= 0.5); {len(out) - n_flag} clean")
    print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
