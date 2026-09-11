"""
Section 4 tail — hERG liability QSAR model
===========================================
Pulls hERG (KCNH2, CHEMBL240) binding-assay bioactivity from the ChEMBL REST
API, trains a RandomForest classifier on 2048-bit Morgan fingerprints
(radius 2), and fits a Tree SHAP explainer for interpretability.

Active/inactive cutoff: pChEMBL >= 5 (IC50/Ki <= 10 uM), the standard
hERG-liability screening threshold used across the QSAR literature (e.g.
Wang et al. 2020 in silico hERG panels; FDA ICH S7B nonclinical guidance
uses 10x margin over free Cmax with the same 10 uM anchor).

Output:
  results/phase7/selectivity/herg_qsar_model.pkl   (sklearn RandomForest)
  results/phase7/selectivity/herg_training_data.csv
  results/phase7/selectivity/herg_qsar_metrics.json (ROC-AUC, CV, n actives/inactives)
  results/phase7/selectivity/herg_shap_summary.png

Run: conda run -n lipogate python code/phase7_robustness/4C_selectivity/01_train_herg_qsar.py
"""
import json
import os
import time
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
OUT = PROJECT / "results/phase7/selectivity"
OUT.mkdir(parents=True, exist_ok=True)

CHEMBL_TARGET = "CHEMBL240"  # KCNH2 / hERG
PCHEMBL_ACTIVE_CUTOFF = 5.0  # IC50/Ki <= 10 uM -> active (hERG liability)
FP_BITS = 2048
FP_RADIUS = 2
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


def fetch_chembl_herg(max_pages=200, page_size=1000):
    rows, offset = [], 0
    for page in range(max_pages):
        url = (f"{BASE_URL}?target_chembl_id={CHEMBL_TARGET}"
               f"&standard_type__in=IC50,Ki&assay_type=B"
               f"&pchembl_value__isnull=false&limit={page_size}&offset={offset}")
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as resp:
                    data = json.load(resp)
                break
            except Exception as e:
                if attempt == 2:
                    raise
                time.sleep(2 * (attempt + 1))
        acts = data.get("activities", [])
        if not acts:
            break
        for a in acts:
            rows.append({
                "molecule_chembl_id": a.get("molecule_chembl_id"),
                "smiles": a.get("canonical_smiles"),
                "pchembl_value": a.get("pchembl_value"),
                "standard_type": a.get("standard_type"),
                "standard_value": a.get("standard_value"),
                "standard_units": a.get("standard_units"),
            })
        print(f"  page {page}: offset={offset} -> {len(acts)} rows (total {len(rows)})")
        nxt = data.get("page_meta", {}).get("next")
        if not nxt:
            break
        offset += page_size
    return pd.DataFrame(rows)


def morgan_fp(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(m, FP_RADIUS, nBits=FP_BITS)
    arr = np.zeros((FP_BITS,), dtype=np.int8)
    for b in fp.GetOnBits():
        arr[b] = 1
    return arr


def main():
    raw_path = OUT / "herg_training_data_raw.csv"
    if raw_path.exists():
        print(f"reusing cached ChEMBL pull: {raw_path}")
        df = pd.read_csv(raw_path)
    else:
        print(f"fetching hERG (target={CHEMBL_TARGET}) bioactivity from ChEMBL...")
        df = fetch_chembl_herg()
        df.to_csv(raw_path, index=False)
        print(f"fetched {len(df)} raw activity rows -> {raw_path}")

    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    df = df.dropna(subset=["smiles", "pchembl_value"])
    # one row per molecule: take the most potent (max pchembl) measurement
    df = df.sort_values("pchembl_value", ascending=False).drop_duplicates("molecule_chembl_id")
    df["active"] = (df["pchembl_value"] >= PCHEMBL_ACTIVE_CUTOFF).astype(int)
    print(f"deduplicated to {len(df)} molecules; "
          f"{df['active'].sum()} active (hERG liability) / {(1 - df['active']).sum()} inactive")

    fps, keep_idx = [], []
    for i, smi in enumerate(df["smiles"]):
        fp = morgan_fp(smi)
        if fp is not None:
            fps.append(fp)
            keep_idx.append(i)
    df = df.iloc[keep_idx].reset_index(drop=True)
    X = np.vstack(fps)
    y = df["active"].values
    print(f"{len(df)} molecules parsed to valid fingerprints")

    df.to_csv(OUT / "herg_training_data.csv", index=False)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=42)
    clf = RandomForestClassifier(
        n_estimators=500, max_depth=None, min_samples_leaf=2,
        n_jobs=-1, class_weight="balanced", random_state=42)
    clf.fit(X_train, y_train)
    test_auc = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_auc = cross_val_score(clf, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)

    # refit on all data for the deployed model
    clf_full = RandomForestClassifier(
        n_estimators=500, max_depth=None, min_samples_leaf=2,
        n_jobs=-1, class_weight="balanced", random_state=42)
    clf_full.fit(X, y)

    import pickle
    with open(OUT / "herg_qsar_model.pkl", "wb") as f:
        pickle.dump({"model": clf_full, "fp_bits": FP_BITS, "fp_radius": FP_RADIUS,
                     "active_cutoff_pchembl": PCHEMBL_ACTIVE_CUTOFF}, f)

    metrics = {
        "n_molecules": int(len(df)),
        "n_active": int(df["active"].sum()),
        "n_inactive": int((1 - df["active"]).sum()),
        "held_out_test_auc": float(test_auc),
        "cv5_auc_mean": float(cv_auc.mean()),
        "cv5_auc_std": float(cv_auc.std()),
        "cv5_auc_folds": [float(v) for v in cv_auc],
    }
    (OUT / "herg_qsar_metrics.json").write_text(json.dumps(metrics, indent=2))
    print(json.dumps(metrics, indent=2))

    # SHAP (Tree explainer) on the held-out test set
    try:
        import shap
        explainer = shap.TreeExplainer(clf)
        sv = explainer.shap_values(X_test)
        sv1 = sv[1] if isinstance(sv, list) else sv[:, :, 1] if sv.ndim == 3 else sv
        fig = plt.figure(figsize=(7, 5))
        shap.summary_plot(sv1, X_test, feature_names=[f"bit_{i}" for i in range(FP_BITS)],
                           show=False, max_display=15)
        plt.tight_layout()
        plt.savefig(OUT / "herg_shap_summary.png", dpi=150)
        plt.close(fig)
        print(f"SHAP summary -> {OUT / 'herg_shap_summary.png'}")
    except ImportError:
        print("WARNING: shap not installed (pip install shap) -- skipping interpretability plot")


if __name__ == "__main__":
    main()
