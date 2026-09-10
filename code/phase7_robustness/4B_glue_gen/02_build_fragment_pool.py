"""
Section 4B.2 — BRICS fragment pool from CNS-approved chemical space
=================================================================
BRICS-decompose the Tier-1 CNS-approved and Tier-2 approved compound sets into
mono-attachment fragments to grow off the CRBN anchor.

Inputs:
  data/compounds/tier1_cns_approved.csv   (col: smiles)
  data/compounds/tier2_approved.csv       (col: smiles)

Output:
  results/phase7/glue_design/fragment_pool.csv
    smiles (with one [*] dummy), mw, n_source_occurrences, ring_count
"""
import os
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors, RDConfig
from rdkit.Chem import FilterCatalog

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
COMP = PROJECT / "data/compounds"
OUT = PROJECT / "results/phase7/glue_design/fragment_pool.csv"

MW_MIN, MW_MAX = 40.0, 250.0
MAX_RINGS = 2
# reactive / unstable groups to exclude from fragments
BAD_SMARTS = [
    "[CX3](=O)[Cl,Br,I]",         # acyl halide
    "[CX3H1](=O)[#6]",            # aldehyde (keep chemistry predictable; drop)
    "C1OC1", "C1NC1", "C1SC1",    # epoxide/aziridine/thiirane
    "[N+]#[C-]", "[N-]=[N+]=[N-]", # isocyanide, azide
    "[S;X2][S;X2]",              # disulfide
    "[C,c][N+](=O)[O-]",         # nitro (metabolic liability)
    "[CX3]=[CX3][CX3]=[O]",      # Michael acceptor enone
]


def clean_dummies(frag):
    """Return SMILES with generic [*] dummies (strip BRICS isotope labels)."""
    m = Chem.MolFromSmiles(frag)
    if m is None:
        return None
    for a in m.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetIsotope(0)
            a.SetAtomMapNum(0)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return Chem.MolToSmiles(m)


def main():
    bad = [Chem.MolFromSmarts(s) for s in BAD_SMARTS]
    pains = FilterCatalog.FilterCatalog(
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)

    smi = []
    for fn in ("tier1_cns_approved.csv", "tier2_approved.csv"):
        p = COMP / fn
        if not p.exists():
            print(f"  missing {p}")
            continue
        df = pd.read_csv(p)
        col = "smiles" if "smiles" in df.columns else df.columns[df.columns.str.contains("smi", case=False)][0]
        smi += df[col].dropna().tolist()
    print(f"{len(smi)} source SMILES")

    counter = Counter()
    for s in smi:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        for frag in BRICS.BRICSDecompose(m, keepNonLeafNodes=False):
            cs = clean_dummies(frag)
            if cs:
                counter[cs] += 1

    rows = []
    for cs, n in counter.items():
        m = Chem.MolFromSmiles(cs)
        if m is None:
            continue
        n_dummy = sum(a.GetAtomicNum() == 0 for a in m.GetAtoms())
        if n_dummy != 1:
            continue
        mw = Descriptors.MolWt(m)                 # includes dummy ~ negligible
        if not (MW_MIN <= mw <= MW_MAX):
            continue
        if Descriptors.RingCount(m) > MAX_RINGS:
            continue
        if any(m.HasSubstructMatch(b) for b in bad if b is not None):
            continue
        if pains.HasMatch(m):
            continue
        rows.append({"smiles": cs, "mw": round(mw, 2),
                     "n_source_occurrences": n,
                     "ring_count": Descriptors.RingCount(m),
                     "n_heavy": m.GetNumHeavyAtoms()})

    pool = pd.DataFrame(rows).sort_values("n_source_occurrences", ascending=False)
    pool.to_csv(OUT, index=False)
    print(f"\n{len(pool)} mono-attachment fragments after filtering")
    print(pool.head(15).to_string(index=False))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
