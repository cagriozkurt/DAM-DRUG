"""
WP2.2 — Per-region CellChat input prep from SEA-AD Multiregion 2026 subclass objects
==================================================================================
For a given brain region, walk every downloaded subclass h5ad
(data/raw/SEA-AD/multiregion_2026/subclass_objects/SEAAD_<SC>_multiregion.h5ad),
subset to that region via obs["Brain Region"], pull RAW counts from
.layers["UMIs"], coarsen Subclass -> broad cell type, cap 5,000 cells/broad-type
(R dgCMatrix int32 safety), concatenate, and write the CellChat inputs.

Same output layout / seed / caps / coarsen() as the accepted MTG prep
(code/phase2_LR/01_prep_mtg_for_cellchat.py) so each region's CellChat run is
parameter-identical.

Usage (inside scenic.sif):
  python s3_prep.py <REGION> <subclass_objects_dir>

Outputs -> results/phase7/cellchat_multiregion/<REGION>/prep/
  counts_raw.h5 (CSC), cell_meta.csv, gene_names.csv
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import h5py

PROJDIR = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
MAX_CELLS_PER_TYPE = 5000
SEED = 42
REGION_COL_CANDIDATES = ["Brain Region", "region", "Region", "brain_region",
                         "region_of_interest_acronym"]
SUBCLASS_COL_CANDIDATES = ["Subclass", "subclass", "Supertype", "supertype",
                           "Class", "class"]


def coarsen(ct: str) -> str:
    ct = str(ct)
    if "Micro" in ct or "PVM" in ct or ct == "Immune":
        return "Microglia"
    if "Astro" in ct:
        return "Astrocyte"
    if "Oligo" in ct and "OPC" not in ct:
        return "Oligodendrocyte"
    if "OPC" in ct:
        return "OPC"
    gaba = ["Sst", "Vip", "Pvalb", "Lamp5", "Sncg", "Pax6", "Chandelier",
            "Sst-Chodl", "Lamp5-Lhx6"]
    if any(ct == g or ct.startswith(g + "_") or ct.startswith(g) for g in gaba):
        return "InhNeuron"
    exc = ["L23-IT", "L4-IT", "L5-IT", "L5-ET", "L56-NP", "L6-CT", "L6-IT",
           "L6-IT-Car3", "L6b", "EC-IT", "CA2-4", "CA1", "DG", "Sub-CA1", "IT"]
    if any(e in ct for e in exc):
        return "ExcNeuron"
    if any(x in ct for x in ["Endo", "VLMC", "Peri"]):
        return "Vascular"
    if "Ependymal" in ct:
        return "Ependymal"
    return ct


def pick(cols, cands):
    return next((c for c in cands if c in cols), None)


def main():
    region = sys.argv[1]
    sc_dir = Path(sys.argv[2])
    out = PROJDIR / f"results/phase7/cellchat_multiregion/{region}/prep"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    objs = sorted(sc_dir.glob("SEAAD_*_multiregion.h5ad"))
    if not objs:
        sys.exit(f"no subclass objects in {sc_dir} — run s3_01 first")
    print(f"[{region}] {len(objs)} subclass objects")

    parts = []      # (csr counts cells x genes, meta DataFrame)
    var_ref = None
    for p in objs:
        ad = sc.read_h5ad(p, backed="r")
        rcol = pick(ad.obs.columns, REGION_COL_CANDIDATES)
        scol = pick(ad.obs.columns, SUBCLASS_COL_CANDIDATES)
        if rcol is None:
            print(f"  {p.name}: no region column, skipping"); ad.file.close(); continue
        mask = ad.obs[rcol].astype(str).str.upper() == region.upper()
        n = int(mask.sum())
        if n == 0:
            ad.file.close(); continue
        sub = ad[mask.values]
        ct = (sub.obs[scol].astype(str) if scol else
              pd.Series(p.stem.split("_")[1], index=sub.obs_names))
        broad = ct.apply(coarsen)
        keep = []
        for bt in broad.unique():
            idx = np.where(broad.values == bt)[0]
            if len(idx) > MAX_CELLS_PER_TYPE:
                idx = rng.choice(idx, MAX_CELLS_PER_TYPE, replace=False)
            keep.append(idx)
        keep = np.sort(np.concatenate(keep))
        sub2 = sub[keep]

        layer = "UMIs" if "UMIs" in sub2.layers else None
        X = sub2.layers[layer] if layer else sub2.X   # prefer raw UMIs
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        if var_ref is None:
            var_ref = list(sub2.var_names)
        elif list(sub2.var_names) != var_ref:
            # align to reference gene order
            gi = [list(sub2.var_names).index(g) if g in set(sub2.var_names) else -1
                  for g in var_ref]
            cols = [i for i in gi if i >= 0]
            X = X[:, cols]

        meta = pd.DataFrame({
            "barcode": sub2.obs_names,
            "subclass": ct.values[keep],
            "cell_type_broad": broad.values[keep],
            "donor": sub2.obs[pick(sub2.obs.columns, ["Donor ID", "donor_id", "donor"])].values
                     if pick(sub2.obs.columns, ["Donor ID", "donor_id", "donor"]) else "NA",
            "region": region,
            "source_object": p.stem,
        }).set_index("barcode")
        parts.append((X, meta))
        ad.file.close()
        print(f"  {p.name:55s} +{X.shape[0]:6d} cells")

    if not parts:
        sys.exit(f"[{region}] no cells in any subclass object — check region name")

    counts = sp.vstack([x for x, _ in parts]).tocsc()   # cells x genes
    meta = pd.concat([m for _, m in parts])
    genes = np.array(var_ref, dtype="S")
    print(f"[{region}] total {counts.shape[0]} cells x {counts.shape[1]} genes; "
          f"nnz {counts.nnz:,} (R int32 safe < 2.1e9)")

    with h5py.File(out / "counts_raw.h5", "w") as f:
        f.create_dataset("data", data=counts.data, compression="gzip")
        f.create_dataset("indices", data=counts.indices, compression="gzip")
        f.create_dataset("indptr", data=counts.indptr, compression="gzip")
        f.attrs["shape"] = counts.shape
        f.attrs["format"] = "csc"
        f.create_dataset("barcodes", data=np.array(meta.index, dtype="S"), compression="gzip")
        f.create_dataset("gene_names", data=genes, compression="gzip")
    meta.to_csv(out / "cell_meta.csv")
    pd.DataFrame({"gene": [g.decode() for g in genes]}).to_csv(out / "gene_names.csv", index=False)
    print(meta["cell_type_broad"].value_counts().to_string())
    print(f"[{region}] wrote prep/")


if __name__ == "__main__":
    main()
