"""
WP2.2 — Per-region CellChat input prep (neurons + microglia)
===========================================================
Generalises code/phase2_LR/01_prep_mtg_for_cellchat.py to an arbitrary SEA-AD
region h5ad. Same output layout, same seed, same 5,000-cells/broad-type cap
(R dgCMatrix int32 safety), same coarsen() mapping — so each region's CellChat
run is parameter-identical to the accepted MTG analysis.

Usage (inside scenic.sif):
  python s3_prep.py <REGION> <path/to/region.h5ad>

Outputs -> results/phase7/cellchat_multiregion/<REGION>/prep/
  counts_raw.h5, cell_meta.csv, gene_names.csv
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


def coarsen(ct: str) -> str:
    ct = str(ct)
    if "Micro" in ct or "PVM" in ct:
        return "Microglia"
    if "Astro" in ct:
        return "Astrocyte"
    if "Oligo" in ct and "OPC" not in ct:
        return "Oligodendrocyte"
    if "OPC" in ct:
        return "OPC"
    if any(x in ct for x in ["Exc", "L2", "L3", "L4", "L5", "L6"]):
        return "ExcNeuron"
    if any(x in ct for x in ["Inh", "PVALB", "SST", "VIP", "LAMP5", "Pvalb",
                             "Sst", "Vip", "Lamp5", "Sncg"]):
        return "InhNeuron"
    if any(x in ct for x in ["Endo", "VLMC", "Peri"]):
        return "Vascular"
    return ct


def main():
    region = sys.argv[1]
    h5ad = Path(sys.argv[2])
    out = PROJDIR / f"results/phase7/cellchat_multiregion/{region}/prep"
    out.mkdir(parents=True, exist_ok=True)

    print(f"[{region}] loading {h5ad.name}")
    ad = sc.read_h5ad(h5ad, backed="r")
    print(f"  {ad.n_obs:,} cells x {ad.n_vars:,} genes")

    celltype_col = next((c for c in ["subclass", "Subclass", "supertype",
                                     "Supertype", "cell_type", "class", "Class"]
                         if c in ad.obs.columns), None)
    donor_col = next((c for c in ["donor_id", "Donor ID", "donor", "sample"]
                      if c in ad.obs.columns), None)
    print(f"  celltype col: {celltype_col}  donor col: {donor_col}")

    ct_fine = (ad.obs[celltype_col].astype(str) if celltype_col
               else pd.Series("Unknown", index=ad.obs.index))
    broad = ct_fine.apply(coarsen)
    # keep only sender/receiver-relevant lineages + context cells
    keep_types = {"Microglia", "InhNeuron", "ExcNeuron", "Astrocyte",
                  "Oligodendrocyte", "OPC", "Vascular"}
    rng = np.random.default_rng(SEED)
    keep_idx = []
    for bt in broad.unique():
        if bt not in keep_types:
            continue
        idx = ad.obs.index[broad == bt].to_numpy()
        if len(idx) > MAX_CELLS_PER_TYPE:
            idx = rng.choice(idx, MAX_CELLS_PER_TYPE, replace=False)
        keep_idx.extend(idx.tolist())
    print(f"  kept {len(keep_idx):,} cells across {broad[broad.isin(keep_types)].nunique()} broad types")

    sub = ad[keep_idx]
    counts = sub.raw.X if sub.raw is not None else sub.X
    genes = sub.raw.var_names if sub.raw is not None else sub.var_names
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    csc = counts.tocsc()
    print(f"  counts {csc.shape}  nnz {csc.nnz:,} (R int32 safe < 2.1e9)")

    with h5py.File(out / "counts_raw.h5", "w") as f:
        f.create_dataset("data", data=csc.data, compression="gzip")
        f.create_dataset("indices", data=csc.indices, compression="gzip")
        f.create_dataset("indptr", data=csc.indptr, compression="gzip")
        f.attrs["shape"] = csc.shape
        f.attrs["format"] = "csc"
        f.create_dataset("barcodes", data=np.array(sub.obs_names, dtype="S"), compression="gzip")
        f.create_dataset("gene_names", data=np.array(genes, dtype="S"), compression="gzip")

    meta_cols = [c for c in [celltype_col, donor_col, "Braak", "sex", "age_at_death"]
                 if c and c in sub.obs.columns]
    meta = sub.obs[meta_cols].copy()
    meta.index.name = "barcode"
    meta["cell_type"] = ct_fine.loc[sub.obs_names].values
    meta["cell_type_broad"] = broad.loc[sub.obs_names].values
    meta["region"] = region
    meta.to_csv(out / "cell_meta.csv")
    pd.DataFrame({"gene": list(genes)}).to_csv(out / "gene_names.csv", index=False)
    ad.file.close()
    print(f"  wrote prep/ for {region}")
    print(meta["cell_type_broad"].value_counts().to_string())


if __name__ == "__main__":
    main()
