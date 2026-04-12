"""
Phase 2.9 — Prepare MTG data for CellChat / NicheNet
=====================================================
Reads SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad (34 GB),
subsets to MTG region, extracts raw counts + metadata, and
writes compact outputs that the R CellChat script can load.

Outputs (all in results/phase2/LR/prep/):
  counts_raw.h5       — sparse count matrix (cells × genes), HDF5
  cell_meta.csv       — cell barcodes + cell type labels + donor info
  gene_names.csv      — gene list matching columns of counts_raw.h5
  microglia_markers.csv — DAM marker genes for NicheNet geneset

Why Python preprocessing:
  The 34 GB h5ad is too large to load in R directly. Python/scanpy
  streams it efficiently; we write only the cells and genes needed.

Run: sbatch code/slurm/15_prep_mtg.slurm
  or interactively:
      apptainer exec containers/scenic.sif python code/phase2_LR/01_prep_mtg_for_cellchat.py
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import h5py
from pathlib import Path

PROJDIR = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
MTG_H5AD = PROJDIR / "data/raw/SEA-AD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
OUT      = PROJDIR / "results/phase2/LR/prep"
OUT.mkdir(parents=True, exist_ok=True)

# Microglial state markers for NicheNet receiver geneset
# (DAM signature vs Homeostatic — from Phase 1 DGE)
DAM_MARKERS = [
    "TREM2", "APOE", "SPP1", "LPL", "GPNMB", "CD9", "LGALS3",
    "IKZF1", "IRF8", "PPARG", "SPI1", "RUNX1", "CEBPB",
    "CST7", "CTSL", "CTSD", "FTL", "FTH1", "AXL", "MYO1E",
    "TMEM119", "P2RY12", "CX3CR1",   # homeostatic (down in DAM)
]


def main():
    print(f"Loading {MTG_H5AD.name}...")
    adata = sc.read_h5ad(MTG_H5AD, backed="r")   # memory-mapped, don't load .X yet
    print(f"  Full object: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"  obs columns: {list(adata.obs.columns[:15])}")

    # ── Detect region and cell-type columns ───────────────────────────────────
    # SEA-AD MTG h5ad column names vary slightly across releases
    region_col = next(
        (c for c in ["region_of_interest_acronym", "region", "ROIGroupCoarse",
                     "donor_region", "Region"]
         if c in adata.obs.columns), None
    )
    celltype_col = next(
        (c for c in ["subclass", "Subclass", "cell_type", "class", "Class",
                     "supertype", "Supertype"]
         if c in adata.obs.columns), None
    )
    donor_col = next(
        (c for c in ["donor_id", "Donor ID", "donor", "sample"]
         if c in adata.obs.columns), None
    )
    print(f"  Region col:   {region_col}")
    print(f"  CellType col: {celltype_col}")
    print(f"  Donor col:    {donor_col}")

    # ── Subset to MTG ─────────────────────────────────────────────────────────
    if region_col:
        mtg_mask = adata.obs[region_col].astype(str).str.upper().str.contains("MTG")
        n_mtg = mtg_mask.sum()
        print(f"  MTG cells: {n_mtg:,} / {adata.n_obs:,}")
        if n_mtg < 1000:
            print("  WARNING: very few MTG cells — using all cells as fallback")
            mtg_mask = pd.Series(True, index=adata.obs.index)
    else:
        print("  No region column found — using all cells")
        mtg_mask = pd.Series(True, index=adata.obs.index)

    # ── Subsample per broad cell type to avoid R int32 overflow ──────────────
    # R sparseMatrix uses 32-bit int indices; >~2.1B non-zeros causes NAs.
    # 1.38M cells × ~2000 nnz/cell ≈ 2.76B — overflows. Cap at 5000/type.
    MAX_CELLS_PER_TYPE = 5000
    np.random.seed(42)

    obs_full = adata[mtg_mask].obs.copy()
    # Apply coarsen mapping on the fine cell-type column to get broad types
    ct_fine = obs_full[celltype_col].astype(str) if celltype_col else pd.Series(
        "Unknown", index=obs_full.index)

    def coarsen_early(ct):
        if "Micro" in ct or "PVM" in ct:   return "Microglia"
        if "Astro" in ct:                   return "Astrocyte"
        if "Oligo" in ct and "OPC" not in ct: return "Oligodendrocyte"
        if "OPC" in ct:                     return "OPC"
        if any(x in ct for x in ["Exc","L2","L3","L4","L5","L6"]): return "ExcNeuron"
        if any(x in ct for x in ["Inh","Vip","Pvalb","Sst","Lamp5","Sncg",
                                   "PVALB","SST","VIP","LAMP5"]): return "InhNeuron"
        if any(x in ct for x in ["Endo","VLMC","Peri"]): return "Vascular"
        return ct

    broad_types = ct_fine.apply(coarsen_early)
    keep_idx = []
    for bt in broad_types.unique():
        idx = obs_full.index[broad_types == bt].tolist()
        if len(idx) > MAX_CELLS_PER_TYPE:
            idx = np.random.choice(idx, MAX_CELLS_PER_TYPE, replace=False).tolist()
        keep_idx.extend(idx)

    print(f"  Subsampled to {len(keep_idx):,} cells "
          f"(max {MAX_CELLS_PER_TYPE}/type, {broad_types.nunique()} types)")

    print(f"  Accessing subset via disk-mapping...")
    adata_mtg = adata[keep_idx]
    print(f"  Loading sparse counts for {len(keep_idx):,} cells...")

    # ── Get raw counts ─────────────────────────────────────────────────────────
    if adata_mtg.raw is not None:
        print("  Using .raw.X for counts")
        counts = adata_mtg.raw.X
        genes  = adata_mtg.raw.var_names
    else:
        print("  No .raw — using .X directly (assuming raw counts)")
        counts = adata_mtg.X
        genes  = adata_mtg.var_names

    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    else:
        counts = counts.tocsr()

    print(f"  Counts shape: {counts.shape}")
    print(f"  Count dtype:  {counts.dtype}")
    print(f"  Total non-zeros: {counts.nnz:,} (max safe for R int32: ~2.1B)")

    # ── Cell metadata ─────────────────────────────────────────────────────────
    keep_cols = list(dict.fromkeys(          # deduplicate, preserve order
        c for c in [celltype_col, donor_col, region_col,
                    "class", "Class", "subclass", "Subclass",
                    "supertype", "Supertype",
                    "Braak", "braak_stage", "APOE_genotype",
                    "age_at_death", "sex"]
        if c is not None and c in adata_mtg.obs.columns
    ))
    meta = adata_mtg.obs[keep_cols].copy()
    meta.index.name = "barcode"

    # Standardize cell type label column name for R script
    if celltype_col and celltype_col in meta.columns:
        meta["cell_type"] = meta[celltype_col].astype(str)
    else:
        meta["cell_type"] = "Unknown"

    # Aggregate broad class for CellChat grouping
    # Map Micro-PVM → Microglia; keep other subclasses as-is
    def coarsen(ct):
        ct = str(ct)
        if "Micro" in ct or "PVM" in ct:
            return "Microglia"
        if "Astro" in ct:
            return "Astrocyte"
        if "Oligo" in ct and "OPC" not in ct:
            return "Oligodendrocyte"
        if "OPC" in ct:
            return "OPC"
        if "Exc" in ct or "L2" in ct or "L3" in ct or "L4" in ct or "L5" in ct or "L6" in ct:
            return "ExcNeuron"
        if "Inh" in ct or "PVALB" in ct or "SST" in ct or "VIP" in ct or "LAMP5" in ct:
            return "InhNeuron"
        if "Endo" in ct or "VLMC" in ct or "Peri" in ct:
            return "Vascular"
        return ct
    meta["cell_type_broad"] = meta["cell_type"].apply(coarsen)

    ct_counts = meta["cell_type_broad"].value_counts()
    print("\n  Cell type composition (broad):")
    for ct, n in ct_counts.items():
        print(f"    {ct:25s} {n:7,}")

    # ── Save ──────────────────────────────────────────────────────────────────
    # counts_raw.h5 — CSC sparse matrix (R dgCMatrix native format)
    # counts is currently cells × genes (CSR). Convert to CSC so that:
    #   indptr → column pointers (length ngenes+1) — matches R sparseMatrix(p=)
    #   indices → row indices (cell indices)         — matches R sparseMatrix(i=)
    counts_csc = counts.tocsc()   # genes × cells after t() in R; here still cells×genes CSC
    h5_path = OUT / "counts_raw.h5"
    print(f"\nWriting {h5_path} (CSC format for R compatibility)...")
    with h5py.File(h5_path, "w") as f:
        f.create_dataset("data",    data=counts_csc.data,    compression="gzip")
        f.create_dataset("indices", data=counts_csc.indices, compression="gzip")
        f.create_dataset("indptr",  data=counts_csc.indptr,  compression="gzip")
        f.attrs["shape"] = counts_csc.shape   # (ncells, ngenes)
        f.attrs["format"] = "csc"
        # Store barcodes and gene names inside the h5
        f.create_dataset("barcodes",
                         data=np.array(adata_mtg.obs_names, dtype="S"),
                         compression="gzip")
        f.create_dataset("gene_names",
                         data=np.array(genes, dtype="S"),
                         compression="gzip")
    adata.file.close()

    # cell_meta.csv
    meta_path = OUT / "cell_meta.csv"
    meta.to_csv(meta_path)
    print(f"Wrote {meta_path} ({len(meta):,} rows)")

    # gene_names.csv (redundant with h5 but convenient for R)
    genes_path = OUT / "gene_names.csv"
    pd.DataFrame({"gene": list(genes)}).to_csv(genes_path, index=False)
    print(f"Wrote {genes_path} ({len(genes):,} genes)")

    # DAM marker list for NicheNet
    markers_present = [g for g in DAM_MARKERS if g in genes]
    pd.DataFrame({"gene": markers_present}).to_csv(
        OUT / "microglia_markers.csv", index=False
    )
    print(f"Wrote microglia_markers.csv ({len(markers_present)} / {len(DAM_MARKERS)} markers present)")

    print("\nDone.")


if __name__ == "__main__":
    main()
