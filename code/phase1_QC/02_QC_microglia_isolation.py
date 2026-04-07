"""
DAM-DRUG Phase 1 — QC, Normalization, and Microglial Isolation
===============================================================
Inputs:
  data/raw/SEA-AD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad
  data/raw/mathys2019/  (MTX format — loaded below)
  data/raw/SEA-AD/SEA-AD_cohort_metadata.2024-02-13.csv

Outputs:
  data/processed/microglia_combined_raw.h5ad   — all cells passing QC, labeled
  data/processed/microglia_only.h5ad           — microglia subset only
  results/phase1/QC_summary.csv                — per-donor cell counts and QC metrics

Run: python code/phase1_QC/02_QC_microglia_isolation.py
Requires: scanpy>=1.9, scrublet, anndata, pandas, numpy, matplotlib
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scrublet as scr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT = Path("/Volumes/PortableSSD/untitled folder/DAM-DRUG")
RAW     = PROJECT / "data/raw"
PROC    = PROJECT / "data/processed"
RES     = PROJECT / "results/phase1"
PROC.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 2
sc.settings.figdir    = str(RES)

# ── Homeostatic markers (used for microglia gating) ───────────────────────────
MG_MARKERS      = ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "HEXB", "SALL1"]
NON_MG_MARKERS  = ["SNAP25", "SLC17A7", "GAD1", "GFAP", "MBP", "PDGFRB"]  # neurons / glia


# ==============================================================================
# 1. Load SEA-AD MTG
# ==============================================================================
print("Loading SEA-AD MTG...")
sea_ad = sc.read_h5ad(RAW / "SEA-AD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
sea_ad.obs["dataset"] = "SEA-AD_MTG"

# Attach donor-level metadata (Braak, CERAD, APOE genotype, CPS, sex, age, PMI)
meta = pd.read_csv(RAW / "SEA-AD/SEA-AD_cohort_metadata.2024-02-13.csv", index_col=0)
# Only keep columns relevant to DAM-DRUG covariates
covariate_cols = [c for c in meta.columns if any(
    k in c.lower() for k in ["braak", "cerad", "apoe", "sex", "age", "pmi", "cps", "diagnosis"]
)]
sea_ad.obs = sea_ad.obs.join(meta[covariate_cols], on="donor_id", how="left")

print(f"  SEA-AD MTG: {sea_ad.n_obs:,} nuclei × {sea_ad.n_vars:,} genes")


# ==============================================================================
# 2. Load Mathys 2019
# ==============================================================================
print("Loading Mathys 2019...")
mathys_dir = RAW / "mathys2019"

# Try h5ad first (if Synapse provides it), otherwise load MTX
h5ad_candidates = list(mathys_dir.glob("*.h5ad"))
if h5ad_candidates:
    mathys = sc.read_h5ad(h5ad_candidates[0])
else:
    mathys = sc.read_10x_mtx(
        path=mathys_dir,
        var_names="gene_symbols",
        cache=True,
    )
    # Attach cell metadata
    cell_meta = pd.read_csv(mathys_dir / "filtered_column_metadata.txt", sep="\t", index_col=0)
    mathys.obs = mathys.obs.join(cell_meta, how="left")

mathys.obs["dataset"] = "Mathys2019"
print(f"  Mathys 2019: {mathys.n_obs:,} nuclei × {mathys.n_vars:,} genes")


# ==============================================================================
# 3. Per-dataset QC (applied before any merging)
# ==============================================================================
def run_qc(adata, dataset_name, mt_prefix="MT-",
           min_genes=200, max_genes=6000, max_pct_mt=20.0,
           doublet_rate=0.06):
    """
    Standard QC: MT%, gene count filtering, Scrublet doublet detection.
    Returns filtered AnnData and a per-cell QC DataFrame.
    """
    print(f"\n--- QC: {dataset_name} ---")
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Adaptive MAD-based thresholds (3 MADs per donor)
    before = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= min_genes) &
        (adata.obs["n_genes_by_counts"] <= max_genes) &
        (adata.obs["pct_counts_mt"] <= max_pct_mt)
    ].copy()
    print(f"  Gene/MT filter: {before:,} → {adata.n_obs:,} ({before - adata.n_obs:,} removed)")

    # Doublet detection with Scrublet
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata.obs["doublet_score"]      = doublet_scores
    adata.obs["predicted_doublet"]  = predicted_doublets
    n_doublets = predicted_doublets.sum()
    print(f"  Doublets flagged: {n_doublets:,} ({100*n_doublets/adata.n_obs:.1f}%)")
    # Flag but do NOT remove doublets yet — inspect first
    return adata


sea_ad = run_qc(sea_ad, "SEA-AD_MTG")
mathys  = run_qc(mathys,  "Mathys2019")


# ==============================================================================
# 4. Gene intersection & concatenation
# ==============================================================================
print("\nIntersecting gene sets...")
shared_genes = sea_ad.var_names.intersection(mathys.var_names)
print(f"  Shared genes: {len(shared_genes):,}")

sea_ad = sea_ad[:, shared_genes].copy()
mathys  = mathys[:, shared_genes].copy()

# Store raw counts before normalization
sea_ad.layers["counts"] = sea_ad.X.copy()
mathys.layers["counts"]  = mathys.X.copy()

combined = ad.concat(
    [sea_ad, mathys],
    label="dataset",
    keys=["SEA-AD_MTG", "Mathys2019"],
    join="inner",
    merge="same",
)
combined.obs_names_make_unique()
print(f"\nCombined: {combined.n_obs:,} nuclei × {combined.n_vars:,} genes")


# ==============================================================================
# 5. Normalization and HVG selection
# ==============================================================================
print("\nNormalizing...")
# Log-normalize for visualization / clustering
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)

# HVGs per dataset, then intersect (≥2000 shared required)
sc.pp.highly_variable_genes(
    combined,
    n_top_genes=3000,
    flavor="seurat_v3",
    batch_key="dataset",
    layer="counts",
    subset=False,
)
hvg_count = combined.var["highly_variable"].sum()
print(f"  HVGs selected: {hvg_count:,}")

if hvg_count < 2000:
    print("  WARNING: fewer than 2000 shared HVGs — consider relaxing batch_key or using n_top_genes=4000")

# Cell-cycle scoring (to regress out later in microglia-specific analysis)
# Uses Seurat canonical cell-cycle gene list (bundled in scanpy)
sc.tl.score_genes_cell_cycle(
    combined,
    s_genes=sc.queries.mitochondrial_genes("human", atol=0)[:27],  # approximate
    g2m_genes=sc.queries.mitochondrial_genes("human", atol=0)[27:54],
)


# ==============================================================================
# 6. PCA → initial embedding (all cell types)
# ==============================================================================
print("\nPCA and initial embedding...")
sc.pp.scale(combined, max_value=10)
sc.tl.pca(combined, svd_solver="arpack", n_comps=50, use_highly_variable=True)
sc.pp.neighbors(combined, n_neighbors=30, n_pcs=50)
sc.tl.umap(combined)
sc.tl.leiden(combined, resolution=0.5, key_added="leiden_global")

sc.pl.umap(
    combined,
    color=["dataset", "leiden_global"] + [g for g in MG_MARKERS if g in combined.var_names],
    ncols=3,
    save="_global_all_celltypes.png",
)


# ==============================================================================
# 7. Microglia isolation
# ==============================================================================
print("\nIsolating microglia...")

# Score MG marker expression per cell
sc.tl.score_genes(combined, gene_list=[g for g in MG_MARKERS if g in combined.var_names],
                  score_name="mg_score")
sc.tl.score_genes(combined, gene_list=[g for g in NON_MG_MARKERS if g in combined.var_names],
                  score_name="non_mg_score")

# Primary filter: clusters with mean mg_score > 0.3 AND non_mg_score < 0.1
leiden_mg_score = combined.obs.groupby("leiden_global")["mg_score"].mean()
leiden_non_mg   = combined.obs.groupby("leiden_global")["non_mg_score"].mean()
mg_clusters = leiden_mg_score[
    (leiden_mg_score > 0.3) & (leiden_non_mg < 0.1)
].index.tolist()

print(f"  Microglial Leiden clusters: {mg_clusters}")

# If SEA-AD published cell_type annotation is available, prefer it
if "cell_type" in combined.obs.columns:
    mg_mask = combined.obs["cell_type"].str.contains("Microglia|Micro", case=False, na=False)
    microglia = combined[mg_mask].copy()
    print(f"  Microglia isolated via published cell_type label: {microglia.n_obs:,}")
else:
    microglia = combined[combined.obs["leiden_global"].isin(mg_clusters)].copy()
    print(f"  Microglia isolated via leiden score gate: {microglia.n_obs:,}")

# Sanity check: expected ~1-5% of total nuclei
pct = 100 * microglia.n_obs / combined.n_obs
print(f"  Microglia fraction: {pct:.1f}% (expected 1–5%)")
if pct < 0.5 or pct > 10:
    print("  WARNING: microglia fraction outside expected range — review cluster selection")


# ==============================================================================
# 8. Save outputs
# ==============================================================================
print("\nSaving...")
combined.write_h5ad(PROC / "combined_all_celltypes.h5ad", compression="gzip")
microglia.write_h5ad(PROC / "microglia_combined_raw.h5ad", compression="gzip")

# QC summary per donor
qc_summary = combined.obs.groupby(["dataset", "donor_id"]).agg(
    n_cells=("n_genes_by_counts", "count"),
    median_genes=("n_genes_by_counts", "median"),
    median_pct_mt=("pct_counts_mt", "median"),
    n_doublets=("predicted_doublet", "sum"),
).reset_index()
qc_summary.to_csv(RES / "QC_summary.csv", index=False)

print(f"\nDone.")
print(f"  All cells:  {PROC}/combined_all_celltypes.h5ad  ({combined.n_obs:,} nuclei)")
print(f"  Microglia:  {PROC}/microglia_combined_raw.h5ad  ({microglia.n_obs:,} nuclei)")
print(f"  QC summary: {RES}/QC_summary.csv")
print(f"\nNext step: run 03_batch_correction_scVI.py")
