"""
DAM-DRUG Phase 1 — Explore SEA-AD Microglia Pre-Release Object
================================================================
Input:  data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad
Output: results/phase1/explore/  (figures + CSV summaries)

Goals:
  1. Validate supertype labels vs. expected 5-state annotation scheme
  2. Inspect APOE genotype distribution across supertypes/donors
  3. Check Braak stage correlation with microglial state
  4. Flag PU.1-low (SPI1-low) / CD28+ cluster (reactive subpopulation)
  5. Check LDAM marker expression (LGALS3, SPP1, CD9, GPNMB)
  6. Inspect batch structure (brain regions, donors)
  7. Confirm scVI latent space + UMAP are present and usable
  8. Identify any annotation gaps before DGE/GRN analysis

Run: python code/phase1_QC/02_explore_microglia_object.py
Requires: scanpy>=1.9, anndata, pandas, numpy, matplotlib, seaborn
"""

import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
RAW     = PROJECT / "data/raw/SEA-AD"
RES     = PROJECT / "results/phase1/explore"
RES.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 2
sc.settings.figdir    = str(RES)

H5AD = Path(os.environ.get(
    "SEA_AD_H5AD",
    str(RAW / "SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad")
))

# ==============================================================================
# 1. Load and basic inspection
# ==============================================================================
print("=" * 70)
print("Loading SEA-AD Microglia pre-release object...")
mg = sc.read_h5ad(H5AD)
print(f"  Shape: {mg.n_obs:,} nuclei × {mg.n_vars:,} genes")
print(f"\n--- obs columns ({len(mg.obs.columns)}) ---")
for col in sorted(mg.obs.columns):
    n_unique = mg.obs[col].nunique()
    sample   = mg.obs[col].dropna().iloc[0] if mg.obs[col].notna().any() else "NA"
    print(f"  {col:<40}  n_unique={n_unique:<6}  sample={sample}")

print(f"\n--- var columns ({len(mg.var.columns)}) ---")
for col in sorted(mg.var.columns):
    print(f"  {col}")

print(f"\n--- obsm keys ---")
for k in mg.obsm.keys():
    print(f"  {k}: shape {mg.obsm[k].shape}")

print(f"\n--- layers ---")
for k in mg.layers.keys():
    print(f"  {k}")

print(f"\n--- uns keys ---")
for k in mg.uns.keys():
    print(f"  {k}")


# ==============================================================================
# 2. Supertype label validation
# ==============================================================================
print("\n" + "=" * 70)
print("Supertype label distribution:")

# Identify the supertype/subtype column
supertype_col = None
for candidate in ["Supertype", "supertype", "cell_type", "subtype", "cluster",
                  "leiden", "Class", "Subclass"]:
    if candidate in mg.obs.columns:
        supertype_col = candidate
        break

if supertype_col:
    st_counts = mg.obs[supertype_col].value_counts()
    print(f"\n  Column: '{supertype_col}'")
    print(st_counts.to_string())

    # Save
    st_counts.reset_index().rename(
        columns={"index": supertype_col, supertype_col: "n_cells"}
    ).to_csv(RES / "supertype_counts.csv", index=False)
else:
    print("  WARNING: No supertype column found. Check obs columns above.")

# Expected states for cross-reference
EXPECTED_STATES = {
    "Homeostatic": ["Micro-PVM_1", "homeostatic", "P2RY12-high"],
    "DAM":         ["Micro-PVM_2", "DAM", "TREM2-high", "SPP1-high"],
    "IRM":         ["Micro-PVM_3", "IRM", "interferon", "ISG"],
    "LAM":         ["Micro-PVM_4", "LAM", "LDAM", "lipid"],
    "Reactive":    ["reactive", "proliferating", "CD28", "SPI1-low"],
}
print("\n  Expected state mapping (Micro-PVM_1–4 scheme):")
for state, aliases in EXPECTED_STATES.items():
    print(f"    {state}: {aliases}")


# ==============================================================================
# 3. APOE genotype distribution
# ==============================================================================
print("\n" + "=" * 70)
print("APOE genotype distribution:")

apoe_col = None
for candidate in ["APOE_genotype", "apoe_genotype", "APOE", "apoe",
                  "apoe4_dosage", "APOE4", "apoe4"]:
    if candidate in mg.obs.columns:
        apoe_col = candidate
        break

if apoe_col:
    apoe_counts = mg.obs[apoe_col].value_counts()
    print(f"\n  Column: '{apoe_col}'")
    print(apoe_counts.to_string())

    if supertype_col:
        apoe_st = pd.crosstab(mg.obs[supertype_col], mg.obs[apoe_col], normalize="index") * 100
        print(f"\n  APOE genotype % per supertype:")
        print(apoe_st.round(1).to_string())
        apoe_st.to_csv(RES / "apoe_by_supertype_pct.csv")
else:
    print("  WARNING: No APOE genotype column found.")
    apoe_candidates = [c for c in mg.obs.columns if "apoe" in c.lower() or "e4" in c.lower()]
    print(f"  APOE-related columns: {apoe_candidates}")


# ==============================================================================
# 4. Braak stage correlation
# ==============================================================================
print("\n" + "=" * 70)
print("Braak stage distribution:")

braak_col = None
for candidate in ["Braak_stage", "braak_stage", "Braak", "braak",
                  "neuropathology_braak_stage"]:
    if candidate in mg.obs.columns:
        braak_col = candidate
        break

if braak_col:
    braak_counts = mg.obs[braak_col].value_counts().sort_index()
    print(f"\n  Column: '{braak_col}'")
    print(braak_counts.to_string())

    if supertype_col:
        braak_st = pd.crosstab(mg.obs[supertype_col], mg.obs[braak_col], normalize="index") * 100
        print(f"\n  Braak stage % per supertype:")
        print(braak_st.round(1).to_string())
        braak_st.to_csv(RES / "braak_by_supertype_pct.csv")
else:
    print("  WARNING: No Braak stage column found.")
    braak_candidates = [c for c in mg.obs.columns if "braak" in c.lower()]
    print(f"  Braak-related columns: {braak_candidates}")


# ==============================================================================
# 5. Diagnosis / AD status
# ==============================================================================
print("\n" + "=" * 70)
print("Diagnosis / disease status:")

dx_col = None
for candidate in ["Diagnosis", "diagnosis", "disease", "AD_status",
                  "Overall_AD_neuropathological_Change", "overall_ad_neuropath"]:
    if candidate in mg.obs.columns:
        dx_col = candidate
        break

if dx_col:
    print(f"\n  Column: '{dx_col}'")
    print(mg.obs[dx_col].value_counts().to_string())
else:
    dx_candidates = [c for c in mg.obs.columns if any(
        k in c.lower() for k in ["diag", "disease", "pathol", "ad_"]
    )]
    print(f"  Diagnosis-related columns: {dx_candidates}")


# ==============================================================================
# 6. Brain region distribution
# ==============================================================================
print("\n" + "=" * 70)
print("Brain region distribution:")

region_col = None
for candidate in ["brain_region", "region", "tissue", "anatomical_region",
                  "tissue_section", "Region"]:
    if candidate in mg.obs.columns:
        region_col = candidate
        break

if region_col:
    print(f"\n  Column: '{region_col}'")
    print(mg.obs[region_col].value_counts().to_string())
else:
    region_candidates = [c for c in mg.obs.columns if "region" in c.lower()
                         or "tissue" in c.lower() or "area" in c.lower()]
    print(f"  Region-related columns: {region_candidates}")


# ==============================================================================
# 7. Donor statistics
# ==============================================================================
print("\n" + "=" * 70)
donor_col = None
for candidate in ["donor_id", "donor", "subject_id", "Donor_ID", "individual"]:
    if candidate in mg.obs.columns:
        donor_col = candidate
        break

if donor_col:
    n_donors = mg.obs[donor_col].nunique()
    cells_per_donor = mg.obs[donor_col].value_counts()
    print(f"Donor statistics:")
    print(f"  N donors: {n_donors}")
    print(f"  Cells/donor: median={cells_per_donor.median():.0f}, "
          f"min={cells_per_donor.min()}, max={cells_per_donor.max()}")

    # Sex distribution
    sex_col = next((c for c in ["sex", "Sex", "gender"] if c in mg.obs.columns), None)
    if sex_col:
        print(f"\n  Sex distribution:")
        print(f"  {mg.obs[[donor_col, sex_col]].drop_duplicates()[sex_col].value_counts().to_string()}")


# ==============================================================================
# 8. Marker gene expression check
# ==============================================================================
print("\n" + "=" * 70)
print("Marker gene expression check:")

# Core microglia markers
MG_HOMEOSTATIC = ["P2RY12", "TMEM119", "CX3CR1", "HEXB", "SALL1", "CSF1R"]
DAM_MARKERS    = ["TREM2", "SPP1", "CD9", "GPNMB", "LGALS3", "LPL", "APOE"]
LAM_MARKERS    = ["LGALS3", "SPP1", "CD9", "GPNMB", "LIPA", "LPL"]
IRM_MARKERS    = ["IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "IRF7"]
SPI1_LOW       = ["SPI1"]  # PU.1 — look for low-expressing cluster
BHLHE_TARGETS  = ["BHLHE40", "BHLHE41"]
TF_TARGETS     = ["SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
                  "RELB", "BHLHE40", "BHLHE41"]
REACTIVE       = ["CD28", "MKI67", "TOP2A"]  # proliferating / reactive

all_markers = list(set(
    MG_HOMEOSTATIC + DAM_MARKERS + LAM_MARKERS + IRM_MARKERS +
    TF_TARGETS + REACTIVE
))

present = [g for g in all_markers if g in mg.var_names]
missing = [g for g in all_markers if g not in mg.var_names]
print(f"  Markers present: {len(present)}/{len(all_markers)}")
if missing:
    print(f"  Missing: {missing}")

# Mean expression per marker (across all cells)
# Use .X (log-normalized) if available, else layers
print("\n  Mean log-normalized expression (top markers):")
for gene in sorted(present):
    gene_idx = list(mg.var_names).index(gene)
    if sp.issparse(mg.X):
        expr = mg.X[:, gene_idx].toarray().flatten()
    else:
        expr = mg.X[:, gene_idx]
    pct_expr = 100 * (expr > 0).mean()
    print(f"    {gene:<12}  mean={expr.mean():.3f}  pct_expr={pct_expr:.1f}%")


# ==============================================================================
# 9. Confirm scVI latent space and UMAP
# ==============================================================================
print("\n" + "=" * 70)
print("Embedding validation:")

scvi_key = next((k for k in mg.obsm.keys() if "scvi" in k.lower() or "scVI" in k), None)
umap_key = next((k for k in mg.obsm.keys() if "umap" in k.lower()), None)
pca_key  = next((k for k in mg.obsm.keys() if "pca"  in k.lower()), None)

print(f"  scVI latent: {scvi_key} {'✓' if scvi_key else '✗ MISSING'}")
print(f"  UMAP:        {umap_key} {'✓' if umap_key else '✗ MISSING'}")
print(f"  PCA:         {pca_key}  {'✓' if pca_key  else '(not required)'}")

if umap_key and supertype_col:
    print(f"\n  Plotting UMAP colored by {supertype_col}...")
    sc.pl.embedding(mg, basis=umap_key, color=[supertype_col],
                    save=f"_supertype_umap.png", show=False, frameon=False)
    print(f"  Saved: {RES}/{umap_key}_supertype_umap.png")

# If APOE and Braak available, add those to UMAP
extra_colors = []
if apoe_col:  extra_colors.append(apoe_col)
if braak_col: extra_colors.append(braak_col)
if region_col: extra_colors.append(region_col)
marker_for_plot = [g for g in ["SPI1", "TREM2", "P2RY12", "SPP1",
                                "BHLHE40", "BHLHE41"] if g in present]

if umap_key and (extra_colors or marker_for_plot):
    color_list = extra_colors + marker_for_plot
    sc.pl.embedding(mg, basis=umap_key, color=color_list,
                    ncols=3, save="_metadata_markers_umap.png",
                    show=False, frameon=False)
    print("  Saved: results/phase1/explore/umap_metadata_markers_umap.png")


# ==============================================================================
# 10. PU.1-low cluster detection (SPI1-low reactive microglia)
# ==============================================================================
print("\n" + "=" * 70)
print("SPI1 (PU.1) expression distribution:")

if "SPI1" in mg.var_names:
    spi1_idx = list(mg.var_names).index("SPI1")
    if sp.issparse(mg.X):
        spi1_expr = mg.X[:, spi1_idx].toarray().flatten()
    else:
        spi1_expr = mg.X[:, spi1_idx]

    # Define SPI1-low as bottom 10th percentile
    spi1_low_thresh = np.percentile(spi1_expr[spi1_expr > 0], 10) if (spi1_expr > 0).any() else 0
    spi1_low_mask   = spi1_expr < spi1_low_thresh
    n_low = spi1_low_mask.sum()
    print(f"  SPI1-low (< {spi1_low_thresh:.3f}, bottom 10th pct of expressing cells): {n_low:,} cells ({100*n_low/mg.n_obs:.1f}%)")

    if supertype_col:
        spi1_low_by_st = mg.obs.loc[spi1_low_mask, supertype_col].value_counts()
        print(f"\n  SPI1-low cells by supertype:")
        print(spi1_low_by_st.to_string())

    mg.obs["SPI1_low"] = spi1_low_mask


# ==============================================================================
# 11. CD28+ reactive cluster check
# ==============================================================================
print("\n" + "=" * 70)
print("CD28+ reactive microglia check:")

if "CD28" in mg.var_names:
    cd28_idx = list(mg.var_names).index("CD28")
    if sp.issparse(mg.X):
        cd28_expr = mg.X[:, cd28_idx].toarray().flatten()
    else:
        cd28_expr = mg.X[:, cd28_idx]

    n_cd28 = (cd28_expr > 0).sum()
    print(f"  CD28+ cells: {n_cd28:,} ({100*n_cd28/mg.n_obs:.1f}%)")
    if supertype_col:
        cd28_pos = mg.obs.loc[cd28_expr > 0, supertype_col].value_counts()
        print(f"\n  CD28+ cells by supertype:")
        print(cd28_pos.to_string())
else:
    print("  CD28 not in var_names — gene not detected in this panel")


# ==============================================================================
# 12. LDAM marker co-expression score
# ==============================================================================
print("\n" + "=" * 70)
print("LDAM marker scoring:")

ldam_present = [g for g in LAM_MARKERS if g in mg.var_names]
if len(ldam_present) >= 3:
    sc.tl.score_genes(mg, gene_list=ldam_present, score_name="LDAM_score")
    if supertype_col:
        ldam_by_st = mg.obs.groupby(supertype_col)["LDAM_score"].agg(["mean", "median"])
        print(f"\n  LDAM score by supertype:")
        print(ldam_by_st.round(3).to_string())
        ldam_by_st.to_csv(RES / "ldam_score_by_supertype.csv")
else:
    print(f"  Insufficient LDAM markers present ({ldam_present}) — skipping")


# ==============================================================================
# 13. Summary report
# ==============================================================================
print("\n" + "=" * 70)
print("EXPLORATION SUMMARY")
print("=" * 70)
print(f"  Object:       {H5AD.name}")
print(f"  Cells:        {mg.n_obs:,}")
print(f"  Genes:        {mg.n_vars:,}")
print(f"  Supertypes:   {supertype_col} ({mg.obs[supertype_col].nunique() if supertype_col else 'N/A'} types)")
print(f"  APOE col:     {apoe_col}")
print(f"  Braak col:    {braak_col}")
print(f"  Diagnosis:    {dx_col}")
print(f"  Region col:   {region_col}")
print(f"  Donor col:    {donor_col} ({mg.obs[donor_col].nunique() if donor_col else 'N/A'} donors)")
print(f"  scVI latent:  {scvi_key}")
print(f"  UMAP:         {umap_key}")
print(f"  TF markers:   {[g for g in TF_TARGETS if g in mg.var_names]}")
print(f"  Missing TFs:  {[g for g in TF_TARGETS if g not in mg.var_names]}")

# Annotation gap assessment
print("\n  ANNOTATION GAP ASSESSMENT:")
gaps = []
if not apoe_col:
    gaps.append("APOE genotype — critical covariate for APOE4 stratification")
if not braak_col:
    gaps.append("Braak stage — needed for AD severity correlation")
if not dx_col:
    gaps.append("Diagnosis/AD status — needed for case/control DGE")
if not scvi_key:
    gaps.append("scVI latent space — must re-run batch correction")
if not umap_key:
    gaps.append("UMAP — must re-run embedding")
missing_tfs = [g for g in ["SPI1", "IRF8", "PPARG", "RUNX1"] if g not in mg.var_names]
if missing_tfs:
    gaps.append(f"Critical TFs missing from gene panel: {missing_tfs}")

if gaps:
    for g in gaps:
        print(f"    ✗ {g}")
else:
    print("    ✓ No critical gaps — proceed to DGE and GRN analysis")

print(f"\n  Output directory: {RES}")
print("\nNext step: run 03_DGE_microglial_states.py")
