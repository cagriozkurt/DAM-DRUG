"""
DAM-DRUG Phase 1 — Batch Correction (scVI primary, Harmony fallback)
=====================================================================
Input:  data/processed/microglia_combined_raw.h5ad
Output: data/processed/microglia_scVI_corrected.h5ad
        data/processed/microglia_harmony_corrected.h5ad  (fallback)
        results/phase1/batch_correction_benchmark.csv    (kBET / iLISI)

Acceptance thresholds (per Research Plan):
  kBET acceptance rate > 0.70
  ARI vs published annotations > 0.75

Run: python code/phase1_QC/03_batch_correction_scVI.py
Requires: scvi-tools>=1.0, harmonypy, scib (for benchmarking)
GPU strongly recommended for scVI at full SEA-AD scale.
"""

import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
PROC    = PROJECT / "data/processed"
RES     = PROJECT / "results/phase1"
sc.settings.figdir = str(RES)

# ==============================================================================
# Load microglia subset
# ==============================================================================
print("Loading microglia subset...")
mg = sc.read_h5ad(PROC / "microglia_combined_raw.h5ad")
print(f"  {mg.n_obs:,} microglia × {mg.n_vars:,} genes")

# Confirm raw counts are in layers["counts"]
assert "counts" in mg.layers, "Raw counts layer missing — re-run 02_QC_microglia_isolation.py"


# ==============================================================================
# Strategy A — scVI (primary; Luecken 2022 winner at >500K cells)
# ==============================================================================
print("\n=== Strategy A: scVI ===")
import scvi

scvi.settings.seed = 42

# Setup: batch = dataset × donor (fine-grained batch correction)
# Covariates: sex, brain_region as categorical; pmi, age as continuous
categorical_covs = [c for c in ["sex"] if c in mg.obs.columns]
continuous_covs  = [c for c in ["pmi", "age_at_death"] if c in mg.obs.columns]

scvi.model.SCVI.setup_anndata(
    mg,
    layer="counts",
    batch_key="dataset",
    categorical_covariate_keys=categorical_covs if categorical_covs else None,
    continuous_covariate_keys=continuous_covs if continuous_covs else None,
)

model = scvi.model.SCVI(
    mg,
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb",        # negative binomial — appropriate for snRNA-seq
    dispersion="gene-batch",     # estimate per-gene dispersion per batch
)

# Train — GPU auto-detected; falls back to CPU if not available
use_gpu = True
try:
    import torch
    use_gpu = torch.cuda.is_available()
except ImportError:
    use_gpu = False

print(f"  Training scVI (GPU={use_gpu}, max_epochs=400)...")
model.train(
    max_epochs=400,
    use_gpu=use_gpu,
    early_stopping=True,
    early_stopping_patience=20,
    plan_kwargs={"lr": 1e-3},
)

mg.obsm["X_scVI"] = model.get_latent_representation()
model.save(str(PROC / "scVI_model_microglia"), overwrite=True)

# UMAP on scVI latent space
sc.pp.neighbors(mg, use_rep="X_scVI", n_neighbors=30, key_added="neighbors_scVI")
sc.tl.umap(mg, neighbors_key="neighbors_scVI")
sc.tl.leiden(mg, resolution=0.8, neighbors_key="neighbors_scVI", key_added="leiden_scVI")

sc.pl.umap(mg, color=["dataset", "leiden_scVI", "P2RY12", "IRF8", "APOE"],
           ncols=3, save="_microglia_scVI.png")

mg.write_h5ad(PROC / "microglia_scVI_corrected.h5ad", compression="gzip")
print("  scVI embedding saved.")


# ==============================================================================
# Strategy B — Harmony (fallback; faster, CPU-only)
# ==============================================================================
print("\n=== Strategy B: Harmony (fallback / cross-check) ===")
import harmonypy as hm

# Run PCA first (on log-normalised HVG expression)
mg_harmony = mg.copy()
sc.pp.scale(mg_harmony, max_value=10)
sc.tl.pca(mg_harmony, svd_solver="arpack", n_comps=50, use_highly_variable=True)

ho = hm.run_harmony(
    mg_harmony.obsm["X_pca"],
    mg_harmony.obs,
    vars_use=["dataset", "donor_id"],
    max_iter_harmony=30,
    random_state=42,
)
mg_harmony.obsm["X_pca_harmony"] = ho.Z_corr.T

sc.pp.neighbors(mg_harmony, use_rep="X_pca_harmony", n_neighbors=30, key_added="neighbors_harmony")
sc.tl.umap(mg_harmony, neighbors_key="neighbors_harmony")
sc.tl.leiden(mg_harmony, resolution=0.8, neighbors_key="neighbors_harmony", key_added="leiden_harmony")

sc.pl.umap(mg_harmony, color=["dataset", "leiden_harmony", "P2RY12", "IRF8", "APOE"],
           ncols=3, save="_microglia_harmony.png")

mg_harmony.write_h5ad(PROC / "microglia_harmony_corrected.h5ad", compression="gzip")
print("  Harmony embedding saved.")


# ==============================================================================
# Benchmark: kBET and iLISI (requires scib)
# ==============================================================================
print("\n=== Benchmarking batch correction ===")
try:
    import scib

    metrics = {}

    for method, adata, rep, clust in [
        ("scVI",    mg,         "X_scVI",       "leiden_scVI"),
        ("Harmony", mg_harmony, "X_pca_harmony","leiden_harmony"),
    ]:
        print(f"  Evaluating {method}...")
        kbet = scib.me.kBET(
            adata,
            batch_key="dataset",
            label_key=clust,
            type_="embed",
            embed=rep,
        )
        ilisi = scib.me.ilisi_graph(
            adata,
            batch_key="dataset",
            type_="embed",
            use_rep=rep,
        )
        metrics[method] = {"kBET": kbet, "iLISI": ilisi}
        print(f"    kBET={kbet:.3f}  iLISI={ilisi:.3f}")

    bench_df = pd.DataFrame(metrics).T
    bench_df.to_csv(RES / "batch_correction_benchmark.csv")
    print(f"\n  Benchmark saved to results/phase1/batch_correction_benchmark.csv")

    # Decision: use scVI if kBET >= 0.70, else Harmony
    scvi_kbet = metrics["scVI"]["kBET"]
    if scvi_kbet >= 0.70:
        print(f"\n  DECISION: scVI accepted (kBET={scvi_kbet:.3f} >= 0.70) → proceed with scVI embedding")
    else:
        print(f"\n  WARNING: scVI kBET={scvi_kbet:.3f} < 0.70 threshold")
        print("  Consider: (1) adding more covariates, (2) increasing n_latent, (3) using Harmony as fallback")

except ImportError:
    print("  scib not installed — skipping kBET/iLISI benchmark")
    print("  Install: pip install scib-metrics")

print("\nDone. Next step: run 04_microglial_substate_annotation.py")
