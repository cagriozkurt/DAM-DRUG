"""
DAM-DRUG Item 5.7: CellOracle In-Silico TF Perturbation
=========================================================
Simulates transcription factor knockdown (IKZF1, IRF8, SPI1) in
the pySCENIC-derived microglial GRN and shows predicted cell state
transitions via vector field visualisation.

Pipeline:
  1. Load + prepare anndata (microglia_trajectory.h5ad, stratified 20K subsample)
  2. Build TFdict from pySCENIC adjacency matrix + validated regulon TF list
  3. Initialise Oracle, load TF data, PCA, KNN imputation
  4. Fit GRN per cluster (ridge regression, alpha=10)
  5. For each KO TF: simulate_shift → transition probs → embedding shift → plots

Inputs:
  results/phase1/trajectory/microglia_trajectory.h5ad
  results/phase2/GRN/adj_matrix_aggregated.tsv
  results/phase2/GRN/regulons_aggregated.csv

Outputs:
  results/phase5/celloracle/KO_<TF>_quiver.pdf
  results/phase5/celloracle/KO_<TF>_grid.pdf
  results/phase5/celloracle/KO_<TF>_transition_table.csv

Run:
  apptainer exec containers/scmultiomegrn.sif python code/phase5_validation/01_celloracle_perturbation.py
  OR: sbatch code/slurm/33_celloracle_perturbation.slurm
"""

import os
import re
import pickle
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import celloracle as co
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sc.settings.verbosity = 1

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT      = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
# Primary: trajectory h5ad (has UMAP + state labels); fallback: original SEA-AD h5ad
MICROGLIA_H5  = PROJECT / "results/phase1/trajectory/microglia_trajectory.h5ad"
SEAAD_H5      = PROJECT / "data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad"
ADJ_TSV       = PROJECT / "results/phase2/GRN/adj_matrix_aggregated.tsv"
REGULONS_CSV  = PROJECT / "results/phase2/GRN/regulons_aggregated.csv"
OUT           = PROJECT / "results/phase5/celloracle"
OUT.mkdir(parents=True, exist_ok=True)

# TFs to knock down
KO_TFS = ["IKZF1", "IRF8", "SPI1"]
N_CELLS = 20_000   # stratified subsample — celloracle is RAM-intensive

# Supertype → readable state name
STATE_MAP = {
    # SEA-AD trajectory h5ad uses exact names; SEA-AD original h5ad adds -SEAAD suffix
    "Micro-PVM_1":         "Homeostatic",
    "Micro-PVM_2":         "DAM",
    "Micro-PVM_2_3":       "DAM_IRM",
    "Micro-PVM_2_3-SEAAD": "DAM_IRM",
    "Micro-PVM_3":         "IRM",
    "Micro-PVM_3-SEAAD":   "IRM",
    "Micro-PVM_4":         "LAM",
    "Micro-PVM_4-SEAAD":   "LAM",
    "Micro-PVM_2_1":       "LateAD_DAM",
    "Micro-PVM_2_1-SEAAD": "LateAD_DAM",
}
STATE_COLORS = {
    "Homeostatic": "#4CAF50", "DAM": "#F44336", "DAM_IRM": "#FF9800",
    "IRM": "#2196F3",         "LAM": "#9C27B0", "LateAD_DAM": "#795548",
}


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Prepare anndata
# ─────────────────────────────────────────────────────────────────────────────
ADATA_FILE = OUT / "adata_prepared.h5ad"

if ADATA_FILE.exists():
    print(f"[Step 1] SKIP — loading {ADATA_FILE}")
    adata = sc.read_h5ad(ADATA_FILE)
else:
    print("[Step 1] Loading microglia h5ad...", flush=True)
    h5_path = MICROGLIA_H5 if MICROGLIA_H5.exists() else SEAAD_H5
    print(f"  Source: {h5_path}")
    adata = sc.read_h5ad(h5_path)
    print(f"  Raw: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Diagnose obs columns + Supertype values to catch naming mismatches
    print(f"  obs columns: {list(adata.obs.columns)}")
    for col in ["Supertype", "subtype", "cell_type", "microglial_state", "cluster"]:
        if col in adata.obs.columns:
            vals = adata.obs[col].value_counts()
            print(f"  {col} values (top 10):\n{vals.head(10).to_string()}")
            break

    # State labels
    if "Supertype" in adata.obs.columns:
        adata.obs["state"] = adata.obs["Supertype"].map(STATE_MAP).fillna("Other")
    elif "microglial_state" in adata.obs.columns:
        adata.obs["state"] = adata.obs["microglial_state"]
    else:
        raise ValueError(f"No state column found. obs columns: {list(adata.obs.columns)}")

    adata = adata[adata.obs["state"] != "Other"].copy()
    n_states = adata.obs["state"].nunique()
    print(f"  After state filter: {adata.n_obs:,} cells, {n_states} states")
    print(f"  State counts: {adata.obs['state'].value_counts().to_dict()}")

    # If only 1-2 states found and we used the trajectory h5ad, try original SEA-AD h5ad
    if n_states < 4 and h5_path == MICROGLIA_H5 and SEAAD_H5.exists():
        print(f"  Only {n_states} states in trajectory h5ad — switching to original SEA-AD h5ad...", flush=True)
        adata = sc.read_h5ad(SEAAD_H5)
        print(f"  SEA-AD raw: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        if "Supertype" in adata.obs.columns:
            adata.obs["state"] = adata.obs["Supertype"].map(STATE_MAP).fillna("Other")
        adata = adata[adata.obs["state"] != "Other"].copy()
        print(f"  After state filter: {adata.n_obs:,} cells, {adata.obs['state'].nunique()} states")
        print(f"  State counts: {adata.obs['state'].value_counts().to_dict()}")

    # Stratified subsample to N_CELLS
    if adata.n_obs > N_CELLS:
        print(f"  Stratified subsample → {N_CELLS:,} cells...", flush=True)
        state_counts = adata.obs["state"].value_counts()
        chosen = []
        rng = np.random.default_rng(42)
        for state, total in state_counts.items():
            n = max(200, int(N_CELLS * total / adata.n_obs))
            idx = adata.obs_names[adata.obs["state"] == state].tolist()
            if len(idx) > n:
                idx = rng.choice(idx, n, replace=False).tolist()
            chosen.extend(idx)
        adata = adata[chosen].copy()
        print(f"  After subsample: {adata.n_obs:,} cells")

    # Raw counts check — determine import method for celloracle
    # Priority: .raw → layers["counts"] → .X as-is (already normalized)
    use_raw = False
    if adata.raw is not None:
        print("  Using adata.raw for raw integer counts...")
        raw_adata = adata.raw.to_adata()
        shared_genes = adata.var_names.intersection(raw_adata.var_names)
        adata = adata[:, shared_genes].copy()
        adata.X = raw_adata[adata.obs_names, shared_genes].X
        use_raw = True
    elif "counts" in adata.layers:
        print("  Using adata.layers['counts'] for raw integer counts...")
        adata.X = adata.layers["counts"]
        use_raw = True
    else:
        X_sample = adata.X[:200].toarray() if sp.issparse(adata.X) else adata.X[:200]
        if np.allclose(X_sample, np.round(X_sample)):
            print("  adata.X appears to be raw integer counts — using directly.")
            use_raw = True
        else:
            print("  adata.X is normalized (non-integer). Will use import_anndata_as_normalized_count.")
    adata.uns["use_raw"] = use_raw

    # Gene filter: keep genes expressed in ≥5% of cells
    sc.pp.filter_genes(adata, min_cells=max(10, int(0.05 * adata.n_obs)))
    print(f"  After gene filter: {adata.n_vars:,} genes", flush=True)

    # UMAP must exist
    if "X_umap" not in adata.obsm:
        print("  X_umap not found — computing UMAP on normalised counts...")
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm, target_sum=1e4)
        sc.pp.log1p(adata_norm)
        sc.pp.highly_variable_genes(adata_norm, n_top_genes=2000)
        sc.pp.pca(adata_norm, n_comps=50)
        sc.pp.neighbors(adata_norm, n_pcs=50)
        sc.tl.umap(adata_norm)
        adata.obsm["X_umap"] = adata_norm.obsm["X_umap"]

    adata.write_h5ad(ADATA_FILE)
    print(f"  Saved: {ADATA_FILE}", flush=True)

print(f"\n  adata.shape: {adata.shape}")
print(f"  States: {adata.obs['state'].value_counts().to_dict()}")
print(f"  X dtype: {adata.X.dtype}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Build TFdict from pySCENIC GRN
# ─────────────────────────────────────────────────────────────────────────────
TFDICT_FILE = OUT / "tfdict.pkl"

if TFDICT_FILE.exists():
    with open(TFDICT_FILE, "rb") as f:
        TFdict = pickle.load(f)
    print(f"\n[Step 2] SKIP — TFdict loaded ({len(TFdict)} TFs)")
else:
    print("\n[Step 2] Building TFdict from pySCENIC adj_matrix...", flush=True)

    adj = pd.read_csv(ADJ_TSV, sep="\t")
    if list(adj.columns[:3]) != ["TF", "target", "importance"]:
        adj.columns = ["TF", "target", "importance"] + list(adj.columns[3:])
    print(f"  adj_matrix: {len(adj):,} edges, {adj['TF'].nunique()} TFs")

    # Get validated regulon TF set (post-RcisTarget motif filter)
    regulon_tfs = set()
    try:
        reg = pd.read_csv(REGULONS_CSV, header=0)
        # pySCENIC CSV: first column contains regulon names like "IKZF1(+)"
        for name in reg.iloc[:, 0].dropna().unique():
            tf = re.sub(r"\([^)]*\)$", "", str(name)).strip()
            if tf and tf != "nan":
                regulon_tfs.add(tf)
        if not regulon_tfs:
            # Try header row itself (some formats have TF names as column headers)
            for col in reg.columns:
                tf = re.sub(r"\([^)]*\)$", "", str(col)).strip()
                regulon_tfs.add(tf)
        print(f"  Validated regulon TFs ({len(regulon_tfs)}): {sorted(regulon_tfs)[:15]}")
    except Exception as e:
        print(f"  regulons CSV parse error: {e}. Using all TFs from adj_matrix.")
        regulon_tfs = set(adj["TF"].unique())

    # Filter adj to: validated TFs, genes present in adata
    genes = set(adata.var_names)
    adj_f = adj[
        adj["TF"].isin(regulon_tfs) &
        adj["target"].isin(genes) &
        adj["TF"].isin(genes)
    ].copy()
    print(f"  Filtered edges: {len(adj_f):,}")

    # Top 500 targets per TF by importance score
    TFdict = {}
    for tf, grp in adj_f.groupby("TF"):
        targets = grp.nlargest(500, "importance")["target"].tolist()
        if targets:
            TFdict[tf] = targets
    print(f"  TFdict: {len(TFdict)} TFs")

    # Ensure each KO TF has entries (fall back to raw adj if not in regulons)
    for tf in KO_TFS:
        if tf not in TFdict:
            fb = adj[(adj["TF"] == tf) & (adj["target"].isin(genes))]
            if not fb.empty:
                TFdict[tf] = fb.nlargest(500, "importance")["target"].tolist()
                print(f"  {tf} not in validated regulons — added from raw adj ({len(TFdict[tf])} targets)")
            else:
                print(f"  WARNING: {tf} not found in adj_matrix — KO will be skipped")

    with open(TFDICT_FILE, "wb") as f:
        pickle.dump(TFdict, f)
    print(f"  Saved: {TFDICT_FILE}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Initialise Oracle + fit GRN (done once, reused for all KOs)
# ─────────────────────────────────────────────────────────────────────────────
ORACLE_FILE = OUT / "oracle_fitted.celloracle.oracle"

print("\n[Step 3] Initialising CellOracle Oracle...", flush=True)
oracle = co.Oracle()
use_raw = adata.uns.get("use_raw", False)
if use_raw:
    print("  Importing as raw count (natural_log transform)...")
    oracle.import_anndata_as_raw_count(
        adata,
        cluster_column_name="state",
        embedding_name="X_umap",
        transform="natural_log",
    )
else:
    print("  Importing as normalized count (already log-normalized)...")
    oracle.import_anndata_as_normalized_count(
        adata,
        cluster_column_name="state",
        embedding_name="X_umap",
    )
oracle.import_TF_data(TFdict=TFdict)

if ORACLE_FILE.exists():
    print(f"  oracle_fitted file exists — still need to PCA+KNN+fit (no load API)...")
    oracle.perform_PCA()
    n_pca = min(50, oracle.pca.n_components_)
    k = max(20, min(50, adata.n_obs // 400))
    oracle.knn_imputation(n_pca_dims=n_pca, k=k, metric="euclidean", n_jobs=20)
    oracle.fit_GRN_for_simulation(alpha=10, GRN_unit="cluster", verbose_level=1)
else:
    print("[Step 4] PCA + KNN imputation...", flush=True)
    oracle.perform_PCA()
    n_pca = min(50, oracle.pca.n_components_)
    k = max(20, min(50, adata.n_obs // 400))
    print(f"  n_pca={n_pca}, k={k}")
    oracle.knn_imputation(n_pca_dims=n_pca, k=k, metric="euclidean", n_jobs=20)

    print("[Step 4] Fitting GRN per cluster (ridge alpha=10)...", flush=True)
    oracle.fit_GRN_for_simulation(alpha=10, GRN_unit="cluster", verbose_level=1)

    oracle.to_hdf5(str(ORACLE_FILE))
    print(f"  Saved fitted oracle: {ORACLE_FILE}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Step 5: Simulate KO for each target TF
# ─────────────────────────────────────────────────────────────────────────────
print("\n[Step 5] Simulating TF knockdowns...", flush=True)

cluster_colors = [STATE_COLORS.get(s, "#AAAAAA") for s in adata.obs["state"]]

for tf in KO_TFS:
    if tf not in TFdict:
        print(f"  SKIP {tf} — not in TFdict")
        continue

    print(f"\n  --- {tf} KO ---", flush=True)
    prefix = OUT / f"KO_{tf}"

    # Simulate KO: set target TF expression to 0 in all cells
    oracle.simulate_shift(
        perturb_condition={tf: 0.0},
        n_propagation=3,
    )

    # Transition probabilities
    oracle.estimate_transition_prob(
        n_neighbors=200,
        knn_random=True,
        sampled_fraction=0.5,
        calculate_randomized=True,
        n_jobs=20,
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # ── Quiver plot (raw per-cell vectors) ────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 7))
    oracle.plot_quiver(ax=ax, scale=25, s=3,
                       args={"lw": 0.3, "rasterized": True})
    ax.set_title(f"CellOracle: {tf} KO — per-cell vectors", fontsize=12)
    ax.axis("off")
    plt.tight_layout()
    fig.savefig(f"{prefix}_quiver.pdf", dpi=150, bbox_inches="tight")
    fig.savefig(f"{prefix}_quiver.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {prefix}_quiver.pdf")

    # ── Grid flow plot ────────────────────────────────────────────────────────
    try:
        oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
        oracle.calculate_mass_filter(min_mass=0.005, plot=False)
        oracle.calculate_grid_arrows(smooth=0.5, steps=(40, 40), n_neighbors=100, n_jobs=20)

        fig, ax = plt.subplots(figsize=(8, 7))
        plt.sca(ax)
        oracle.plot_simulation_flow_on_grid(ax=ax, scale=25, s=3,
                                             args={"linewidths": 0.3, "width": 0.004})
        ax.set_title(f"CellOracle: {tf} KO — grid flow", fontsize=12)
        ax.axis("off")
        plt.tight_layout()
        fig.savefig(f"{prefix}_grid.pdf", dpi=150, bbox_inches="tight")
        fig.savefig(f"{prefix}_grid.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {prefix}_grid.pdf")
    except Exception as e:
        print(f"  Grid plot failed: {e}")

    # ── Save per-cell delta summary (direction of shift by state) ─────────────
    try:
        delta_df = pd.DataFrame(
            oracle.delta_embedding,
            index=oracle.adata.obs_names,
            columns=["dUMAP1", "dUMAP2"]
        )
        delta_df["state"] = oracle.adata.obs["state"].values
        summary = delta_df.groupby("state")[["dUMAP1", "dUMAP2"]].mean()
        summary["magnitude"] = np.sqrt(summary["dUMAP1"]**2 + summary["dUMAP2"]**2)
        summary.to_csv(f"{prefix}_delta_summary.csv")
        print(f"  Mean shift by state:\n{summary.round(4).to_string()}")
    except Exception as e:
        print(f"  Delta summary failed: {e}")

print(f"\n=== CellOracle perturbation complete ===")
print(f"Outputs in: {OUT}")
