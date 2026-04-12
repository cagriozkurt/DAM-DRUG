"""
DAM-DRUG Phase 1 — Differential Gene Expression: Microglial States
===================================================================
Input:  data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad
Output: results/phase1/DGE/
  - dge_{contrast}.csv              — full DEG tables (logFC, p-adj, pct)
  - dge_{contrast}_top50.csv        — top 50 per contrast for reporting
  - dge_apoe4_vs_apoe3_{state}.csv  — APOE4-stratified DEGs per state
  - dge_heatmap_{contrast}.png      — expression heatmaps
  - volcano_{contrast}.png          — volcano plots

Contrasts (per Research Plan):
  1. Micro-PVM_4-SEAAD (LAM) vs Micro-PVM_1 (Homeostatic)   [primary DAM]
  2. Micro-PVM_2 (DAM)        vs Micro-PVM_1 (Homeostatic)   [DAM vs HM]
  3. Micro-PVM_3-SEAAD (IRM)  vs Micro-PVM_1 (Homeostatic)   [IRM vs HM]
  4. Braak High (V-VI) vs Low (0-II) — within Micro-PVM_2     [AD progression]
  5. APOE4/4 vs APOE3/3 — across all microglia               [APOE effect]

State→DAM-DRUG mapping:
  Micro-PVM_1           → Homeostatic
  Micro-PVM_2           → DAM
  Micro-PVM_2_3-SEAAD   → DAM/IRM transition
  Micro-PVM_2_1-SEAAD   → Late-AD DAM
  Micro-PVM_3-SEAAD     → IRM
  Micro-PVM_4-SEAAD     → LAM/LDAM  ← primary therapeutic target state

Run: python code/phase1_QC/03_DGE_microglial_states.py
Requires: scanpy>=1.9, anndata, pandas, numpy, matplotlib, seaborn
GPU not required (rank_genes_groups uses CPU Wilcoxon).
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
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
RAW     = PROJECT / "data/raw/SEA-AD"
RES     = PROJECT / "results/phase1/DGE"
RES.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(RES)

H5AD = Path(os.environ.get(
    "SEA_AD_H5AD",
    str(RAW / "SEA-AD_Microglia_multi-regional_final-nuclei.2025-07-24.h5ad")
))

# Column names (confirmed from exploration)
SUPERTYPE_COL = "Supertype"
APOE_COL      = "APOE Genotype"
BRAAK_COL     = "Braak"
DX_COL        = "Overall AD neuropathological Change"
REGION_COL    = "Brain Region"
DONOR_COL     = "Donor ID"
SEX_COL       = "Sex"

# State mapping for readability
STATE_MAP = {
    "Micro-PVM_1":         "Homeostatic",
    "Micro-PVM_2":         "DAM",
    "Micro-PVM_2_3-SEAAD": "DAM-IRM",
    "Micro-PVM_2_1-SEAAD": "LateAD-DAM",
    "Micro-PVM_3-SEAAD":   "IRM",
    "Micro-PVM_4-SEAAD":   "LAM",
}

# TF targets of interest
TF_TARGETS = ["SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
              "RELB", "BHLHE40", "BHLHE41"]

# Drug-relevant genes (targets + markers)
DRUG_GENES = TF_TARGETS + [
    "TREM2", "APOE", "SPP1", "LGALS3", "CD9", "GPNMB", "LPL", "LIPA",
    "P2RY12", "CX3CR1", "TMEM119", "CSF1R", "HEXB",
    "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "IRF7",
    "CD28", "MKI67",
]


# ==============================================================================
# 1. Load object
# ==============================================================================
print("Loading SEA-AD Microglia pre-release object...")
mg = sc.read_h5ad(H5AD)
print(f"  {mg.n_obs:,} nuclei × {mg.n_vars:,} genes")

# Rename columns for scanpy compatibility (spaces → no issue in obs, but let's alias)
mg.obs["apoe_gt"]   = mg.obs[APOE_COL].astype(str)
mg.obs["braak"]     = mg.obs[BRAAK_COL].astype(str)
mg.obs["dx"]        = mg.obs[DX_COL].astype(str)
mg.obs["region"]    = mg.obs[REGION_COL].astype(str)
mg.obs["donor"]     = mg.obs[DONOR_COL].astype(str)
mg.obs["sex"]       = mg.obs[SEX_COL].astype(str)
mg.obs["state"]     = mg.obs[SUPERTYPE_COL].map(STATE_MAP).fillna(mg.obs[SUPERTYPE_COL])

# UMI counts are in layers["UMIs"] — use for DGE (raw counts preferred)
# X appears to be log-normalized (log1p key in uns). Keep layers["UMIs"] as raw.
print(f"  APOE genotypes: {mg.obs['apoe_gt'].value_counts().to_dict()}")
print(f"  Braak stages:   {sorted(mg.obs['braak'].unique())}")
print(f"  AD dx:          {mg.obs['dx'].value_counts().to_dict()}")


# ==============================================================================
# 2. Filter: keep only microglia (exclude Lymphocyte + Monocyte)
# ==============================================================================
mg_only = mg[mg.obs[SUPERTYPE_COL].isin(list(STATE_MAP.keys()))].copy()
print(f"\n  After excluding Lymphocyte + Monocyte: {mg_only.n_obs:,} cells")

# Set raw counts for DGE (Wilcoxon on log-norm X is standard in scanpy)
# X already log-normalized (confirmed by uns['log1p'])
# For rank_genes_groups, scanpy expects log-normalized X — confirmed.


# ==============================================================================
# Helper: run DGE, save results
# ==============================================================================
def run_dge(adata, groupby, groups, reference, contrast_name, n_genes=500):
    """
    Run Wilcoxon rank-sum DGE. Save full table and top-50 table.
    Returns DataFrame with all results.
    """
    print(f"\n  [{contrast_name}] {groups} vs {reference} (n={adata.n_obs:,})")

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups if isinstance(groups, list) else [groups],
        reference=reference,
        method="wilcoxon",
        key_added=f"rank_{contrast_name}",
        pts=True,         # percent expressing
        n_genes=n_genes,
    )

    results = []
    for g in (groups if isinstance(groups, list) else [groups]):
        df = sc.get.rank_genes_groups_df(
            adata,
            group=g,
            key=f"rank_{contrast_name}",
            pval_cutoff=None,   # get all genes
            log2fc_min=None,
        )
        df["group"]    = g
        df["contrast"] = contrast_name
        results.append(df)

    df_all = pd.concat(results, ignore_index=True)

    # Save full results
    df_all.to_csv(RES / f"dge_{contrast_name}.csv", index=False)

    # Save top 50 by abs(log2FC), padj < 0.05
    top50 = (
        df_all[df_all["pvals_adj"] < 0.05]
        .assign(abs_logfc=lambda x: x["logfoldchanges"].abs())
        .sort_values("abs_logfc", ascending=False)
        .head(50)
    )
    top50.to_csv(RES / f"dge_{contrast_name}_top50.csv", index=False)

    n_sig = (df_all["pvals_adj"] < 0.05).sum()
    n_up  = ((df_all["pvals_adj"] < 0.05) & (df_all["logfoldchanges"] > 0.5)).sum()
    n_dn  = ((df_all["pvals_adj"] < 0.05) & (df_all["logfoldchanges"] < -0.5)).sum()
    print(f"    Sig DEGs (padj<0.05): {n_sig:,}  |  UP(logFC>0.5): {n_up:,}  |  DN(logFC<-0.5): {n_dn:,}")

    # TF targets in significant DEGs
    tf_cols = ["names", "logfoldchanges", "pvals_adj"] + [
        c for c in ["pts", "pts_rest"] if c in df_all.columns
    ]
    tf_sig = df_all[df_all["names"].isin(TF_TARGETS) & (df_all["pvals_adj"] < 0.05)][
        tf_cols
    ].sort_values("logfoldchanges", ascending=False)
    if not tf_sig.empty:
        print(f"    TF targets significant:")
        print(tf_sig.to_string(index=False))

    return df_all


def volcano_plot(df, contrast_name, logfc_col="logfoldchanges", padj_col="pvals_adj",
                 name_col="names", highlight=None, fc_thresh=0.5, padj_thresh=0.05):
    """Simple volcano plot; highlight = list of gene names to label."""
    df = df.copy()
    fig, ax = plt.subplots(figsize=(7, 6))

    # Color by significance
    sig_up = (df[padj_col] < padj_thresh) & (df[logfc_col] > fc_thresh)
    sig_dn = (df[padj_col] < padj_thresh) & (df[logfc_col] < -fc_thresh)
    ns     = ~(sig_up | sig_dn)

    df["neg_log10_padj"] = -np.log10(df[padj_col].clip(lower=1e-300))

    ax.scatter(df.loc[ns,     logfc_col], df.loc[ns,     "neg_log10_padj"],
               s=2, c="#bbbbbb", alpha=0.4, rasterized=True)
    ax.scatter(df.loc[sig_dn, logfc_col], df.loc[sig_dn, "neg_log10_padj"],
               s=4, c="#2166ac", alpha=0.6, rasterized=True)
    ax.scatter(df.loc[sig_up, logfc_col], df.loc[sig_up, "neg_log10_padj"],
               s=4, c="#d6604d", alpha=0.6, rasterized=True)

    # Label highlighted genes
    if highlight:
        for _, row in df[df[name_col].isin(highlight)].iterrows():
            ax.annotate(
                row[name_col],
                xy=(row[logfc_col], row["neg_log10_padj"]),
                fontsize=7, ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-", lw=0.5),
                xytext=(row[logfc_col], row["neg_log10_padj"] + 1.5),
            )

    ax.axvline(fc_thresh,  ls="--", lw=0.8, c="grey")
    ax.axvline(-fc_thresh, ls="--", lw=0.8, c="grey")
    ax.axhline(-np.log10(padj_thresh), ls="--", lw=0.8, c="grey")
    ax.set_xlabel("log2 Fold Change", fontsize=11)
    ax.set_ylabel("-log10 adj. p-value", fontsize=11)
    ax.set_title(contrast_name.replace("_", " "), fontsize=11)
    ax.text(0.02, 0.97, f"UP={sig_up.sum():,}  DN={sig_dn.sum():,}",
            transform=ax.transAxes, va="top", fontsize=9, color="#555555")
    plt.tight_layout()
    fig.savefig(RES / f"volcano_{contrast_name}.png", dpi=150)
    plt.close()


# ==============================================================================
# 3. Contrast 1: LAM (Micro-PVM_4-SEAAD) vs Homeostatic (Micro-PVM_1)
# ==============================================================================
print("\n" + "=" * 70)
print("Contrast 1: LAM vs Homeostatic")

c1_cells = mg_only[mg_only.obs[SUPERTYPE_COL].isin(["Micro-PVM_4-SEAAD", "Micro-PVM_1"])].copy()
df_c1 = run_dge(c1_cells, SUPERTYPE_COL,
                "Micro-PVM_4-SEAAD", "Micro-PVM_1",
                "LAM_vs_Homeostatic")
volcano_plot(df_c1, "LAM_vs_Homeostatic", highlight=DRUG_GENES)


# ==============================================================================
# 4. Contrast 2: DAM (Micro-PVM_2) vs Homeostatic (Micro-PVM_1)
# ==============================================================================
print("\n" + "=" * 70)
print("Contrast 2: DAM vs Homeostatic")

c2_cells = mg_only[mg_only.obs[SUPERTYPE_COL].isin(["Micro-PVM_2", "Micro-PVM_1"])].copy()
df_c2 = run_dge(c2_cells, SUPERTYPE_COL,
                "Micro-PVM_2", "Micro-PVM_1",
                "DAM_vs_Homeostatic")
volcano_plot(df_c2, "DAM_vs_Homeostatic", highlight=DRUG_GENES)


# ==============================================================================
# 5. Contrast 3: IRM (Micro-PVM_3-SEAAD) vs Homeostatic (Micro-PVM_1)
# ==============================================================================
print("\n" + "=" * 70)
print("Contrast 3: IRM vs Homeostatic")

c3_cells = mg_only[mg_only.obs[SUPERTYPE_COL].isin(["Micro-PVM_3-SEAAD", "Micro-PVM_1"])].copy()
df_c3 = run_dge(c3_cells, SUPERTYPE_COL,
                "Micro-PVM_3-SEAAD", "Micro-PVM_1",
                "IRM_vs_Homeostatic")
volcano_plot(df_c3, "IRM_vs_Homeostatic", highlight=DRUG_GENES)


# ==============================================================================
# 6. Contrast 4: Braak High vs Low within DAM cells
# ==============================================================================
print("\n" + "=" * 70)
print("Contrast 4: Braak V-VI vs Braak 0-II (within DAM = Micro-PVM_2)")

# Map Braak stages to high/low
braak_high = {"Braak V", "Braak VI"}
braak_low  = {"Braak 0", "Braak II"}

dam_cells = mg_only[mg_only.obs[SUPERTYPE_COL] == "Micro-PVM_2"].copy()
dam_cells = dam_cells[dam_cells.obs["braak"].isin(braak_high | braak_low)].copy()
dam_cells.obs["braak_group"] = dam_cells.obs["braak"].apply(
    lambda x: "BraakHighV_VI" if x in braak_high else "BraakLow0_II"
)
print(f"  DAM cells: {dam_cells.n_obs:,}  |  High={dam_cells.obs['braak_group'].eq('BraakHighV_VI').sum():,}  Low={dam_cells.obs['braak_group'].eq('BraakLow0_II').sum():,}")

df_c4 = run_dge(dam_cells, "braak_group",
                "BraakHighV_VI", "BraakLow0_II",
                "DAM_BraakHigh_vs_Low")
volcano_plot(df_c4, "DAM_BraakHigh_vs_Low", highlight=DRUG_GENES)


# ==============================================================================
# 7. Contrast 5: APOE4/4 vs APOE3/3 (all microglia)
# ==============================================================================
print("\n" + "=" * 70)
print("Contrast 5: APOE4/4 vs APOE3/3 (all microglial states)")

apoe_cells = mg_only[mg_only.obs["apoe_gt"].isin(["4/4", "3/3"])].copy()
print(f"  APOE4/4: {(apoe_cells.obs['apoe_gt']=='4/4').sum():,}  |  APOE3/3: {(apoe_cells.obs['apoe_gt']=='3/3').sum():,}")

if (apoe_cells.obs["apoe_gt"] == "4/4").sum() >= 50 and \
   (apoe_cells.obs["apoe_gt"] == "3/3").sum() >= 50:
    df_c5 = run_dge(apoe_cells, "apoe_gt",
                    "4/4", "3/3",
                    "APOE4_vs_APOE3_allMicroglia")
    volcano_plot(df_c5, "APOE4_vs_APOE3_allMicroglia", highlight=DRUG_GENES)
else:
    print("  Insufficient cells for APOE4/4 vs 3/3 — skipping")

# APOE4 within LAM state specifically
lam_apoe = mg_only[
    (mg_only.obs[SUPERTYPE_COL] == "Micro-PVM_4-SEAAD") &
    (mg_only.obs["apoe_gt"].isin(["4/4", "3/3", "3/4"]))
].copy()
print(f"\n  APOE within LAM: {lam_apoe.obs['apoe_gt'].value_counts().to_dict()}")
if (lam_apoe.obs["apoe_gt"] == "4/4").sum() >= 30:
    df_c5b = run_dge(lam_apoe, "apoe_gt",
                     "4/4", "3/3",
                     "APOE4_vs_APOE3_LAM")
    volcano_plot(df_c5b, "APOE4_vs_APOE3_LAM", highlight=DRUG_GENES)
else:
    print("  Insufficient APOE4/4 cells in LAM — skipping stratified analysis")


# ==============================================================================
# 8. Summary heatmap — TF targets + key markers across all states
# ==============================================================================
print("\n" + "=" * 70)
print("Generating summary heatmap — TF targets across states...")

# Mean log-normalized expression of drug-relevant genes per state
drug_present = [g for g in DRUG_GENES if g in mg_only.var_names]
state_means  = {}

for state, label in STATE_MAP.items():
    cells = mg_only[mg_only.obs[SUPERTYPE_COL] == state]
    if cells.n_obs == 0:
        continue
    gene_idx = [list(mg_only.var_names).index(g) for g in drug_present]
    X_sub = cells.X[:, gene_idx]
    if sp.issparse(X_sub):
        X_sub = X_sub.toarray()
    state_means[f"{label}\n(n={cells.n_obs:,})"] = X_sub.mean(axis=0)

heatmap_df = pd.DataFrame(state_means, index=drug_present)

# Z-score per gene for visualization
heatmap_z = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0).div(
    heatmap_df.std(axis=1).clip(lower=0.01), axis=0
)

fig, ax = plt.subplots(figsize=(10, 12))
sns.heatmap(
    heatmap_z,
    ax=ax,
    cmap="RdBu_r",
    center=0,
    vmin=-2, vmax=2,
    linewidths=0.3,
    linecolor="#cccccc",
    cbar_kws={"label": "Z-score (mean log-norm expr)"},
    annot=heatmap_df.round(2),
    fmt=".2f",
    annot_kws={"size": 7},
)
ax.set_title("DAM-DRUG Target Genes: Mean Expression by Microglial State\n(SEA-AD 240K nuclei)", fontsize=11)
ax.set_xlabel("")
ax.set_ylabel("")
plt.xticks(rotation=30, ha="right", fontsize=9)
plt.yticks(fontsize=8)
plt.tight_layout()
fig.savefig(RES / "heatmap_TF_markers_by_state.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {RES}/heatmap_TF_markers_by_state.png")

# Save the mean expression table
heatmap_df.to_csv(RES / "mean_expr_by_state.csv")


# ==============================================================================
# 9. Print consolidated TF target summary
# ==============================================================================
print("\n" + "=" * 70)
print("TF TARGET SUMMARY ACROSS ALL CONTRASTS")
print("=" * 70)

contrast_files = {
    "LAM_vs_HM":    RES / "dge_LAM_vs_Homeostatic.csv",
    "DAM_vs_HM":    RES / "dge_DAM_vs_Homeostatic.csv",
    "IRM_vs_HM":    RES / "dge_IRM_vs_Homeostatic.csv",
    "BraakH_vs_L":  RES / "dge_DAM_BraakHigh_vs_Low.csv",
}

tf_summary = {}
for contrast, fpath in contrast_files.items():
    if not fpath.exists():
        continue
    df = pd.read_csv(fpath)
    tf_rows = df[df["names"].isin(TF_TARGETS)][["names", "logfoldchanges", "pvals_adj"]].copy()
    tf_rows = tf_rows.set_index("names")
    tf_rows.columns = [f"{contrast}_logFC", f"{contrast}_padj"]
    tf_summary[contrast] = tf_rows

if tf_summary:
    summary_df = pd.concat(tf_summary.values(), axis=1).reindex(TF_TARGETS)
    print(summary_df.round(3).to_string())
    summary_df.to_csv(RES / "TF_target_summary.csv")

print(f"\nDone. Output: {RES}")
print("Next step: run 01_pySCENIC_GRN.py (GRN/regulon inference)")
