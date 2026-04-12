#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=pseudobulk_deseq2
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 1 — Pseudobulk DESeq2 (donor-level replicates, sex + age covariates)
#
# Replaces cell-level Wilcoxon DGE (03_DGE_microglial_states.py) with the
# methodologically correct pseudobulk approach:
#   1. Aggregate raw UMI counts per (donor × microglial state)
#   2. Run pydeseq2 with design ~ group + sex + age_z
#   3. Output mirrors 03_DGE_microglial_states.py file naming
#      (dge_pb_{contrast}.csv) so downstream scripts can use either.
#
# Contrasts (same as Wilcoxon):
#   LAM_vs_Homeostatic
#   DAM_vs_Homeostatic
#   IRM_vs_Homeostatic
#   DAM_BraakHigh_vs_Low    (Braak V-VI vs 0-II, within DAM cells)
#   APOE4_vs_APOE3_allMicroglia
#   APOE4_vs_APOE3_LAM
#
# Input:  data/raw/SEA-AD/SEA-AD_Microglia_multi-regional_final-nuclei.*.h5ad
# Output: results/phase1/DGE_pseudobulk/
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

OUT_DIR=$PROJDIR/results/phase1/DGE_pseudobulk
mkdir -p "$PROJDIR/logs" "$OUT_DIR"

export DAM_DRUG_DIR=$PROJDIR

echo "Start: $(date)"
echo "Node:  $SLURMD_NODENAME"

apptainer exec \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  conda run -n scanpy_env python - << 'PYEOF'
import os, warnings, sys
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJDIR = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
H5AD    = next((PROJDIR / "data/raw/SEA-AD").glob(
                "SEA-AD_Microglia_multi-regional_final-nuclei.*.h5ad"))
OUT_DIR = PROJDIR / "results/phase1/DGE_pseudobulk"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Loading {H5AD.name} ...")
import anndata as ad
mg_full = ad.read_h5ad(H5AD)
print(f"  {mg_full.n_obs:,} nuclei × {mg_full.n_vars:,} genes")

# ── Column names ──────────────────────────────────────────────────────────────
SUPERTYPE_COL = "Supertype"
DONOR_COL     = "Donor ID"
SEX_COL       = "Sex"
BRAAK_COL     = "Braak"
APOE_COL      = "APOE Genotype"

# Find age column (several possible names in SEA-AD releases)
AGE_COL = next(
    (c for c in mg_full.obs.columns
     if "age" in c.lower() and "death" in c.lower()),
    next((c for c in mg_full.obs.columns if c.lower().startswith("age")), None)
)
print(f"  Age column: {AGE_COL}")

STATE_MAP = {
    "Micro-PVM_1":         "Homeostatic",
    "Micro-PVM_2":         "DAM",
    "Micro-PVM_2_3-SEAAD": "DAM-IRM",
    "Micro-PVM_2_1-SEAAD": "LateAD-DAM",
    "Micro-PVM_3-SEAAD":   "IRM",
    "Micro-PVM_4-SEAAD":   "LAM",
}
TF_TARGETS = ["SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
              "RELB", "BHLHE40", "BHLHE41"]

# ── Subset to true microglia, alias columns ───────────────────────────────────
mg = mg_full[mg_full.obs[SUPERTYPE_COL].isin(STATE_MAP)].copy()
mg.obs["state"]  = mg.obs[SUPERTYPE_COL].map(STATE_MAP)
mg.obs["donor"]  = mg.obs[DONOR_COL].astype(str)
mg.obs["sex"]    = mg.obs[SEX_COL].astype(str)
mg.obs["braak"]  = mg.obs[BRAAK_COL].astype(str)
mg.obs["apoe"]   = mg.obs[APOE_COL].astype(str)
if AGE_COL:
    mg.obs["age_raw"] = pd.to_numeric(mg.obs[AGE_COL], errors="coerce")
print(f"  Microglia only: {mg.n_obs:,} cells  |  {mg.obs['donor'].nunique()} donors")

# ── Raw count layer ───────────────────────────────────────────────────────────
if "UMIs" in mg.layers:
    LAYER = "UMIs"
elif "counts" in mg.layers:
    LAYER = "counts"
else:
    # Fall back to X if already integer counts
    LAYER = None
    print("  WARNING: no 'UMIs' or 'counts' layer — using .X (may be log-normalized)")
print(f"  Using layer: {LAYER or 'X'}")


# ── Pseudobulk helper ─────────────────────────────────────────────────────────
def pseudobulk(adata, group_col, min_cells=10):
    """
    Aggregate raw UMI counts per (donor, group_col).
    Returns (counts_df  [samples × genes, integer],
             meta_df    [samples × covariates]).
    Drops samples with < min_cells cells.
    """
    counts_rows, meta_rows = [], []

    for (donor, group), cell_idx in adata.obs.groupby(["donor", group_col]).groups.items():
        if len(cell_idx) < min_cells:
            continue
        cells = adata[cell_idx]
        X = cells.layers[LAYER] if LAYER else cells.X
        if sp.issparse(X):
            X = np.array(X.sum(axis=0)).flatten()
        else:
            X = X.sum(axis=0).flatten()
        counts_rows.append(X.astype(int))

        ref_obs = adata.obs.loc[cell_idx[0]]
        row = {
            "sample_id": f"{donor}__{group}",
            "donor":     donor,
            "group":     group,
            "sex":       ref_obs["sex"],
            "n_cells":   len(cell_idx),
        }
        if AGE_COL:
            row["age_raw"] = ref_obs["age_raw"]
        meta_rows.append(row)

    counts_df = pd.DataFrame(
        np.vstack(counts_rows),
        columns=adata.var_names,
        index=[r["sample_id"] for r in meta_rows],
    )
    meta_df = pd.DataFrame(meta_rows).set_index("sample_id")

    # Standardize age
    if AGE_COL and "age_raw" in meta_df.columns:
        mu, sd = meta_df["age_raw"].mean(), meta_df["age_raw"].std()
        meta_df["age_z"] = (meta_df["age_raw"] - mu) / sd if sd > 0 else 0.0
    else:
        meta_df["age_z"] = 0.0

    # Encode sex as binary int (DESeq2 needs numeric or consistent strings)
    meta_df["sex_bin"] = meta_df["sex"].str.upper().str.startswith("F").astype(int)

    return counts_df, meta_df


def run_deseq2(counts_df, meta_df, test_col, case, ctrl, contrast_name,
               extra_factors=None):
    """
    Run pydeseq2 on pseudobulk counts.
    design ~ test_col + sex_bin [+ age_z if variance > 0]
    """
    # Keep only samples belonging to case / ctrl groups
    keep = meta_df[test_col].isin([case, ctrl])
    c_sub = counts_df.loc[keep]
    m_sub = meta_df.loc[keep].copy()

    n_case = (m_sub[test_col] == case).sum()
    n_ctrl = (m_sub[test_col] == ctrl).sum()
    print(f"\n  [{contrast_name}]  case={case} (n={n_case})  ctrl={ctrl} (n={n_ctrl})")

    if n_case < 2 or n_ctrl < 2:
        print("    SKIP — fewer than 2 samples in one group")
        return None

    # Build design factors
    design_factors = [test_col]
    if m_sub["sex_bin"].nunique() > 1:
        design_factors.append("sex_bin")
    if m_sub["age_z"].std() > 0.01 and m_sub["age_z"].nunique() > 2:
        design_factors.append("age_z")
    print(f"    Design factors: {design_factors}")

    # Remove genes with 0 counts in all samples (speeds up fitting)
    keep_genes = c_sub.sum(axis=0) > 0
    c_sub = c_sub.loc[:, keep_genes]

    try:
        dds = DeseqDataSet(
            counts=c_sub,
            metadata=m_sub,
            design_factors=design_factors,
            quiet=True,
        )
        dds.deseq2()

        stat = DeseqStats(dds, contrast=[test_col, case, ctrl], quiet=True)
        stat.summary()
        res = stat.results_df.reset_index().rename(columns={"index": "names"})

        # Add pct expressing (fraction of pseudobulk samples with count > 0)
        case_idx  = m_sub[m_sub[test_col] == case].index
        ctrl_idx  = m_sub[m_sub[test_col] == ctrl].index
        res["pts"]      = (c_sub.loc[case_idx] > 0).mean(axis=0).reindex(res["names"]).values
        res["pts_rest"] = (c_sub.loc[ctrl_idx] > 0).mean(axis=0).reindex(res["names"]).values
        res["n_samples_case"] = n_case
        res["n_samples_ctrl"] = n_ctrl
        res["contrast"] = contrast_name

        res.to_csv(OUT_DIR / f"dge_pb_{contrast_name}.csv", index=False)

        n_sig = res["padj"].lt(0.05).sum()
        n_up  = (res["padj"].lt(0.05) & res["log2FoldChange"].gt(0.5)).sum()
        n_dn  = (res["padj"].lt(0.05) & res["log2FoldChange"].lt(-0.5)).sum()
        print(f"    Sig (padj<0.05): {n_sig}  UP(logFC>0.5): {n_up}  DN(logFC<-0.5): {n_dn}")

        # TF targets
        tf_res = res[res["names"].isin(TF_TARGETS)][
            ["names", "log2FoldChange", "padj"]
        ].sort_values("log2FoldChange", ascending=False)
        if not tf_res.empty:
            print("    TF targets:")
            print(tf_res.to_string(index=False))

        return res

    except Exception as e:
        print(f"    ERROR: {e}")
        return None


# ── Contrast 1-3: State vs Homeostatic ───────────────────────────────────────
print("\n" + "="*70)
print("Building pseudobulk matrix (all microglia × donor × state) ...")

mg.obs["pb_group"] = mg.obs["state"]
counts_all, meta_all = pseudobulk(mg, "pb_group", min_cells=10)
print(f"  Pseudobulk matrix: {counts_all.shape[0]} samples  "
      f"({meta_all['group'].value_counts().to_dict()})")

for case_state in ["LAM", "DAM", "IRM"]:
    run_deseq2(counts_all, meta_all, "group",
               case_state, "Homeostatic",
               f"{case_state}_vs_Homeostatic")

# ── Contrast 4: Braak High vs Low within DAM ─────────────────────────────────
print("\n" + "="*70)
print("Contrast 4: Braak High (V-VI) vs Low (0-II) within DAM ...")

braak_high = {"Braak V", "Braak VI"}
braak_low  = {"Braak 0", "Braak II"}
dam = mg[(mg.obs["state"] == "DAM") &
         (mg.obs["braak"].isin(braak_high | braak_low))].copy()
dam.obs["braak_group"] = dam.obs["braak"].apply(
    lambda x: "BraakHigh" if x in braak_high else "BraakLow"
)
counts_b, meta_b = pseudobulk(dam, "braak_group", min_cells=10)
print(f"  Pseudobulk: {counts_b.shape[0]} samples")
run_deseq2(counts_b, meta_b, "group",
           "BraakHigh", "BraakLow",
           "DAM_BraakHigh_vs_Low")

# ── Contrast 5-6: APOE4/4 vs APOE3/3 ────────────────────────────────────────
print("\n" + "="*70)
print("Contrasts 5-6: APOE4/4 vs APOE3/3 ...")

for label, subset in [("allMicroglia", mg),
                       ("LAM",          mg[mg.obs["state"] == "LAM"])]:
    apoe_sub = subset[subset.obs["apoe"].isin(["4/4", "3/3"])].copy()
    counts_a, meta_a = pseudobulk(apoe_sub, "apoe", min_cells=10)
    print(f"  [{label}] Pseudobulk: {counts_a.shape[0]} samples")
    run_deseq2(counts_a, meta_a, "group",
               "4/4", "3/3",
               f"APOE4_vs_APOE3_{label}")

# ── Consolidated TF summary ───────────────────────────────────────────────────
print("\n" + "="*70)
print("TF TARGET SUMMARY ACROSS ALL PSEUDOBULK CONTRASTS")
print("="*70)

tf_rows = []
for csv in sorted(OUT_DIR.glob("dge_pb_*.csv")):
    df = pd.read_csv(csv)
    contrast = csv.stem.replace("dge_pb_", "")
    tf_df = df[df["names"].isin(TF_TARGETS)][
        ["names", "log2FoldChange", "padj"]
    ].copy()
    tf_df["contrast"] = contrast
    tf_rows.append(tf_df)

if tf_rows:
    summary = pd.concat(tf_rows, ignore_index=True)
    pivot = summary.pivot_table(
        index="names", columns="contrast",
        values=["log2FoldChange", "padj"], aggfunc="first"
    )
    print(pivot.round(3).to_string())
    summary.to_csv(OUT_DIR / "TF_pseudobulk_summary.csv", index=False)

print(f"\nDone. Output: {OUT_DIR}")
PYEOF

echo ""
echo "Finished: $(date)"
echo "Output: $OUT_DIR"
echo ""
echo "Key output files:"
echo "  dge_pb_LAM_vs_Homeostatic.csv"
echo "  dge_pb_DAM_vs_Homeostatic.csv"
echo "  dge_pb_IRM_vs_Homeostatic.csv"
echo "  dge_pb_DAM_BraakHigh_vs_Low.csv"
echo "  dge_pb_APOE4_vs_APOE3_allMicroglia.csv"
echo "  dge_pb_APOE4_vs_APOE3_LAM.csv"
echo "  TF_pseudobulk_summary.csv"
