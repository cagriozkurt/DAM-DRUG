"""
DAM-DRUG Table 1 — Cohort Demographics
=======================================
Reads microglia_trajectory.h5ad obs (no counts loaded) and produces
a per-donor demographics table grouped by cognitive status (AD vs Control).

Outputs:
  results/tables/table1_cohort_demographics.csv
  results/tables/table1_cohort_demographics.md

Run on TRUBA (needs h5ad):
  apptainer exec --cleanenv --bind /arf:/arf \\
    ~/containers/scanpy-env.sif \\
    conda run -n scenic python code/tables/make_table1_cohort.py
"""

import os
import numpy as np
import pandas as pd
import scipy.stats as stats
from pathlib import Path

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
H5AD    = PROJECT / "results/phase1/trajectory/microglia_trajectory.h5ad"
OUT     = PROJECT / "results/tables"
OUT.mkdir(parents=True, exist_ok=True)

# ── Load obs only (backed mode avoids loading 3 GB count matrix) ──────────
import scanpy as sc
print(f"Reading obs from {H5AD.name}...")
adata = sc.read_h5ad(H5AD, backed="r")
obs   = adata.obs.copy()
print(f"  {len(obs):,} cells, {obs['Donor ID'].nunique()} donors")
print(f"  Available columns: {list(obs.columns)}")

# ── Per-donor metadata (one row per donor) ────────────────────────────────
# Columns of interest — try multiple name variants for robustness
COL = {
    "donor":    "Donor ID",
    "sex":      next((c for c in ["Sex", "sex", "Gender", "gender"] if c in obs.columns), None),
    "age":      next((c for c in ["Age at Death", "age_at_death", "Age"] if c in obs.columns), None),
    "braak":    next((c for c in ["Braak", "Braak stage", "braak_stage"] if c in obs.columns), None),
    "apoe":     next((c for c in ["APOE Genotype", "apoe_genotype", "APOE4 Status"] if c in obs.columns), None),
    "cogstatus":next((c for c in ["Cognitive Status", "cognitive_status", "Overall AD neuropathological Change"]
                      if c in obs.columns), None),
    "pmi":      next((c for c in ["PMI", "pmi"] if c in obs.columns), None),
    "neurotyp": next((c for c in ["Neurotypical reference", "neurotypical"] if c in obs.columns), None),
}
print("\nColumn mapping:")
for k, v in COL.items():
    print(f"  {k:12s} → {v}")

# Cell counts per donor
cell_counts = obs.groupby("Donor ID", observed=True).size().rename("n_cells")

# One row per donor — take first value (demographic fields are constant per donor)
donor_cols = [v for v in COL.values() if v is not None]
donors = obs.groupby("Donor ID", observed=True)[donor_cols].first().copy()
donors = donors.join(cell_counts)

print(f"\nPer-donor table: {len(donors)} donors × {donors.shape[1]} cols")

# ── Define AD vs Control groups ───────────────────────────────────────────
# Use Cognitive Status if available; fall back to Neurotypical reference
if COL["cogstatus"]:
    cog = donors[COL["cogstatus"]].astype(str).str.strip()
    print(f"\nCognitive status values: {sorted(cog.unique())}")
    # SEA-AD uses: "Dementia", "No Dementia", "MCI"
    ad_mask  = cog.str.contains("Dementia", case=False, na=False) & \
               ~cog.str.contains("No Dementia", case=False, na=False)
    ctrl_mask = cog.str.contains("No Dementia|Control|Normal", case=False, na=False)
elif COL["neurotyp"]:
    nt = donors[COL["neurotyp"]].astype(str).str.strip()
    ctrl_mask = nt.str.lower().isin(["true", "yes", "1"])
    ad_mask   = ~ctrl_mask
else:
    raise ValueError("No cognitive status or neurotypical column found — check obs columns above")

donors["group"] = "Other"
donors.loc[ad_mask,   "group"] = "AD/Dementia"
donors.loc[ctrl_mask, "group"] = "Control"
print(f"\nGroup counts: {donors['group'].value_counts().to_dict()}")

# ── Helper functions ──────────────────────────────────────────────────────
def fmt_mean_sd(series):
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0: return "—"
    return f"{s.mean():.1f} ± {s.std():.1f}"

def fmt_n_pct(series, value):
    n_total = len(series.dropna())
    n_match = (series.astype(str).str.strip() == str(value)).sum()
    if n_total == 0: return "—"
    return f"{n_match} ({100*n_match/n_total:.0f}%)"

def pval_fmt(p):
    if np.isnan(p): return "—"
    if p < 0.001:   return "<0.001"
    return f"{p:.3f}"

def ttest_p(a, b):
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2: return np.nan
    return stats.ttest_ind(a, b).pvalue

def chi2_p(series, groups):
    ct = pd.crosstab(series, groups)
    if ct.shape[0] < 2 or ct.shape[1] < 2: return np.nan
    return stats.chi2_contingency(ct, correction=False).pvalue

ad   = donors[donors["group"] == "AD/Dementia"]
ctrl = donors[donors["group"] == "Control"]

# ── Build Table 1 ─────────────────────────────────────────────────────────
rows = []

def add_row(feature, ad_val, ctrl_val, pval="—"):
    rows.append({"Feature": feature,
                 f"AD/Dementia (n={len(ad)})": ad_val,
                 f"Control (n={len(ctrl)})": ctrl_val,
                 "p-value": pval})

add_row("N donors", str(len(ad)), str(len(ctrl)))
add_row("N cells (total)", f"{ad['n_cells'].sum():,}", f"{ctrl['n_cells'].sum():,}")
add_row("N cells per donor (mean±SD)",
        fmt_mean_sd(ad["n_cells"]), fmt_mean_sd(ctrl["n_cells"]),
        pval_fmt(ttest_p(ad["n_cells"], ctrl["n_cells"])))

if COL["sex"]:
    sex_col = donors[COL["sex"]].astype(str).str.strip()
    sex_vals = [v for v in sex_col.unique() if v not in ("nan", "")]
    # Detect female label
    f_label = next((v for v in sex_vals if v.lower() in ("female","f","woman")), sex_vals[0] if sex_vals else None)
    if f_label:
        add_row(f"Sex — {f_label} (%)",
                fmt_n_pct(ad[COL["sex"]], f_label),
                fmt_n_pct(ctrl[COL["sex"]], f_label),
                pval_fmt(chi2_p(sex_col[donors["group"].isin(["AD/Dementia","Control"])],
                                donors.loc[donors["group"].isin(["AD/Dementia","Control"]), "group"])))

if COL["age"]:
    add_row("Age at death (mean±SD, years)",
            fmt_mean_sd(ad[COL["age"]]), fmt_mean_sd(ctrl[COL["age"]]),
            pval_fmt(ttest_p(ad[COL["age"]], ctrl[COL["age"]])))

if COL["pmi"]:
    add_row("PMI (mean±SD, hours)",
            fmt_mean_sd(ad[COL["pmi"]]), fmt_mean_sd(ctrl[COL["pmi"]]),
            pval_fmt(ttest_p(ad[COL["pmi"]], ctrl[COL["pmi"]])))

if COL["braak"]:
    braak_col = donors[COL["braak"]].astype(str).str.strip()
    # Braak distributions
    braak_in_scope = donors[donors["group"].isin(["AD/Dementia","Control"])]
    for stage in ["Braak 0", "Braak I", "Braak II", "Braak III", "Braak IV", "Braak V", "Braak VI"]:
        n_ad   = (ad[COL["braak"]].astype(str).str.strip() == stage).sum()
        n_ctrl = (ctrl[COL["braak"]].astype(str).str.strip() == stage).sum()
        if n_ad + n_ctrl > 0:
            add_row(f"  {stage}",
                    f"{n_ad} ({100*n_ad/max(len(ad),1):.0f}%)",
                    f"{n_ctrl} ({100*n_ctrl/max(len(ctrl),1):.0f}%)")
    add_row("Braak stage (χ² p-value)", "", "",
            pval_fmt(chi2_p(braak_in_scope[COL["braak"]].astype(str).str.strip(),
                             braak_in_scope["group"])))

if COL["apoe"]:
    apoe_col = donors[COL["apoe"]].astype(str).str.strip()
    apoe_vals = sorted([v for v in apoe_col.unique() if v not in ("nan","","None")])
    print(f"\nAPOE values: {apoe_vals}")
    # ε4 carrier = any genotype containing "4"
    e4_ad   = apoe_col[donors["group"] == "AD/Dementia"].str.contains("4", na=False).sum()
    e4_ctrl = apoe_col[donors["group"] == "Control"].str.contains("4", na=False).sum()
    add_row("APOE ε4 carrier (%)",
            f"{e4_ad} ({100*e4_ad/max(len(ad),1):.0f}%)",
            f"{e4_ctrl} ({100*e4_ctrl/max(len(ctrl),1):.0f}%)",
            pval_fmt(chi2_p(apoe_col[donors["group"].isin(["AD/Dementia","Control"])].str.contains("4", na=False),
                             donors.loc[donors["group"].isin(["AD/Dementia","Control"]), "group"])))
    # Full APOE breakdown
    for geno in apoe_vals:
        n_ad   = (apoe_col[donors["group"] == "AD/Dementia"] == geno).sum()
        n_ctrl = (apoe_col[donors["group"] == "Control"] == geno).sum()
        if n_ad + n_ctrl > 0:
            add_row(f"  APOE {geno}",
                    f"{n_ad} ({100*n_ad/max(len(ad),1):.0f}%)",
                    f"{n_ctrl} ({100*n_ctrl/max(len(ctrl),1):.0f}%)")

# ── Output ────────────────────────────────────────────────────────────────
df = pd.DataFrame(rows)
csv_path = OUT / "table1_cohort_demographics.csv"
df.to_csv(csv_path, index=False)
print(f"\nWrote {csv_path}")

# Markdown
md_lines = ["# Table 1 — Cohort Demographics\n"]
# Manual markdown table (avoids tabulate dependency)
cols = list(df.columns)
header = "| " + " | ".join(cols) + " |"
sep    = "| " + " | ".join(["---"] * len(cols)) + " |"
body   = "\n".join("| " + " | ".join(str(df.iloc[i][c]) for c in cols) + " |"
                   for i in range(len(df)))
md_lines.append(header + "\n" + sep + "\n" + body)
md_lines.append("\n\nStatistical tests: continuous variables — two-sample t-test; "
                "categorical variables — Pearson χ² test. "
                "SEA-AD multi-region dataset (Allen Institute for Brain Science). "
                "Cells: microglial nuclei from snRNA-seq after quality filtering.")
md_path = OUT / "table1_cohort_demographics.md"
md_path.write_text("\n".join(md_lines))
print(f"Wrote {md_path}")

print("\n--- Table preview ---")
print(df.to_string(index=False))
