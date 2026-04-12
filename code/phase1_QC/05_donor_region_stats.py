"""
DAM-DRUG — Donor/Region Summary Statistics
===========================================
Addresses ChatGPT reviewer Comment 5:
  - Substate cell counts by donor and by region
  - Donor-level pseudotime summaries
  - IKZF1 AUCell by donor (if available in obs)
  - Region composition per microglial substate

Outputs:
  results/phase1/donor_region/substate_by_donor.csv
  results/phase1/donor_region/substate_by_region.csv
  results/phase1/donor_region/pseudotime_by_donor.csv
  results/phase1/donor_region/ikzf1_auc_by_donor.csv   (if AUC in obs)
  results/phase1/donor_region/donor_summary.csv
"""

import os
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
H5AD    = PROJECT / "results/phase1/trajectory/microglia_trajectory.h5ad"
OUT     = PROJECT / "results/phase1/donor_region"
OUT.mkdir(parents=True, exist_ok=True)

print("Loading obs only from microglia_trajectory.h5ad ...")
adata = ad.read_h5ad(H5AD, backed="r")
obs   = adata.obs.copy()
adata.file.close()
print(f"  {obs.shape[0]:,} cells × obs columns: {list(obs.columns)}")

# ── Identify key columns ────────────────────────────────────────────────────
# Try common naming patterns for donor and region
DONOR_COL  = next((c for c in obs.columns if c in ("Donor ID", "donor_id", "donor")), None)
REGION_COL = next((c for c in obs.columns if c.lower() in ("region", "brain_region",
                                                             "region_of_interest",
                                                             "roi", "anatomical_division_label")), None)
STATE_COL  = next((c for c in obs.columns if c in ("microglial_state", "state",
                                                    "Supertype_mapped", "substate")), None)

# pySCENIC AUCell scores may be stored as obs columns like "IKZF1(+)"
AUC_COL = next((c for c in obs.columns if "IKZF1" in c), None)

# dpt_pseudotime
PT_COL = next((c for c in obs.columns if "pseudotime" in c.lower()), None)

print(f"  Donor col:  {DONOR_COL}")
print(f"  Region col: {REGION_COL}")
print(f"  State col:  {STATE_COL}")
print(f"  AUC col:    {AUC_COL}")
print(f"  PT col:     {PT_COL}")

# Fallback: if Supertype present, remap to state names
SUPERTYPE_TO_STATE = {
    "Micro-PVM_1":         "Homeostatic",
    "Micro-PVM_2":         "DAM",
    "Micro-PVM_2_3-SEAAD": "DAM-IRM",
    "Micro-PVM_2_1-SEAAD": "LateAD-DAM",
    "Micro-PVM_3-SEAAD":   "IRM",
    "Micro-PVM_4-SEAAD":   "LAM",
}

if STATE_COL is None and "Supertype" in obs.columns:
    obs["microglial_state"] = obs["Supertype"].map(SUPERTYPE_TO_STATE)
    STATE_COL = "microglial_state"

if STATE_COL is None:
    # Try to find any column with known state names
    for c in obs.columns:
        if obs[c].astype(str).isin(SUPERTYPE_TO_STATE.values()).any():
            STATE_COL = c
            break

print(f"  Final state col: {STATE_COL}")
if STATE_COL:
    print(f"  State counts:\n{obs[STATE_COL].value_counts()}")

# ── 1. Substate × donor cell counts ────────────────────────────────────────
if DONOR_COL and STATE_COL:
    ct_donor = (obs.groupby([STATE_COL, DONOR_COL])
                   .size()
                   .rename("n_cells")
                   .reset_index())
    ct_donor_pivot = ct_donor.pivot(index=DONOR_COL, columns=STATE_COL, values="n_cells").fillna(0).astype(int)
    ct_donor_pivot["Total"] = ct_donor_pivot.sum(axis=1)
    ct_donor_pivot.to_csv(OUT / "substate_by_donor.csv")
    print(f"\nSubstate × donor: {ct_donor_pivot.shape}")
    print(ct_donor_pivot.describe())

    # Summary: n donors contributing each state
    donor_per_state = ct_donor.groupby(STATE_COL).agg(
        n_donors=(DONOR_COL, "nunique"),
        mean_cells=("n_cells", "mean"),
        median_cells=("n_cells", "median"),
        min_cells=("n_cells", "min"),
        max_cells=("n_cells", "max"),
    )
    donor_per_state.to_csv(OUT / "donor_summary.csv")
    print(f"\nDonor summary per state:\n{donor_per_state.to_string()}")

# ── 2. Substate × region cell counts ───────────────────────────────────────
if REGION_COL and STATE_COL:
    ct_region = (obs.groupby([STATE_COL, REGION_COL])
                    .size()
                    .rename("n_cells")
                    .reset_index())
    ct_region_pivot = ct_region.pivot(index=REGION_COL, columns=STATE_COL, values="n_cells").fillna(0).astype(int)
    ct_region_pivot["Total"] = ct_region_pivot.sum(axis=1)
    ct_region_pivot.to_csv(OUT / "substate_by_region.csv")
    print(f"\nSubstate × region:\n{ct_region_pivot.to_string()}")
else:
    print(f"\n[WARN] Region column not found — skipping region breakdown.")
    print(f"  Available columns: {[c for c in obs.columns if 'region' in c.lower() or 'area' in c.lower()]}")

# ── 3. Pseudotime by donor ──────────────────────────────────────────────────
if PT_COL and DONOR_COL and STATE_COL:
    pt_donor = (obs.groupby([DONOR_COL, STATE_COL])[PT_COL]
                   .agg(["median", "mean", "std", "count"])
                   .round(4)
                   .reset_index())
    pt_donor.to_csv(OUT / "pseudotime_by_donor.csv", index=False)
    print(f"\nPseudotime by donor (first 10 rows):\n{pt_donor.head(10).to_string()}")

    # Per-donor overall median (test if IKZF1 pattern holds across donors)
    pt_state = (obs.groupby(STATE_COL)[PT_COL]
                   .agg(["median", "mean", "std", "count"])
                   .round(4))
    print(f"\nMedian pseudotime per state:\n{pt_state.to_string()}")
else:
    print(f"\n[WARN] Pseudotime column not found.")
    print(f"  PT_COL={PT_COL}, DONOR_COL={DONOR_COL}, STATE_COL={STATE_COL}")

# ── 4. IKZF1 AUCell by donor ───────────────────────────────────────────────
if AUC_COL and DONOR_COL and STATE_COL:
    auc_donor = (obs.groupby([DONOR_COL, STATE_COL])[AUC_COL]
                    .agg(["median", "mean", "std", "count"])
                    .round(5)
                    .reset_index())
    auc_donor.to_csv(OUT / "ikzf1_auc_by_donor.csv", index=False)

    # Donor-level: is IKZF1 highest in LateAD-DAM for each donor?
    auc_pivot = auc_donor.pivot(index=DONOR_COL, columns=STATE_COL, values="median")
    if "LateAD-DAM" in auc_pivot.columns:
        lateAD_highest = (auc_pivot.idxmax(axis=1) == "LateAD-DAM").sum()
        total_donors = auc_pivot.dropna(subset=["LateAD-DAM"]).shape[0]
        print(f"\nIKZF1 AUCell highest in LateAD-DAM: {lateAD_highest}/{total_donors} donors")
    print(f"\nIKZF1 AUCell by donor (first 10):\n{auc_donor.head(10).to_string()}")
else:
    print(f"\n[WARN] IKZF1 AUCell column not found in obs. AUC_COL={AUC_COL}")
    print(f"  Columns containing 'AUC' or 'auc' or 'IKZF1': "
          f"{[c for c in obs.columns if 'auc' in c.lower() or 'IKZF1' in c]}")

# ── 5. Disease group breakdown ──────────────────────────────────────────────
DX_COL = next((c for c in obs.columns if "diagnosis" in c.lower() or
               c in ("Overall AD neuropathological Change", "Cognitive Status")), None)
if DX_COL and DONOR_COL:
    dx_donors = obs[[DONOR_COL, DX_COL]].drop_duplicates()
    print(f"\nDisease group counts:\n{dx_donors[DX_COL].value_counts().to_string()}")

BRAAK_COL = next((c for c in obs.columns if "braak" in c.lower()), None)
if BRAAK_COL and DONOR_COL:
    braak_donors = obs[[DONOR_COL, BRAAK_COL]].drop_duplicates()
    print(f"\nBraak distribution:\n{braak_donors[BRAAK_COL].value_counts().sort_index().to_string()}")

print("\nDone. Outputs written to:", OUT)
