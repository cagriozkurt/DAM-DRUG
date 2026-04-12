"""
DAM-DRUG Phase 2 — Regulon–Pseudotime Spearman Correlation (item 2.7)
======================================================================
Correlates per-cell AUCell regulon activity with diffusion pseudotime.

Inputs:
  results/phase2/GRN/scenic_auc_aggregated.loom   (100K cells × 46 regulons)
  results/phase1/trajectory/pseudotime.csv         (236K cells)

Outputs:
  results/phase2/GRN/regulon_pseudotime_corr.csv
    columns: regulon, TF, rho, pval, padj, n_cells, peak_state, direction
  results/phase2/GRN/regulon_pseudotime_plot.pdf
    Smoothed AUC vs pseudotime for significant regulons (|rho| > 0.3)
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.ndimage import uniform_filter1d
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import loompy

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
GRN      = PROJECT / "results/phase2/GRN"
TRAJ     = PROJECT / "results/phase1/trajectory"

AUC_LOOM = GRN / "scenic_auc_aggregated.loom"
PT_CSV   = TRAJ / "pseudotime.csv"
OUT_CORR = GRN / "regulon_pseudotime_corr.csv"
OUT_PLOT = GRN / "regulon_pseudotime_plot.pdf"

RHO_THRESH = 0.3     # effect-size filter (100K cells → any ρ is p<0.001)
N_BINS     = 50      # pseudotime bins for smoothed trajectory plot

TF_TARGETS = {"SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
              "RELB", "BHLHE40", "BHLHE41"}


def main():
    # ── Load AUC loom ──────────────────────────────────────────────────────────
    log.info(f"Loading AUC loom: {AUC_LOOM.name}")
    with loompy.connect(str(AUC_LOOM), mode="r", validate=False) as ds:
        auc_raw       = ds.ca["RegulonsAUC"]
        regulon_names = list(auc_raw.dtype.names)
        cell_ids      = list(ds.ca["CellID"])
        states        = list(ds.ca["state"])
        auc_df = pd.DataFrame(
            {r: auc_raw[r] for r in regulon_names},
            index=cell_ids,
        )
        auc_df["state"] = states
    log.info(f"  AUC matrix: {auc_df.shape[0]:,} cells × {len(regulon_names)} regulons")

    # ── Load pseudotime ────────────────────────────────────────────────────────
    log.info(f"Loading pseudotime: {PT_CSV.name}")
    pt = pd.read_csv(PT_CSV, index_col="cell_id")
    log.info(f"  Pseudotime: {len(pt):,} cells")

    # ── Inner join ─────────────────────────────────────────────────────────────
    merged = auc_df.join(pt[["dpt_pseudotime"]], how="inner")
    n_cells = len(merged)
    log.info(f"  After inner join: {n_cells:,} cells")
    if n_cells < 1000:
        log.error("Fewer than 1000 cells in overlap — check cell ID format mismatch")
        sys.exit(1)

    pseudotime = merged["dpt_pseudotime"].values

    # ── Spearman correlations ──────────────────────────────────────────────────
    log.info("Computing Spearman correlations...")
    rows = []
    for reg in regulon_names:
        auc_vals = merged[reg].values
        rho, pval = stats.spearmanr(auc_vals, pseudotime)
        tf = reg.split("(")[0]
        rows.append({"regulon": reg, "TF": tf, "rho": rho, "pval": pval})

    corr_df = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(corr_df["pval"], method="fdr_bh")
    corr_df["padj"]   = padj
    corr_df["n_cells"] = n_cells
    corr_df["direction"] = np.where(corr_df["rho"] > 0, "disease_up", "homeostatic_up")

    # Add peak state from per-state AUC means
    state_mean = merged.groupby("state")[regulon_names].mean()
    corr_df["peak_state"] = corr_df["regulon"].map(
        lambda r: state_mean[r].idxmax() if r in state_mean.columns else "NA"
    )

    corr_df = corr_df.sort_values("rho", ascending=False).reset_index(drop=True)
    corr_df.to_csv(OUT_CORR, index=False)
    log.info(f"  Saved: {OUT_CORR}")

    # ── Report ─────────────────────────────────────────────────────────────────
    sig = corr_df[corr_df["rho"].abs() >= RHO_THRESH]
    log.info(f"\n{'='*65}")
    log.info(f"REGULON–PSEUDOTIME CORRELATION  (|rho| >= {RHO_THRESH}, n={n_cells:,} cells)")
    log.info(f"{'='*65}")
    log.info(f"{'Regulon':<35} {'rho':>6} {'padj':>9} {'peak_state':<14} {'target?'}")
    log.info("-" * 65)
    for _, row in corr_df.iterrows():
        flag  = "*** TARGET" if row["TF"] in TF_TARGETS else ""
        mark  = ">" if abs(row["rho"]) >= RHO_THRESH else " "
        log.info(f"{mark} {row['regulon']:<33} {row['rho']:>6.3f} {row['padj']:>9.2e} "
                 f"{row['peak_state']:<14} {flag}")
    log.info(f"\n{len(sig)} regulons with |rho| >= {RHO_THRESH}")

    # Summary for TF targets
    log.info(f"\nTF TARGET SUMMARY:")
    for tf in sorted(TF_TARGETS):
        matches = corr_df[corr_df["TF"] == tf]
        if matches.empty:
            log.info(f"  {tf:<12}  no regulon")
        else:
            for _, r in matches.iterrows():
                log.info(f"  {r['regulon']:<35} rho={r['rho']:+.3f}  peak={r['peak_state']}")

    # ── Plot: smoothed AUC vs pseudotime ───────────────────────────────────────
    if sig.empty:
        log.info("No significant regulons to plot.")
        return

    sig_regs = sig["regulon"].tolist()
    n_regs   = len(sig_regs)

    # Sort cells by pseudotime and smooth in bins
    merged_sorted = merged.sort_values("dpt_pseudotime")
    pt_sorted     = merged_sorted["dpt_pseudotime"].values

    ncols = min(4, n_regs)
    nrows = int(np.ceil(n_regs / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows),
                             squeeze=False)

    for idx, reg in enumerate(sig_regs):
        ax  = axes[idx // ncols][idx % ncols]
        auc = merged_sorted[reg].values
        # Bin-smooth
        bin_edges = np.linspace(pt_sorted[0], pt_sorted[-1], N_BINS + 1)
        bin_mid, bin_auc = [], []
        for i in range(N_BINS):
            mask = (pt_sorted >= bin_edges[i]) & (pt_sorted < bin_edges[i + 1])
            if mask.sum() > 0:
                bin_mid.append((bin_edges[i] + bin_edges[i + 1]) / 2)
                bin_auc.append(auc[mask].mean())

        rho = corr_df.loc[corr_df["regulon"] == reg, "rho"].values[0]
        col = "#CC3333" if rho > 0 else "#3366CC"
        ax.plot(bin_mid, bin_auc, color=col, linewidth=2)
        ax.set_title(f"{reg}\nρ={rho:+.3f}", fontsize=9)
        ax.set_xlabel("Pseudotime", fontsize=8)
        ax.set_ylabel("AUC", fontsize=8)
        ax.tick_params(labelsize=7)

    # Hide empty panels
    for idx in range(n_regs, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    plt.suptitle(f"Regulon activity vs diffusion pseudotime  (|ρ| ≥ {RHO_THRESH})",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    fig.savefig(OUT_PLOT, bbox_inches="tight")
    plt.close()
    log.info(f"  Plot saved: {OUT_PLOT}")
    log.info("\nDone.")


if __name__ == "__main__":
    main()
