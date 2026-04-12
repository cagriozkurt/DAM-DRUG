"""
Supplementary Figure S6 — BHLHE40/41 Motif-Agnostic Co-expression Rescue
=========================================================================
BHLHE40 and BHLHE41 lack HOCOMOCO v11 motifs and cannot be evaluated by
pySCENIC. This script computes pseudobulk Spearman co-expression across
84 donor samples as a surrogate regulon and tests enrichment of DAM/LAM
marker genes in the top co-expressed gene set.

Outputs:
  results/figures/supp_fig_S6_bhlhe_coexpr.pdf
  results/figures/supp_fig_S6_bhlhe_coexpr.png
  results/phase2/GRN/bhlhe40_coexpr_top200.csv
  results/phase2/GRN/bhlhe41_coexpr_top200.csv
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr, fisher_exact
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT  = Path(__file__).resolve().parents[2]
RAW_DIR  = PROJECT / "data/raw/SEA-AD"
OUT_DIR  = PROJECT / "results/figures"
GRN_DIR  = PROJECT / "results/phase2/GRN"
OUT_DIR.mkdir(parents=True, exist_ok=True)
GRN_DIR.mkdir(parents=True, exist_ok=True)

# ── DAM/LAM marker gene sets (from Table 2 / pseudobulk DGE) ──────────────────
DAM_MARKERS = {
    "TREM2", "TYROBP", "APOE", "CST7", "LPL", "CLEC7A", "CD68",
    "ADAMDEC1", "ITGAX", "SPP1", "CCL3", "CCL4"
}
LAM_MARKERS = {
    "LGALS3", "GPNMB", "CD9", "FABP5", "LIPA", "NPC2",
    "CCL3L1", "CCL4L2", "RGCC", "PLIN2"
}
DISEASE_MARKERS = DAM_MARKERS | LAM_MARKERS

TOP_N = 200  # surrogate regulon size


def find_seaad_h5ad(raw_dir: Path) -> Path:
    """Locate the SEA-AD microglia h5ad, tolerating filename date changes."""
    candidates = sorted(raw_dir.glob("SEA-AD_Microglia_multi-regional_final-nuclei.*.h5ad"))
    if not candidates:
        raise FileNotFoundError(
            f"No SEA-AD microglia h5ad found in {raw_dir}. "
            "Expected: SEA-AD_Microglia_multi-regional_final-nuclei.*.h5ad"
        )
    if len(candidates) > 1:
        print(f"  WARNING: multiple SEA-AD h5ad files found; using most recent: {candidates[-1].name}")
    return candidates[-1]


def main():
    data_path = find_seaad_h5ad(RAW_DIR)

    # ── Load data and build pseudobulk ────────────────────────────────────
    print(f"Loading AnnData (backed mode): {data_path.name}")
    adata = sc.read_h5ad(data_path, backed="r")

    # Build pseudobulk: sum raw counts per donor
    print("Building pseudobulk matrix...")
    donors = adata.obs["Donor ID"].unique()
    genes  = adata.var_names.tolist()

    pb_rows   = []
    donor_ids = []
    for donor in donors:
        mask   = adata.obs["Donor ID"] == donor
        counts = np.asarray(adata[mask].X.sum(axis=0)).flatten()
        pb_rows.append(counts)
        donor_ids.append(donor)

    pb = pd.DataFrame(pb_rows, index=donor_ids, columns=genes)

    # Log-normalise (library-size normalise then log1p)
    pb_norm = pb.div(pb.sum(axis=1), axis=0) * 1e6
    pb_norm = np.log1p(pb_norm)

    # Filter: keep genes expressed in ≥ 10 donors
    expressed = (pb > 0).sum(axis=0) >= 10
    pb_norm   = pb_norm.loc[:, expressed]
    print(f"Expressed genes after filtering: {expressed.sum()}")

    all_genes = pb_norm.columns.tolist()

    # ── Compute Spearman co-expression for BHLHE40 and BHLHE41 ───────────
    results = {}
    for tf in ["BHLHE40", "BHLHE41"]:
        if tf not in pb_norm.columns:
            print(f"WARNING: {tf} not found in pseudobulk matrix — skipping")
            continue
        print(f"Computing Spearman correlations for {tf}...")
        tf_expr = pb_norm[tf].values
        corrs   = {}
        for gene in all_genes:
            if gene == tf:
                continue
            r, _ = spearmanr(tf_expr, pb_norm[gene].values)
            corrs[gene] = r
        corr_series = pd.Series(corrs).sort_values(ascending=False)
        results[tf] = corr_series

        # Save top 200
        top200 = corr_series.head(TOP_N).reset_index()
        top200.columns = ["gene", "spearman_r"]
        top200.to_csv(GRN_DIR / f"{tf.lower()}_coexpr_top200.csv", index=False)
        print(f"  Saved top {TOP_N} co-expressed genes for {tf}")

    # ── Fisher's exact test: enrichment of disease markers in top 200 ────
    print("\nFisher's exact test for DAM/LAM marker enrichment:")
    for tf, corr_series in results.items():
        top200_genes = set(corr_series.head(TOP_N).index)
        background   = set(all_genes) - {tf}
        a = len(top200_genes & DISEASE_MARKERS)
        b = len(top200_genes - DISEASE_MARKERS)
        c = len((background - top200_genes) & DISEASE_MARKERS)
        d = len((background - top200_genes) - DISEASE_MARKERS)
        or_, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        print(f"  {tf}: {a}/{TOP_N} disease markers in top200 | OR={or_:.2f} | p={p:.3e}")

    # ── Plot ──────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    for idx, tf in enumerate(["BHLHE40", "BHLHE41"]):
        if tf not in results:
            continue
        corr_series  = results[tf]
        top200_genes = corr_series.head(TOP_N)

        # Panel A/C: top 30 co-expressed genes bar plot
        ax_bar = fig.add_subplot(gs[idx, 0])
        top30  = top200_genes.head(30)
        colors = ["#E64B35" if g in DISEASE_MARKERS else "#4DBBD5" for g in top30.index]
        ax_bar.barh(range(len(top30)), top30.values[::-1], color=colors[::-1])
        ax_bar.set_yticks(range(len(top30)))
        ax_bar.set_yticklabels(top30.index[::-1], fontsize=7)
        ax_bar.set_xlabel("Spearman ρ", fontsize=9)
        ax_bar.set_title(f"{tf} — Top 30 co-expressed genes\n(red = DAM/LAM marker)", fontsize=9)
        ax_bar.axvline(0, color="black", linewidth=0.5)

        # Panel B/D: genome-wide correlation distribution
        ax_dist = fig.add_subplot(gs[idx, 1])
        ax_dist.hist(corr_series.values, bins=80, color="#8491B4", alpha=0.8, edgecolor="none")
        ax_dist.axvline(top200_genes.iloc[-1], color="#E64B35", linestyle="--", linewidth=1,
                        label=f"Top {TOP_N} cutoff (ρ={top200_genes.iloc[-1]:.3f})")
        ax_dist.set_xlabel("Spearman ρ", fontsize=9)
        ax_dist.set_ylabel("Gene count", fontsize=9)
        ax_dist.set_title(f"{tf} — Genome-wide co-expression distribution", fontsize=9)
        ax_dist.legend(fontsize=8)

    fig.suptitle(
        "Supplementary Figure S6 — BHLHE40/41 Motif-Agnostic Co-expression Rescue\n"
        "Pseudobulk Spearman correlations across 84 donors (surrogate for pySCENIC regulon)",
        fontsize=10, y=1.01
    )

    for ext in ("pdf", "png"):
        out_path = OUT_DIR / f"supp_fig_S6_bhlhe_coexpr.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved {out_path}")

    plt.close()
    print("Done.")


if __name__ == "__main__":
    main()
