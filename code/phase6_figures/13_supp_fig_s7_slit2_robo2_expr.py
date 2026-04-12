"""
Supplementary Figure S8 — SLIT2 and ROBO2 expression across MTG cell types
============================================================================
Dot plot showing:
  - SLIT2 expression across all 15 CellChat cell types (sender perspective)
  - ROBO2 expression across all 15 CellChat cell types (receiver perspective)

Data: counts_raw.h5 + gene_names.csv + cell_meta.csv from CellChat prep
Outputs:
  results/figures/supp_fig_S7_slit2_robo2_expr.pdf
  results/figures/supp_fig_S7_slit2_robo2_expr.png
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

PROJECT  = Path(__file__).resolve().parents[2]
LR_PREP  = PROJECT / "results/phase2/LR/prep"
OUT_DIR  = PROJECT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GENES_CSV = LR_PREP / "gene_names.csv"
CELL_META = LR_PREP / "cell_meta.csv"
COUNTS_H5 = LR_PREP / "counts_raw.h5"

GENES_OF_INTEREST = ["SLIT2", "ROBO2", "GAS6", "MERTK"]

INTERNEURON_TYPES = ["Chandelier", "Lamp5", "Lamp5 Lhx6", "Pvalb", "Sncg",
                     "Sst", "Sst Chodl", "Vip"]


def main():
    # ── Load gene names ────────────────────────────────────────────────────
    print("Loading gene names...")
    gene_names = pd.read_csv(GENES_CSV)["gene"].tolist()
    print(f"  {len(gene_names):,} genes")

    target_idx = {}
    for g in GENES_OF_INTEREST:
        if g in gene_names:
            target_idx[g] = gene_names.index(g)
        else:
            print(f"  WARNING: {g} not found in gene list")

    print(f"  Gene indices: {target_idx}")

    # ── Load cell metadata ─────────────────────────────────────────────────
    print("Loading cell metadata...")
    meta = pd.read_csv(CELL_META)
    print(f"  {len(meta):,} cells")
    print(f"  Cell types: {meta['cell_type'].value_counts().to_dict()}")

    # ── Load count matrix (sparse H5) ─────────────────────────────────────
    print("Loading count matrix from H5...")
    with h5py.File(COUNTS_H5, "r") as f:
        print(f"  H5 keys: {list(f.keys())}")
        # Flat CSC format: data, indices, indptr at top level
        data    = f["data"][:]
        indices = f["indices"][:]
        indptr  = f["indptr"][:]
        n_genes = len(gene_names)
        n_cells = len(meta)
        # indptr length = n_cols + 1; determine orientation
        if len(indptr) == n_genes + 1:
            # genes × cells (CSC with genes as columns)
            X = sp.csc_matrix((data, indices, indptr), shape=(n_cells, n_genes))
        elif len(indptr) == n_cells + 1:
            # cells × genes stored as CSC (cells as columns) — transpose
            X = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T
        else:
            raise ValueError(
                f"Cannot infer shape: indptr len={len(indptr)}, "
                f"n_cells={n_cells}, n_genes={n_genes}"
            )

    print(f"  Matrix shape: {X.shape} (cells × genes)")

    # Transpose if needed: should be (n_cells, n_genes)
    if X.shape[0] == len(gene_names) and X.shape[1] == len(meta):
        X = X.T
        print(f"  Transposed to: {X.shape}")

    # ── Extract target gene expression ─────────────────────────────────────
    print("Extracting target gene expression...")
    expr = {}
    for g, idx in target_idx.items():
        col   = np.asarray(X[:, idx].todense()).ravel()
        total = np.asarray(X.sum(axis=1)).ravel()
        expr[g] = np.log1p(col / np.clip(total, 1, None) * 1e6)

    expr_df = pd.DataFrame(expr, index=meta.index)
    expr_df["cell_type"] = meta["cell_type"].values

    # ── Compute per-cell-type stats ────────────────────────────────────────
    print("Computing per-cell-type statistics...")
    cell_types = meta["cell_type"].value_counts().index.tolist()

    rows = []
    for ct in cell_types:
        mask = expr_df["cell_type"] == ct
        for g in GENES_OF_INTEREST:
            if g not in expr_df.columns:
                continue
            vals = expr_df.loc[mask, g].values
            rows.append({
                "cell_type": ct,
                "gene":      g,
                "mean_expr": vals.mean(),
                "pct_expr":  (vals > 0).mean() * 100,
                "n_cells":   mask.sum(),
            })

    stats = pd.DataFrame(rows)
    print(stats.groupby("gene")[["mean_expr", "pct_expr"]].describe())

    # ── Plot: dot plot ─────────────────────────────────────────────────────
    genes_plot = [g for g in ["SLIT2", "ROBO2", "GAS6", "MERTK"] if g in target_idx]

    # Order cell types: interneurons first (SLIT2 senders), microglia last
    CT_ORDER = (
        [ct for ct in INTERNEURON_TYPES if ct in cell_types] +
        [ct for ct in cell_types if ct not in INTERNEURON_TYPES]
    )
    CT_ORDER = [ct for ct in CT_ORDER if ct in cell_types]

    stats_pivot_mean = stats.pivot(
        index="cell_type", columns="gene", values="mean_expr"
    ).reindex(CT_ORDER)
    stats_pivot_pct  = stats.pivot(
        index="cell_type", columns="gene", values="pct_expr"
    ).reindex(CT_ORDER)

    fig, axes = plt.subplots(1, len(genes_plot), figsize=(3.5 * len(genes_plot), 7))
    if len(genes_plot) == 1:
        axes = [axes]

    for ax, gene in zip(axes, genes_plot):
        mean_vals = stats_pivot_mean[gene].fillna(0).values
        pct_vals  = stats_pivot_pct[gene].fillna(0).values

        vmax  = np.percentile(mean_vals[mean_vals > 0], 95) if (mean_vals > 0).any() else 1
        norm  = mcolors.Normalize(vmin=0, vmax=vmax)
        cmap  = plt.cm.YlOrRd
        x     = np.zeros(len(CT_ORDER))
        y     = np.arange(len(CT_ORDER))
        sizes = (pct_vals / 100) ** 0.5 * 350  # sqrt scaling

        sc = ax.scatter(x, y, s=sizes, c=mean_vals, cmap=cmap, vmin=0, vmax=vmax,
                        edgecolors="grey", linewidths=0.3)

        ax.set_yticks(y)
        if gene == genes_plot[0]:
            ax.set_yticklabels(CT_ORDER, fontsize=8.5)
            for i, ct in enumerate(CT_ORDER):
                if ct in INTERNEURON_TYPES:
                    ax.get_yticklabels()[i].set_color("#1565C0")
                elif ct == "Microglia":
                    ax.get_yticklabels()[i].set_color("#B71C1C")
        else:
            ax.set_yticklabels([])

        ax.set_xticks([])
        ax.set_title(f"*{gene}*\n(log-norm CPM)", fontsize=10, style="italic")
        ax.invert_yaxis()
        ax.set_xlim(-0.5, 0.5)
        ax.spines[["top", "right", "bottom"]].set_visible(False)

        cb = plt.colorbar(sc, ax=ax, shrink=0.4, pad=0.08)
        cb.set_label("Mean\nlog-norm expr", fontsize=7)
        cb.ax.tick_params(labelsize=7)

    # Legend for dot size
    legend_pct     = [10, 25, 50]
    legend_handles = [
        plt.scatter([], [], s=(p/100)**0.5 * 350, color="grey", edgecolors="grey",
                    label=f"{p}%") for p in legend_pct
    ]
    fig.legend(handles=legend_handles, title="% cells\nexpressing",
               loc="lower right", fontsize=7, title_fontsize=8,
               frameon=True, bbox_to_anchor=(0.98, 0.02))

    fig.text(0.5, 0.97,
             "SLIT2/ROBO2 expression across MTG cell types (SEA-AD; n=45,000 cells)\n"
             "Blue labels = inhibitory interneurons (SLIT2 senders); Red = Microglia (ROBO2 receiver)",
             ha="center", va="top", fontsize=9, style="italic")

    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(OUT_DIR / "supp_fig_S7_slit2_robo2_expr.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "supp_fig_S7_slit2_robo2_expr.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: supp_fig_S7_slit2_robo2_expr.pdf/.png")


if __name__ == "__main__":
    main()
