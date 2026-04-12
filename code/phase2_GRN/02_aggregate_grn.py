"""
DAM-DRUG Phase 2 — Aggregate GRNBoost2 adjacency matrices across seeds
=======================================================================
Reads adj_matrix.tsv (seed=42) + adj_matrix_seed{1..4}.tsv and computes
mean edge weight across all 5 seeds.

Edges present in fewer than min_seeds seeds are discarded (default: 3/5)
to remove seed-specific noise.

Output:
  results/phase2/GRN/adj_matrix_aggregated.tsv
    columns: TF, target, importance_mean, importance_std, n_seeds
    sorted by importance_mean descending

Run (after all 4 array jobs complete):
  conda run -n scenic python code/phase2_GRN/25_aggregate_grn.py
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG"))
GRN_DIR = PROJECT / "results/phase2/GRN"

# All 5 seed files
SEED_FILES = {
    42: GRN_DIR / "adj_matrix.tsv",
    1:  GRN_DIR / "adj_matrix_seed1.tsv",
    2:  GRN_DIR / "adj_matrix_seed2.tsv",
    3:  GRN_DIR / "adj_matrix_seed3.tsv",
    4:  GRN_DIR / "adj_matrix_seed4.tsv",
}
OUT = GRN_DIR / "adj_matrix_aggregated.tsv"
MIN_SEEDS = 3   # edge must appear in at least 3/5 seeds

# ── Load all available adj matrices ──────────────────────────────────────────
frames = []
for seed, path in SEED_FILES.items():
    if not path.exists():
        print(f"WARNING: seed={seed} not found ({path.name}) — skipping")
        continue
    df = pd.read_csv(path, sep="\t")
    # pySCENIC GRN output columns: TF, target, importance
    df = df[["TF", "target", "importance"]].copy()
    df["seed"] = seed
    frames.append(df)
    print(f"  seed={seed:2d}: {len(df):,} edges loaded from {path.name}")

if not frames:
    print("ERROR: no adjacency files found")
    sys.exit(1)

n_seeds = len(frames)
print(f"\n{n_seeds} seed file(s) loaded")

# ── Aggregate ─────────────────────────────────────────────────────────────────
all_edges = pd.concat(frames, ignore_index=True)

agg = (
    all_edges
    .groupby(["TF", "target"])["importance"]
    .agg(
        importance_mean="mean",
        importance_std="std",
        n_seeds="count",
    )
    .reset_index()
)

before = len(agg)
agg = agg[agg["n_seeds"] >= MIN_SEEDS].copy()
print(f"Edges before filtering: {before:,}")
print(f"Edges with n_seeds >= {MIN_SEEDS}: {len(agg):,}  "
      f"({100*len(agg)/before:.1f}% retained)")

agg = agg.sort_values("importance_mean", ascending=False).reset_index(drop=True)

# ── Save ──────────────────────────────────────────────────────────────────────
agg.to_csv(OUT, sep="\t", index=False)
print(f"\nAggregated adjacency matrix → {OUT}")
print(f"  {len(agg):,} edges  ({agg['TF'].nunique()} TFs, {agg['target'].nunique()} targets)")

# ── TF target summary ─────────────────────────────────────────────────────────
TF_TARGETS = ["SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
              "RELB", "BHLHE40", "BHLHE41"]
print("\nTop aggregated edges for TF targets:")
print(f"{'TF':<12} {'target':<20} {'mean_imp':>9} {'std':>7} {'n_seeds':>8}")
print("-" * 62)
for tf in TF_TARGETS:
    top = agg[agg["TF"] == tf].head(3)
    if top.empty:
        print(f"{tf:<12}  (no edges in aggregated matrix)")
    else:
        for _, row in top.iterrows():
            print(f"{row.TF:<12} {row.target:<20} {row.importance_mean:9.4f} "
                  f"{(row.importance_std or 0):7.4f} {int(row.n_seeds):8d}")

print(f"\nNext: re-run pyscenic ctx using {OUT.name} as input adjacency matrix")
