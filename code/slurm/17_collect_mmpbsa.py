#!/usr/bin/env python3
"""
17_collect_mmpbsa.py
--------------------
Aggregate mmpbsa_result.txt files from all 26 MM-GBSA workdirs into a
ranked CSV.  Run from the project root or pass PROJDIR as first argument.

Output: results/phase4/mmpbsa/mmpbsa_summary.csv
Columns: target, chembl_id, vina_score, dg_gbsa, sem, rank_vina, rank_gbsa
"""

import sys
import csv
from pathlib import Path

PROJDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/arf/scratch/mozkurt/DAM-DRUG"
)
MMPBSA_ROOT = PROJDIR / "results" / "phase4" / "mmpbsa"
OUT_CSV = MMPBSA_ROOT / "mmpbsa_summary.csv"


def main():
    rows = []
    missing = []

    for result_file in sorted(MMPBSA_ROOT.rglob("mmpbsa_result.txt")):
        txt = result_file.read_text().strip()
        if not txt:
            missing.append(result_file)
            continue
        parts = txt.split(",")
        if len(parts) != 5:
            print(f"WARN: unexpected format in {result_file}: {txt!r}")
            continue
        target, chembl, vina, dg, sem = parts
        rows.append(
            {
                "target": target,
                "chembl_id": chembl,
                "vina_score": float(vina),
                "dg_gbsa": float(dg),
                "sem": float(sem),
            }
        )

    if missing:
        print(f"WARNING: {len(missing)} workdir(s) have empty/missing result files:")
        for p in missing:
            print(f"  {p}")

    if not rows:
        print("ERROR: No results found — check that mmpbsa_result.txt files exist.")
        sys.exit(1)

    # Rank within each target by vina (ascending = more negative = better)
    # and by dg_gbsa (ascending)
    from collections import defaultdict

    by_target = defaultdict(list)
    for r in rows:
        by_target[r["target"]].append(r)

    ranked = []
    for target, group in sorted(by_target.items()):
        for rank_v, r in enumerate(
            sorted(group, key=lambda x: x["vina_score"]), start=1
        ):
            r["rank_vina"] = rank_v
        for rank_g, r in enumerate(
            sorted(group, key=lambda x: x["dg_gbsa"]), start=1
        ):
            r["rank_gbsa"] = rank_g
        ranked.extend(group)

    # Global sort: by dg_gbsa ascending (most negative first)
    ranked.sort(key=lambda x: x["dg_gbsa"])

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "target", "chembl_id", "vina_score",
                "dg_gbsa", "sem", "rank_vina", "rank_gbsa",
            ],
        )
        writer.writeheader()
        writer.writerows(ranked)

    print(f"Wrote {len(ranked)} results to {OUT_CSV}")
    print()
    print(f"{'Target':<30} {'CHEMBL':<15} {'Vina':>8} {'ΔG_GBSA':>10} {'SEM':>6}  rank_v rank_g")
    print("-" * 85)
    for r in ranked:
        print(
            f"{r['target']:<30} {r['chembl_id']:<15} {r['vina_score']:>8.3f} "
            f"{r['dg_gbsa']:>10.2f} {r['sem']:>6.2f}  "
            f"{r['rank_vina']:>6} {r['rank_gbsa']:>6}"
        )


if __name__ == "__main__":
    main()
