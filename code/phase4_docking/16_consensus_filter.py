"""
DAM-DRUG Phase 4 — Consensus filter & candidate shortlist
==========================================================
Combines Vina + Gnina CNN scores across all targets, applies
target-specific thresholds, flags cross-target promiscuous binders,
and outputs a ranked shortlist for MM-GBSA (top 5–10 per target).

Runs locally (no cluster required).

Usage:
    python 16_consensus_filter.py

Outputs:
    results/phase4/consensus_hits.csv        — all passing compounds per target
    results/phase4/candidates_shortlist.csv  — top 5-10 per target for MM-GBSA
    results/phase4/consensus_summary.txt     — human-readable report
"""

import csv
import os
from pathlib import Path
from collections import defaultdict

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
PHASE4   = PROJECT / "results/phase4"

# ── Per-target thresholds ─────────────────────────────────────────────────────
# Standard targets: Vina ≤ -6.5 AND CNN ≥ 0.70
# BHLHE41: shallow bHLH surface — lower Vina threshold to -5.5
# IKZF1:   CRBN pocket; CNN is reliable but Vina underestimates — CNN ≥ 0.85 only
# MAF:     top score -5.33 (Vina), shallow bZIP — use CNN ≥ 0.72 as sole criterion
THRESHOLDS = {
    "PPARG_1FM9_LBD_prep":   {"vina": -6.5,  "cnn": 0.70},
    "IRF8_AF2_DBD_prep":     {"vina": -6.5,  "cnn": 0.70},
    "MAF_4EOT_bZIP_prep":    {"vina":  0.0,  "cnn": 0.75},   # CNN-only; Vina uninformative
    "RUNX1_1LJM_Runt_prep":  {"vina": -5.5,  "cnn": 0.70},
    "BHLHE41_AF2_bHLH_prep": {"vina": -5.5,  "cnn": 0.70},
    # IKZF1: CNN hits identified as false positives (no glutarimide pharmacophore).
    # CNS drug library incompatible with CRBN-IKZF1 molecular glue pocket.
    # Excluded from MM-GBSA. Requires dedicated IMiD-analogue library. See Methods note.
    "IKZF1_8RQC_CRBN":       {"vina":  0.0,  "cnn": 99.0},   # effectively excluded
}

# Shortlist size per target (for MM-GBSA)
SHORTLIST_N = {
    "PPARG_1FM9_LBD_prep":   10,
    "IRF8_AF2_DBD_prep":      5,
    "MAF_4EOT_bZIP_prep":     5,
    "RUNX1_1LJM_Runt_prep":   5,
    "BHLHE41_AF2_bHLH_prep":  5,
    "IKZF1_8RQC_CRBN":        5,
}

# Cross-target promiscuity: consensus in ≥ N targets flags a compound
# (still included but annotated; final call left to researcher)
PROMISCUITY_N = 3


def load_scores(target_stem):
    """Return list of dicts from gnina_scores CSV, or empty if missing."""
    f = PHASE4 / f"gnina_scores_{target_stem}.csv"
    if not f.exists():
        return []
    return list(csv.DictReader(open(f)))


def consensus_score(vina, cnn):
    """
    Composite rank score for sorting within target.
    Normalise both metrics to [0,1] range and combine:
      - Vina contribution: scale -10→0 to 1→0
      - CNN contribution: already 0→1
    Weight: 40% Vina, 60% CNN (CNN more reliable per Wong 2022)
    """
    vina_norm = min(max(-vina / 10.0, 0), 1)   # -10 → 1.0, 0 → 0.0
    return 0.4 * vina_norm + 0.6 * cnn


def main():
    # ── Load all scores ───────────────────────────────────────────────────────
    all_scores = {}   # target → {chembl_id → {vina, cnn, cnn_affinity}}
    for tgt in THRESHOLDS:
        rows = load_scores(tgt)
        all_scores[tgt] = {
            r["chembl_id"]: {
                "vina":         float(r["vina_score"]),
                "cnn":          float(r["cnn_score"]),
                "cnn_affinity": float(r["cnn_affinity"]),
            }
            for r in rows
        }

    # ── Apply thresholds ──────────────────────────────────────────────────────
    passing = {}   # target → [(chembl_id, vina, cnn, composite)]
    for tgt, thresh in THRESHOLDS.items():
        hits = []
        for cid, s in all_scores[tgt].items():
            if s["vina"] <= thresh["vina"] and s["cnn"] >= thresh["cnn"]:
                comp = consensus_score(s["vina"], s["cnn"])
                hits.append((cid, s["vina"], s["cnn"], s["cnn_affinity"], comp))
        hits.sort(key=lambda x: -x[4])  # sort by composite score
        passing[tgt] = hits

    # ── Cross-target promiscuity ──────────────────────────────────────────────
    target_count = defaultdict(int)   # chembl_id → number of targets it passes
    for hits in passing.values():
        for cid, *_ in hits:
            target_count[cid] += 1

    promiscuous = {cid for cid, n in target_count.items() if n >= PROMISCUITY_N}

    # ── Write consensus_hits.csv ──────────────────────────────────────────────
    hits_csv = PHASE4 / "consensus_hits.csv"
    with open(hits_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "chembl_id", "vina_score", "cnn_score",
                    "cnn_affinity", "composite_score", "promiscuous"])
        for tgt, hits in passing.items():
            for cid, vina, cnn, cnnaff, comp in hits:
                w.writerow([tgt, cid, vina, cnn, cnnaff,
                            round(comp, 4), cid in promiscuous])

    # ── Write candidates_shortlist.csv (top-N per target, non-promiscuous first) ──
    shortlist_rows = []
    for tgt, hits in passing.items():
        n = SHORTLIST_N[tgt]
        # Prefer non-promiscuous; pad with promiscuous if needed
        non_prom = [(cid, v, c, ca, comp) for cid, v, c, ca, comp in hits if cid not in promiscuous]
        prom     = [(cid, v, c, ca, comp) for cid, v, c, ca, comp in hits if cid in promiscuous]
        selected = (non_prom + prom)[:n]
        for rank, (cid, vina, cnn, cnnaff, comp) in enumerate(selected, 1):
            shortlist_rows.append({
                "rank":            rank,
                "target":          tgt,
                "chembl_id":       cid,
                "vina_score":      vina,
                "cnn_score":       cnn,
                "cnn_affinity":    cnnaff,
                "composite_score": round(comp, 4),
                "promiscuous":     cid in promiscuous,
                "tier":            "Tier1" if "PPARG" in tgt or "IRF8" in tgt or "MAF" in tgt
                                   else "Tier2",
            })

    shortlist_csv = PHASE4 / "candidates_shortlist.csv"
    with open(shortlist_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=shortlist_rows[0].keys())
        w.writeheader()
        w.writerows(shortlist_rows)

    # ── Human-readable report ─────────────────────────────────────────────────
    lines = []
    lines.append("DAM-DRUG Phase 4 — Consensus Filter Report")
    lines.append(f"Generated: 2026-03-24")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Thresholds applied:")
    for tgt, t in THRESHOLDS.items():
        vina_str = f"Vina ≤ {t['vina']}" if t["vina"] > -9 else "Vina: disabled"
        lines.append(f"  {tgt:35s}  {vina_str}  |  CNN ≥ {t['cnn']}")
    lines.append("")
    lines.append("Scoring: composite = 0.4×(Vina/10) + 0.6×CNN  [CNN-weighted per Wong 2022]")
    lines.append(f"Promiscuity flag: consensus hits in ≥ {PROMISCUITY_N} targets")
    lines.append("")

    total_shortlist = 0
    for tgt, hits in passing.items():
        n = SHORTLIST_N[tgt]
        thresh = THRESHOLDS[tgt]
        vina_str = f"Vina ≤ {thresh['vina']}" if thresh["vina"] > -9 else "Vina: disabled"
        lines.append(f"{'─'*70}")
        lines.append(f"  {tgt}")
        lines.append(f"  Threshold: {vina_str}  |  CNN ≥ {thresh['cnn']}")
        lines.append(f"  Passing: {len(hits)}   Shortlisted for MM-GBSA: {min(n, len(hits))}")
        lines.append("")

        non_prom = [(c, v, cnn, ca, comp) for c, v, cnn, ca, comp in hits if c not in promiscuous]
        prom     = [(c, v, cnn, ca, comp) for c, v, cnn, ca, comp in hits if c in promiscuous]
        selected = (non_prom + prom)[:n]
        total_shortlist += len(selected)

        lines.append(f"  {'#':>3}  {'ChEMBL':20s}  {'Vina':>7}  {'CNN':>7}  {'Comp':>7}  {'Flag'}")
        lines.append(f"  {'---':>3}  {'-'*20}  {'-'*7}  {'-'*7}  {'-'*7}  {'----'}")
        for i, (cid, vina, cnn, _, comp) in enumerate(selected, 1):
            flag = "⚠ promiscuous" if cid in promiscuous else ""
            lines.append(f"  {i:>3}  {cid:20s}  {vina:>7.2f}  {cnn:>7.4f}  {comp:>7.4f}  {flag}")

        # Remaining passing but not shortlisted
        remaining = (non_prom + prom)[n:]
        if remaining:
            lines.append(f"       ... {len(remaining)} more passing compounds not shortlisted")
        lines.append("")

    lines.append("=" * 70)
    lines.append(f"Total shortlisted for MM-GBSA: {total_shortlist} compounds")
    lines.append("")
    lines.append("Notes:")
    lines.append("  - BHLHE41/RUNX1/MAF: low Vina scores reflect shallow pockets; treat")
    lines.append("    with caution — MM-GBSA will arbitrate.")
    lines.append("  - IKZF1: EXCLUDED from MM-GBSA. All 4 CNN hits (CHEMBL278819/minaprine,")
    lines.append("    CHEMBL1214124/perampanel, CHEMBL220492, CHEMBL1201/thiothixene) lack the")
    lines.append("    glutarimide pharmacophore required for CRBN binding. CNS drug library is")
    lines.append("    pharmacophore-incompatible with CRBN-IKZF1 molecular glue pocket.")
    lines.append("    Requires dedicated IMiD-analogue library or de novo design (DeLinker).")
    lines.append("  - Promiscuous compounds flagged but not excluded; selectivity")
    lines.append("    analysis (Step 4.12) required before advancement.")
    lines.append("  - Next step: 17_run_mmpbsa.py (MM-GBSA on shortlisted compounds)")

    report = PHASE4 / "consensus_summary.txt"
    report.write_text("\n".join(lines) + "\n")

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n".join(lines))
    print(f"\nWritten:")
    print(f"  {hits_csv}")
    print(f"  {shortlist_csv}")
    print(f"  {report}")


if __name__ == "__main__":
    main()
