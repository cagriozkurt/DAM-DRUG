"""
DAM-DRUG Phase 4 — Gnina CNN rescoring of Vina poses
=====================================================
Takes top-N Vina poses per target, filters promiscuous binders,
rescores with gnina CNN (--cnn_scoring rescore, no new search).

Usage:
    python 05_run_gnina.py <target_stem> [--cpus N] [--top N]

Example:
    python 05_run_gnina.py PPARG_1FM9_LBD_prep --cpus 20 --top 30

Output per compound:
    results/phase4/gnina/<target_stem>/<chembl_id>.pdbqt   (rescored poses)
    results/phase4/gnina_scores_<target_stem>.csv           (ranked by CNN score)
"""

import os
import sys
import csv
import re
import argparse
import logging
import subprocess
import multiprocessing
import tempfile
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT      = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
VINA_OUT     = PROJECT / "results/phase4/vina"
GNINA_OUT    = PROJECT / "results/phase4/gnina"
RECEPTOR_DIR = PROJECT / "data/docking/receptors"
CONF_DIR     = PROJECT / "data/docking/configs"
SCORE_DIR    = PROJECT / "results/phase4"
GNINA_BIN    = os.environ.get("GNINA", "gnina")

# Receptor PDBQT stem → filename mapping
RECEPTOR_MAP = {
    "PPARG_1FM9_LBD_prep":   "PPARG_1FM9_LBD_prep.pdbqt",
    "IRF8_AF2_DBD_prep":     "IRF8_AF2_DBD_prep.pdbqt",
    "MAF_4EOT_bZIP_prep":    "MAF_4EOT_bZIP_prep.pdbqt",
    "RUNX1_1LJM_Runt_prep":  "RUNX1_1LJM_Runt_prep.pdbqt",
    "BHLHE41_AF2_bHLH_prep": "BHLHE41_AF2_bHLH_prep.pdbqt",
    "IKZF1_8RQC_CRBN":       "IKZF1_8RQC_CRBN_prep.pdbqt",
}

# Promiscuous binders: top-20 in ≥4 of 6 targets → reject
BLACKLIST = {
    "CHEMBL85",       # risperidone — 6/6 targets
    "CHEMBL1108",     # 5/6
    "CHEMBL3039520",  # 5/6
    "CHEMBL4105630",  # 5/6
    "CHEMBL502",      # 5/6
    "CHEMBL54",       # 4/6
    "CHEMBL1621",     # 4/6
    "CHEMBL708",      # 4/6
    "CHEMBL3306803",  # 4/6
    "CHEMBL1510",     # 4/6
}


def parse_gnina_scores(stdout: str):
    """
    Parse gnina stdout table (format varies by build):

        mode |  affinity  |    CNN     |   CNN
             | (kcal/mol) | pose score | affinity
        -----+------------+------------+----------
            1       14.89       0.8944      5.431

    Trigger on the separator line (-----+) to avoid header-text variation.
    Returns list of (mode, affinity, cnn_pose_score, cnn_affinity) tuples.
    We use cnn_pose_score (col 3) as the primary metric.
    """
    results = []
    in_table = False
    for line in stdout.splitlines():
        # Separator line marks start of data rows
        if re.match(r"\s*-{3,}\+", line):
            in_table = True
            continue
        if not in_table:
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                mode          = int(parts[0])
                affinity      = float(parts[1])
                cnn_score     = float(parts[2])
                cnn_affinity  = float(parts[3])
                results.append((mode, affinity, cnn_score, cnn_affinity))
            except ValueError:
                break  # end of table
        else:
            break
    return results


def parse_conf_box(conf_path: Path):
    """Extract center_x/y/z and size_x/y/z from a Vina conf file."""
    box = {}
    for line in conf_path.read_text().splitlines():
        line = line.strip()
        for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
            if line.startswith(key):
                box[key] = line.split("=")[1].strip()
    return box


def extract_model1(pdbqt_path: Path) -> str:
    """
    Vina 1.2.x wraps each pose in MODEL/ENDMDL with REMARK lines.
    Gnina expects a bare PDBQT (no MODEL/ENDMDL).
    Extract the first model, strip MODEL/ENDMDL/REMARK lines.
    """
    lines = []
    in_model = False
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("MODEL 1"):
            in_model = True
            continue
        if line.startswith("ENDMDL") and in_model:
            break
        if in_model and not line.startswith("REMARK"):
            lines.append(line)
    return "\n".join(lines) + "\n"


def rescore_one(args):
    chembl_id, receptor, in_pdbqt, out_pdbqt, box = args

    # Extract model 1 from Vina multi-model PDBQT into a temp file
    pose_content = extract_model1(in_pdbqt)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False, mode="w")
    tmp.write(pose_content)
    tmp.flush()
    tmp_path = tmp.name
    tmp.close()

    cmd = [
        GNINA_BIN,
        "--receptor",    str(receptor),
        "--ligand",      tmp_path,
        "--cnn_scoring", "rescore",
        "--cnn",         "crossdock_default2018",  # single model, fast on CPU
        "--no_gpu",
        "--center_x",    box["center_x"],
        "--center_y",    box["center_y"],
        "--center_z",    box["center_z"],
        "--size_x",      box["size_x"],
        "--size_y",      box["size_y"],
        "--size_z",      box["size_z"],
        "--out",         str(out_pdbqt),
        "--cpu",         "1",
        "--seed",        "42",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        scores = parse_gnina_scores(result.stdout)
        if not scores:
            return chembl_id, None, None, f"no scores parsed: {result.stderr[:80]}"
        # Best mode = highest CNN score
        best = max(scores, key=lambda x: x[2])
        return chembl_id, best[2], best[3], "ok"
    except subprocess.TimeoutExpired:
        return chembl_id, None, None, "timeout"
    except Exception as e:
        return chembl_id, None, None, str(e)[:80]
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_stem", help="e.g. PPARG_1FM9_LBD_prep")
    parser.add_argument("--cpus", type=int, default=20)
    parser.add_argument("--top",  type=int, default=30,
                        help="Number of top Vina hits to rescore per target")
    args = parser.parse_args()

    stem = args.target_stem

    # Receptor PDBQT
    rec_name = RECEPTOR_MAP.get(stem)
    if rec_name is None:
        log.error(f"Unknown target stem: {stem}")
        sys.exit(1)
    receptor = RECEPTOR_DIR / rec_name
    if not receptor.exists():
        log.error(f"Receptor PDBQT not found: {receptor}")
        sys.exit(1)

    # Vina scores CSV
    vina_csv = SCORE_DIR / f"vina_scores_{stem}.csv"
    if not vina_csv.exists():
        log.error(f"Vina scores not found: {vina_csv}")
        sys.exit(1)

    # Load top-N, filter blacklist
    with open(vina_csv) as f:
        rows = list(csv.DictReader(f))

    candidates = [
        r for r in rows
        if r["chembl_id"] not in BLACKLIST
    ][:args.top]

    log.info(f"Target:    {stem}")
    log.info(f"Receptor:  {receptor.name}")
    log.info(f"Rescoring: {len(candidates)} compounds (top {args.top} minus {len(rows[:args.top]) - len(candidates)} blacklisted)")
    log.info(f"CPUs:      {args.cpus}")

    out_dir = GNINA_OUT / stem
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load box from Vina conf
    conf = CONF_DIR / f"{stem}.conf"
    if not conf.exists():
        log.error(f"Vina conf not found: {conf}")
        sys.exit(1)
    box = parse_conf_box(conf)
    log.info(f"Box:       center=({box['center_x']}, {box['center_y']}, {box['center_z']})  "
             f"size={box['size_x']}Å")

    # Check Vina pose files exist
    vina_dir = VINA_OUT / stem
    task_args = []
    missing = []
    for row in candidates:
        cid = row["chembl_id"]
        in_pdbqt = vina_dir / f"{cid}.pdbqt"
        if not in_pdbqt.exists() or in_pdbqt.stat().st_size == 0:
            missing.append(cid)
            continue
        out_pdbqt = out_dir / f"{cid}.pdbqt"
        task_args.append((cid, receptor, in_pdbqt, out_pdbqt, box))

    if missing:
        log.warning(f"  {len(missing)} Vina pose files missing/empty — skipping: {missing[:5]}")

    # Build vina score lookup for output CSV
    vina_lookup = {r["chembl_id"]: float(r["vina_score_kcal_mol"]) for r in rows}

    # Run gnina rescoring in parallel
    results = []
    n_ok = n_fail = 0
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for i, (cid, cnn_score, cnn_affinity, msg) in enumerate(
            pool.imap_unordered(rescore_one, task_args), 1
        ):
            if cnn_score is not None:
                n_ok += 1
                results.append({
                    "chembl_id":     cid,
                    "vina_score":    vina_lookup.get(cid, ""),
                    "cnn_score":     round(cnn_score, 4),
                    "cnn_affinity":  round(cnn_affinity, 3),
                })
            else:
                n_fail += 1
                log.warning(f"  FAILED {cid}: {msg}")
            if i % 10 == 0:
                log.info(f"  Progress: {i}/{len(task_args)} — ok={n_ok} fail={n_fail}")

    log.info(f"\n=== Gnina complete: {stem} ===")
    log.info(f"  Rescored: {n_ok}  Failed: {n_fail}")

    if not results:
        log.error("No results — check gnina binary and receptor PDBQT")
        sys.exit(1)

    # Sort by CNN score descending
    results.sort(key=lambda x: -x["cnn_score"])

    out_csv = SCORE_DIR / f"gnina_scores_{stem}.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["chembl_id", "vina_score", "cnn_score", "cnn_affinity"])
        w.writeheader()
        w.writerows(results)
    log.info(f"  Scores saved: {out_csv.name}  ({len(results)} entries)")

    log.info(f"\n  Top 10 by CNN score:")
    for r in results[:10]:
        log.info(f"    {r['chembl_id']:20s}  CNN={r['cnn_score']:.4f}  Vina={r['vina_score']:6.2f}  CNNaff={r['cnn_affinity']:.2f}")

    log.info("Next step: 21_run_gnina.slurm → 06_consensus_filter.py")


if __name__ == "__main__":
    main()
