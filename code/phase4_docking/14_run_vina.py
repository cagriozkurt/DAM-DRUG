"""
DAM-DRUG Phase 4 — AutoDock Vina docking (single target)
=========================================================
Called by 14_run_vina.slurm for each docking target.
Runs Vina on all ligand PDBQTs in parallel using multiprocessing.

Usage:
    python 14_run_vina.py <target_stem> [--cpus N]

Example:
    python 14_run_vina.py PPARG_1FM9_LBD_prep --cpus 20

Output per compound:
    results/phase4/vina/<target_stem>/<chembl_id>.pdbqt   (poses)
    results/phase4/vina/<target_stem>/<chembl_id>.log     (scores)
"""

import os
import sys
import argparse
import logging
import subprocess
import multiprocessing
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT   = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
PDBQT_DIR = PROJECT / "data/compounds/pdbqt"
CONF_DIR  = PROJECT / "data/docking/configs"
VINA_OUT  = PROJECT / "results/phase4/vina"
VINA_BIN  = os.environ.get("VINA", str(Path.home() / "apps/vina/vina"))


def run_one(args):
    ligand, conf, out_pdbqt, log_path = args
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
        return ligand.stem, 0, "skip"

    cmd = [VINA_BIN,
           "--config",    str(conf),
           "--ligand",    str(ligand),
           "--out",       str(out_pdbqt),
           "--cpu",       "1"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        # Vina 1.2.x writes scores to stdout; save for later parsing
        log_path.write_text(result.stdout)
        return ligand.stem, result.returncode, result.stderr.strip()[:120]
    except subprocess.TimeoutExpired:
        return ligand.stem, -1, "timeout"
    except Exception as e:
        return ligand.stem, -2, str(e)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_stem", help="Prepared receptor stem (e.g. PPARG_1FM9_LBD_prep)")
    parser.add_argument("--cpus", type=int, default=20)
    args = parser.parse_args()

    stem = args.target_stem
    conf = CONF_DIR / f"{stem}.conf"
    if not conf.exists():
        # Try without _prep suffix for PROTAC target
        conf_alt = CONF_DIR / f"{stem}.conf"
        if not conf_alt.exists():
            log.error(f"Config not found: {conf}")
            sys.exit(1)

    out_dir = VINA_OUT / stem
    out_dir.mkdir(parents=True, exist_ok=True)

    ligands = sorted(PDBQT_DIR.glob("*.pdbqt"))
    if not ligands:
        log.error(f"No ligand PDBQTs found in {PDBQT_DIR}")
        sys.exit(1)

    log.info(f"Target:  {stem}")
    log.info(f"Config:  {conf}")
    log.info(f"Ligands: {len(ligands)}")
    log.info(f"CPUs:    {args.cpus}")
    log.info(f"Output:  {out_dir}")

    task_args = [
        (lig, conf, out_dir / f"{lig.stem}.pdbqt", out_dir / f"{lig.stem}.log")
        for lig in ligands
    ]

    n_ok = n_skip = n_fail = 0
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for i, (chembl_id, rc, msg) in enumerate(pool.imap_unordered(run_one, task_args), 1):
            if rc == 0 and msg == "skip":
                n_skip += 1
            elif rc == 0:
                n_ok += 1
            else:
                n_fail += 1
                log.warning(f"  FAILED {chembl_id}: {msg}")
            if i % 50 == 0:
                log.info(f"  Progress: {i}/{len(ligands)} — ok={n_ok} skip={n_skip} fail={n_fail}")

    log.info(f"\n=== Vina complete: {stem} ===")
    log.info(f"  Docked: {n_ok}  Skipped: {n_skip}  Failed: {n_fail}")

    # Extract best scores from log files
    scores = {}
    for log_file in out_dir.glob("*.log"):
        for line in log_file.read_text().splitlines():
            # Vina log: "   1       -8.5      0.000      0.000"
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "1":
                try:
                    scores[log_file.stem] = float(parts[1])
                    break
                except ValueError:
                    pass

    if scores:
        import csv
        score_csv = PROJECT / "results/phase4" / f"vina_scores_{stem}.csv"
        with open(score_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["chembl_id", "vina_score_kcal_mol"])
            for cid, sc in sorted(scores.items(), key=lambda x: x[1]):
                w.writerow([cid, sc])
        log.info(f"  Scores saved: {score_csv.name}  ({len(scores)} entries)")

        top10 = sorted(scores.items(), key=lambda x: x[1])[:10]
        log.info(f"\n  Top 10 by Vina score:")
        for cid, sc in top10:
            log.info(f"    {cid:20s}  {sc:7.2f} kcal/mol")

    log.info("Next step: 15_run_gnina.py (CNN rescoring of top poses)")


if __name__ == "__main__":
    main()
