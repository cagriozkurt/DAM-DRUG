"""
DAM-DRUG Phase 4 — Docking receptor preparation
=================================================
For each docking target (Tier 1 + Tier 2 + IKZF1 PROTAC):
  1. Extract pocket centroid from fpocket pocket1_atm.pdb
  2. Compute Vina grid box size from pocket volume
  3. Convert prepared PDB → rigid PDBQT (obabel, Gasteiger charges)
  4. Write AutoDock Vina config file

IKZF1 PROTAC track (8RQC):
  5. Download 8RQC (CRBN-mezigdomide-IKZF1_ZF2, 2.15 Å X-ray)
  6. Extract CRBN chain + IKZF1_ZF2 chain as receptor (ternary complex approach)
  7. Define box from mezigdomide (QFC) ligand centroid
  8. Strip ligand from receptor PDBQT (re-dock new compounds)

Outputs:
  data/docking/receptors/<stem>.pdbqt
  data/docking/configs/<stem>.conf
  results/phase4/docking_targets.csv
"""

import os
import math
import logging
import subprocess
import urllib.request
import pandas as pd
import numpy as np
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT   = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
FPOCKET   = PROJECT / "results/phase3/fpocket"
PREPARED  = PROJECT / "data/structures/prepared"
RECV_DIR  = PROJECT / "data/docking/receptors"
CONF_DIR  = PROJECT / "data/docking/configs"
RESULT    = PROJECT / "results/phase4"
POCKET_CSV = PROJECT / "results/phase3/pocket_summary.csv"

for d in [RECV_DIR, CONF_DIR, RESULT]:
    d.mkdir(parents=True, exist_ok=True)

OBABEL = os.environ.get("OBABEL", "obabel")

# ── Docking targets ───────────────────────────────────────────────────────────
# (prepared_stem, tier, docking_track)
TARGETS = [
    ("PPARG_1FM9_LBD_prep",   "Tier1", "standard"),
    ("IRF8_AF2_DBD_prep",     "Tier1", "standard"),
    ("MAF_4EOT_bZIP_prep",    "Tier1", "standard"),
    ("RUNX1_1LJM_Runt_prep",  "Tier2", "standard"),
    ("BHLHE41_AF2_bHLH_prep", "Tier2", "standard"),
]
# IKZF1 PROTAC track handled separately below (8RQC)
# PPARG_4EMA and AF2 structures: secondary; skip for initial screen


# ── Helper functions ──────────────────────────────────────────────────────────

def get_pocket_centroid(pocket_atm: Path):
    """Centroid of pocket lining atoms from fpocket pocket1_atm.pdb."""
    coords = []
    for line in pocket_atm.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except (ValueError, IndexError):
                pass
    if not coords:
        return None
    xs, ys, zs = zip(*coords)
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))


def box_size_from_volume(vol_A3: float) -> int:
    """Cube side (Å) for docking grid. Radius from sphere-equiv volume + 10Å padding."""
    r = (3 * vol_A3 / (4 * math.pi)) ** (1 / 3)
    side = 2 * r + 10          # 5Å padding each side
    side = max(15.0, min(30.0, side))
    return int(round(side / 2) * 2)   # round to nearest even Å


def prepare_receptor_pdbqt(pdb_in: Path, pdbqt_out: Path) -> bool:
    """Convert prepared PDB → rigid PDBQT (obabel, Gasteiger charges)."""
    cmd = [OBABEL, str(pdb_in), "-O", str(pdbqt_out),
           "-xr",                         # rigid receptor
           "--partialcharge", "gasteiger",
           "-h"]                          # keep/add hydrogens
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not pdbqt_out.exists():
        log.warning(f"  obabel failed for {pdb_in.name}: {result.stderr.strip()[:200]}")
        return False
    return True


def write_vina_config(conf_path: Path, pdbqt_path: Path,
                      cx: float, cy: float, cz: float, box: int):
    conf_path.write_text(
        f"receptor = {pdbqt_path}\n"
        f"center_x = {cx:.3f}\n"
        f"center_y = {cy:.3f}\n"
        f"center_z = {cz:.3f}\n"
        f"size_x = {box}\n"
        f"size_y = {box}\n"
        f"size_z = {box}\n"
        f"exhaustiveness = 32\n"
        f"num_modes = 10\n"
        f"energy_range = 4\n"
    )


# ── Load pocket summary for volumes ──────────────────────────────────────────
pocket_df = pd.read_csv(POCKET_CSV).set_index("structure_id")


# ── Standard targets ─────────────────────────────────────────────────────────
records = []

for stem, tier, track in TARGETS:
    log.info(f"\n--- {stem} ({tier}) ---")

    pdb_in   = PREPARED / f"{stem}.pdb"
    pdbqt_out = RECV_DIR / f"{stem}.pdbqt"
    conf_out  = CONF_DIR / f"{stem}.conf"
    pocket_atm = FPOCKET / f"{stem}_out" / "pockets" / "pocket1_atm.pdb"

    if not pdb_in.exists():
        log.error(f"  Prepared PDB not found: {pdb_in}")
        continue

    # Centroid
    if not pocket_atm.exists():
        log.warning(f"  pocket1_atm.pdb not found for {stem} — skipping")
        continue
    centroid = get_pocket_centroid(pocket_atm)
    if centroid is None:
        log.warning(f"  Could not compute centroid for {stem}")
        continue
    cx, cy, cz = centroid

    # Box size
    if stem in pocket_df.index:
        vol = pocket_df.loc[stem, "volume_A3"]
        box = box_size_from_volume(float(vol))
    else:
        box = 20
        log.warning(f"  Volume not in pocket_summary — using default 20Å box")

    log.info(f"  Centroid: ({cx:.2f}, {cy:.2f}, {cz:.2f})  Box: {box}Å")

    # Receptor PDBQT
    if pdbqt_out.exists():
        log.info(f"  PDBQT exists, skipping obabel")
    else:
        ok = prepare_receptor_pdbqt(pdb_in, pdbqt_out)
        if ok:
            log.info(f"  Receptor PDBQT: {pdbqt_out.name}  ({pdbqt_out.stat().st_size // 1024} KB)")
        else:
            log.error(f"  PDBQT conversion failed for {stem}")
            continue

    # Config
    write_vina_config(conf_out, pdbqt_out, cx, cy, cz, box)
    log.info(f"  Config: {conf_out.name}")

    records.append({
        "stem": stem, "tier": tier, "track": track,
        "center_x": round(cx, 3), "center_y": round(cy, 3), "center_z": round(cz, 3),
        "box_size": box,
        "volume_A3": float(pocket_df.loc[stem, "volume_A3"]) if stem in pocket_df.index else None,
        "drug_score": float(pocket_df.loc[stem, "drug_score"]) if stem in pocket_df.index else None,
    })


# ── IKZF1 PROTAC track — 8RQC ────────────────────────────────────────────────
log.info("\n--- IKZF1 PROTAC track (8RQC) ---")

pdb_8rqc = PROJECT / "data/structures/pdb" / "8RQC.pdb"
pdb_8rqc.parent.mkdir(parents=True, exist_ok=True)
if not pdb_8rqc.exists():
    log.info("  Downloading 8RQC from RCSB...")
    try:
        url = "https://files.rcsb.org/download/8RQC.pdb"
        urllib.request.urlretrieve(url, pdb_8rqc)
        log.info(f"  Downloaded: {pdb_8rqc.name}  ({pdb_8rqc.stat().st_size // 1024} KB)")
    except Exception as e:
        log.error(f"  Download failed: {e}")
        pdb_8rqc = None

if pdb_8rqc and pdb_8rqc.exists():
    # Identify mezigdomide (QFC) ligand centroid → defines box center
    lig_coords = []
    crbn_chains = set()
    ikzf1_chains = set()
    for line in pdb_8rqc.read_text().splitlines():
        if line.startswith("HETATM") and "QFC" in line[17:20]:
            try:
                lig_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except (ValueError, IndexError):
                pass
        if line.startswith("COMPND"):
            if "CEREBLON" in line.upper():
                pass  # chain ID parsed from ATOM records below
        if line.startswith("ATOM"):
            chain = line[21]
            # Will identify CRBN and IKZF1 chains from SEQRES/COMPND records after loop

    if not lig_coords:
        log.warning("  QFC ligand not found in 8RQC — check chain/residue name")
    else:
        xs, ys, zs = zip(*lig_coords)
        lig_cx = sum(xs)/len(xs)
        lig_cy = sum(ys)/len(ys)
        lig_cz = sum(zs)/len(zs)
        lig_box = 20  # mezigdomide binding site is compact
        log.info(f"  QFC ligand centroid: ({lig_cx:.2f}, {lig_cy:.2f}, {lig_cz:.2f})  Box: {lig_box}Å")

        # Write config (receptor PDBQT for 8RQC prepared separately — see note)
        pdbqt_8rqc = RECV_DIR / "IKZF1_8RQC_CRBN_prep.pdbqt"
        conf_8rqc  = CONF_DIR  / "IKZF1_8RQC_CRBN.conf"
        write_vina_config(conf_8rqc, pdbqt_8rqc, lig_cx, lig_cy, lig_cz, lig_box)
        log.info(f"  Config written: {conf_8rqc.name}")
        log.info("  NOTE: Run 02_ikzf1_protac_prep.py to extract CRBN receptor and prepare PDBQT")

        records.append({
            "stem": "IKZF1_8RQC_CRBN", "tier": "PROTAC", "track": "molecular_glue",
            "center_x": round(lig_cx, 3), "center_y": round(lig_cy, 3),
            "center_z": round(lig_cz, 3), "box_size": lig_box,
            "volume_A3": None, "drug_score": None,
        })


# ── Summary ───────────────────────────────────────────────────────────────────
df = pd.DataFrame(records)
out_csv = RESULT / "docking_targets.csv"
df.to_csv(out_csv, index=False)

log.info(f"\n{'='*70}")
log.info("DOCKING TARGETS SUMMARY")
log.info(f"{'='*70}")
log.info(df[["stem", "tier", "center_x", "center_y", "center_z", "box_size", "drug_score"]].to_string(index=False))
log.info(f"\nSaved → {out_csv}")
log.info("\nNext step: 03_prep_ligands.py (download + filter + 3D-conform Repurposing Hub library)")
