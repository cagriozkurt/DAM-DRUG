"""
DAM-DRUG Phase 4 — Ligand library preparation
===============================================
Downloads and prepares the Tier 1 CNS drug library for docking.

Source: ChEMBL REST API
  - Approved (max_phase = 4) small molecules
  - ATC level-1 code = N (nervous system)
  - MW <= 500, HBD <= 5, HBA <= 10, logP <= 5 (Lipinski Ro5)
  - MW >= 150 (remove fragments)

Pipeline per compound:
  1. SMILES → 3D conformer (obabel --gen3d -p 7.4)
  2. 3D SDF → PDBQT (meeko mk_prepare_ligand.py)

Outputs:
  data/compounds/tier1_cns_approved.smi      (raw SMILES from ChEMBL)
  data/compounds/tier1_cns_3d.sdf            (3D conformers, pH 7.4)
  data/compounds/pdbqt/                       (one .pdbqt per compound)
  data/compounds/tier1_cns_approved.csv      (metadata: name, MW, logP, etc.)

Run on login node (internet required for ChEMBL API).
Conformer generation via obabel (available at ~/.local/bin/obabel).
PDBQT conversion via meeko (pip install --user meeko).
"""

import os
import time
import json
import logging
import subprocess
import urllib.request
import urllib.parse
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT   = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
COMP_DIR  = PROJECT / "data/compounds"
PDBQT_DIR = COMP_DIR / "pdbqt"
for d in [COMP_DIR, PDBQT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

OBABEL = os.environ.get("OBABEL", "obabel")
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

# ── ChEMBL download ───────────────────────────────────────────────────────────

def chembl_get(url: str, retries: int = 3) -> dict:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read())
        except Exception as e:
            if attempt == retries - 1:
                raise
            log.warning(f"  Retry {attempt+1}: {e}")
            time.sleep(2)


def fetch_cns_drugs() -> list[dict]:
    """
    Fetch approved CNS drugs (ATC N*) from ChEMBL.
    Returns list of dicts with: chembl_id, name, smiles, mw, hbd, hba, logp, alogp.
    """
    log.info("Fetching approved CNS drugs from ChEMBL...")

    # Step 1: get all ATC codes starting with N (nervous system)
    # These are ATC level-1 = N; we query level1_desc = NERVOUS SYSTEM
    records = []
    offset  = 0
    limit   = 1000

    while True:
        params = urllib.parse.urlencode({
            "max_phase": 4,
            "molecule_type": "Small molecule",
            "limit": limit,
            "offset": offset,
            "format": "json",
        })
        url = f"{CHEMBL_API}/molecule.json?{params}"
        data = chembl_get(url)
        molecules = data.get("molecules", [])
        if not molecules:
            break

        for mol in molecules:
            # Filter: must have ATC N* code
            atcs = mol.get("atc_classifications") or []
            has_cns = any(a.startswith("N") for a in atcs)
            if not has_cns:
                continue

            props = mol.get("molecule_properties") or {}
            mw    = float(props.get("mw_freebase") or 0)
            hbd   = int(props.get("hbd") or 0)
            hba   = int(props.get("hba") or 0)
            logp  = float(props.get("alogp") or 0)

            # Lipinski Ro5 filter
            if not (150 <= mw <= 500 and hbd <= 5 and hba <= 10 and logp <= 5):
                continue

            smiles = (mol.get("molecule_structures") or {}).get("canonical_smiles")
            if not smiles:
                continue

            records.append({
                "chembl_id": mol["molecule_chembl_id"],
                "name":      mol.get("pref_name") or mol["molecule_chembl_id"],
                "smiles":    smiles,
                "mw":        round(mw, 2),
                "hbd":       hbd,
                "hba":       hba,
                "alogp":     round(logp, 2),
                "atc":       ";".join(atcs),
            })

        total = data.get("page_meta", {}).get("total_count", 0)
        offset += limit
        log.info(f"  Fetched {offset}/{total} molecules; {len(records)} CNS approved so far")

        if offset >= total:
            break
        time.sleep(0.3)   # be polite to ChEMBL API

    log.info(f"Total CNS approved drugs after Ro5 filter: {len(records)}")
    return records


# ── Conformer generation ──────────────────────────────────────────────────────

def smiles_to_3d(smiles: str, name: str, out_sdf: Path) -> bool:
    """SMILES → 3D conformer at pH 7.4 via obabel."""
    cmd = [OBABEL, "-ismi", f"-:{smiles}", "-osdf", "-O", str(out_sdf),
           "--gen3d", "-p", "7.4", "--title", name, "-h"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0 and out_sdf.exists() and out_sdf.stat().st_size > 0


def sdf_to_pdbqt(sdf_path: Path, pdbqt_path: Path) -> bool:
    """3D SDF → PDBQT via meeko mk_prepare_ligand.py."""
    mk = Path.home() / ".local/bin/mk_prepare_ligand.py"
    if not mk.exists():
        mk = Path("mk_prepare_ligand.py")  # fallback if in PATH
    cmd = ["python3", str(mk), "-i", str(sdf_path), "-o", str(pdbqt_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0


# ── Main ──────────────────────────────────────────────────────────────────────

# 1. Fetch from ChEMBL
smi_out = COMP_DIR / "tier1_cns_approved.smi"
csv_out = COMP_DIR / "tier1_cns_approved.csv"

if csv_out.exists():
    log.info(f"Already have {csv_out.name} — loading")
    import csv
    records = []
    with open(csv_out) as f:
        records = list(csv.DictReader(f))
    log.info(f"  Loaded {len(records)} compounds")
else:
    records = fetch_cns_drugs()
    if not records:
        log.error("No compounds fetched — check ChEMBL API connectivity")
        raise SystemExit(1)

    # Write SMILES file and CSV
    with open(smi_out, "w") as f:
        for r in records:
            f.write(f"{r['smiles']}\t{r['name']}\n")

    import csv
    with open(csv_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=records[0].keys())
        w.writeheader()
        w.writerows(records)

    log.info(f"Saved {len(records)} compounds → {csv_out.name}")


# 2. 3D conformers via obabel
log.info("\n=== Generating 3D conformers (obabel) ===")
sdf_dir  = COMP_DIR / "sdf_3d"
sdf_dir.mkdir(exist_ok=True)

n_done = n_fail = 0
for rec in records:
    cid  = rec["chembl_id"]
    name = rec["name"]
    smi  = rec["smiles"]
    sdf  = sdf_dir / f"{cid}.sdf"

    if sdf.exists() and sdf.stat().st_size > 0:
        n_done += 1
        continue

    ok = smiles_to_3d(smi, name, sdf)
    if ok:
        n_done += 1
    else:
        n_fail += 1
        log.warning(f"  3D gen failed: {cid} ({name})")

    if (n_done + n_fail) % 50 == 0:
        log.info(f"  3D conformers: {n_done} done, {n_fail} failed")

log.info(f"3D conformers complete: {n_done} ok, {n_fail} failed")


# 3. PDBQT conversion via meeko
log.info("\n=== Converting to PDBQT (meeko) ===")
n_pdbqt_ok = n_pdbqt_fail = 0
for sdf in sorted(sdf_dir.glob("*.sdf")):
    pdbqt = PDBQT_DIR / f"{sdf.stem}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        n_pdbqt_ok += 1
        continue

    ok = sdf_to_pdbqt(sdf, pdbqt)
    if ok:
        n_pdbqt_ok += 1
    else:
        n_pdbqt_fail += 1
        log.warning(f"  PDBQT failed: {sdf.stem}")

    if (n_pdbqt_ok + n_pdbqt_fail) % 50 == 0:
        log.info(f"  PDBQT: {n_pdbqt_ok} done, {n_pdbqt_fail} failed")

log.info(f"PDBQT complete: {n_pdbqt_ok} ok, {n_pdbqt_fail} failed")
log.info(f"\nLigand library ready: {n_pdbqt_ok} compounds in {PDBQT_DIR}")
log.info("Next step: 14_run_vina.slurm")
