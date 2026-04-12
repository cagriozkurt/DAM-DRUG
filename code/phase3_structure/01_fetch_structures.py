"""
DAM-DRUG Phase 3 — Fetch protein structures
============================================
Downloads AlphaFold2 predicted structures from EBI and key experimental
structures from RCSB PDB for the 7 priority TF targets + 3 additional
targets added in step 3.2b: SPI1, ACSL1, PIK3CA.

Run on login node (no SLURM needed — small files, fast):
  conda run -n scenic python code/phase3_structure/06_fetch_structures.py
"""

import os
import json
import logging
import urllib.request
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
STRUCT  = PROJECT / "data/structures"
AF2_DIR = STRUCT / "af2"
PDB_DIR = STRUCT / "pdb"

for d in [AF2_DIR, PDB_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Target definitions ─────────────────────────────────────────────────────────
# UniProt IDs for AF2 download; PDB IDs for experimental structures
TARGETS = {
    # PPARG: 3 high-quality LBD structures (all confirmed correct)
    "PPARG":   {"uniprot": "P37231",  "pdb": ["2PRG", "1FM9", "4EMA"]},
    # IKZF1: no DNA-binding domain PDB; 6H0F = ZF2 in CRBN/pomalidomide complex (PROTAC track)
    "IKZF1":   {"uniprot": "Q13422",  "pdb": ["6H0F"]},
    # IRF8: no crystal structure in PDB — AF2 only
    "IRF8":    {"uniprot": "Q02556",  "pdb": []},
    # BHLHE41: no usable PDB — AF2 only
    "BHLHE41": {"uniprot": "O14503",  "pdb": []},
    # RUNX1: 4L0Z = RUNX1+ETS1+DNA (likely mouse but highly conserved; caveat in Methods)
    #        1LJM = human RUNX1 Runt domain, 2.50 Å — ideal structure
    "RUNX1":   {"uniprot": "Q01196",  "pdb": ["4L0Z", "1LJM"]},
    # RUNX2: best PDB structures are 4.2 Å (> 3.0 Å threshold) — AF2 only
    "RUNX2":   {"uniprot": "Q13580",  "pdb": []},
    # MAF: no human c-Maf PDB; 4EOT = MafA homodimer (closest large-Maf family PDB, conserved bZIP)
    "MAF":     {"uniprot": "O15525",  "pdb": ["4EOT"]},
    # ── Step 3.2b additions ────────────────────────────────────────────────────
    # SPI1 (PU.1): causally validated DAM driver (Ayata 2025); ETS domain well-defined
    # 8EE9 = human PU.1 ETS domain + DNA (1.22 Å, Cell Reports 2023); chain F = protein
    "SPI1":   {"uniprot": "P17947",  "pdb": ["8EE9"]},
    # ACSL1: LDAM metabolic target (Haney 2024); long-chain acyl-CoA synthetase
    # No mammalian crystal structure in PDB — AF2 model only (P33121)
    "ACSL1":  {"uniprot": "P33121",  "pdb": []},
    # PIK3CA: PI3K-α catalytic subunit; LDAM lipid metabolism (Haney 2024)
    # 4OVU = PIK3CA p110α (chain A) + p85α-nSH2 (chain B), 2.96 Å
    "PIK3CA": {"uniprot": "P42336",  "pdb": ["4OVU"]},
}

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
PDB_URL  = "https://files.rcsb.org/download/{pdbid}.pdb"


def fetch(url: str, dest: Path, label: str):
    if dest.exists() and dest.stat().st_size > 1000:
        log.info(f"  Already exists: {dest.name}")
        return True
    log.info(f"  Downloading {label} ...")
    try:
        urllib.request.urlretrieve(url, dest)
        size_kb = dest.stat().st_size // 1024
        log.info(f"  Saved: {dest.name} ({size_kb} KB)")
        return True
    except Exception as e:
        log.error(f"  FAILED {label}: {e}")
        if dest.exists():
            dest.unlink()
        return False


def get_af2_url(uniprot: str) -> str:
    """Query AFDB API to get the current PDB download URL for a UniProt accession."""
    api_url = AFDB_API.format(uniprot=uniprot)
    try:
        with urllib.request.urlopen(api_url, timeout=15) as resp:
            data = json.loads(resp.read())
        pdb_url = data[0].get("pdbUrl", "")
        if not pdb_url:
            raise ValueError("pdbUrl field missing in AFDB API response")
        log.info(f"  AFDB API → {pdb_url}")
        return pdb_url
    except Exception as e:
        log.error(f"  AFDB API failed for {uniprot}: {e}")
        return ""


failed = []

for tf, info in TARGETS.items():
    log.info(f"\n=== {tf} ===")

    # AF2 — resolve current URL via AFDB API
    af2_path = AF2_DIR / f"{tf}_AF2.pdb"
    af2_url = get_af2_url(info["uniprot"])
    if af2_url:
        ok = fetch(af2_url, af2_path, f"{tf} AF2 (UniProt {info['uniprot']})")
    else:
        ok = False
    if not ok:
        failed.append(f"{tf} AF2")

    # PDB experimentals
    for pdbid in info["pdb"]:
        pdb_path = PDB_DIR / f"{pdbid}.pdb"
        url = PDB_URL.format(pdbid=pdbid)
        ok = fetch(url, pdb_path, f"{pdbid} ({tf})")
        if not ok:
            failed.append(pdbid)

log.info("\n=== Summary ===")
all_files = list(AF2_DIR.glob("*.pdb")) + list(PDB_DIR.glob("*.pdb"))
log.info(f"  Total files: {len(all_files)}")
for f in sorted(all_files):
    log.info(f"  {f.relative_to(STRUCT)}  {f.stat().st_size//1024} KB")

if failed:
    log.error(f"  FAILED: {failed}")
    raise SystemExit(1)
else:
    log.info("  All downloads succeeded.")
