"""
DAM-DRUG Phase 4 — IKZF1 PROTAC track: 8RQC receptor preparation
==================================================================
8RQC = CRBN-mezigdomide-IKZF1_ZF2 complex (2.15 Å X-ray, 2024 Nat Commun)

Strategy: ternary complex docking
  - Receptor = CRBN + IKZF1_ZF2 (both chains retained; molecular glue binds at interface)
  - Ligand site = mezigdomide (QFC) binding pocket
  - Remove QFC ligand from receptor before docking (re-dock new compounds)
  - Keep Zn2+ ions (structural; important for ZF2 fold)

Steps:
  1. Parse 8RQC chain assignments from COMPND records
  2. Extract CRBN chain(s) + IKZF1 ZF2 chain(s); remove DDB1
  3. Remove QFC ligand; keep ZN ions and other cofactors
  4. Run PDBFixer (add missing atoms, H at pH 7.4)
  5. Convert to PDBQT (obabel)
  6. Update IKZF1_8RQC_CRBN.conf with final receptor PDBQT path

Outputs:
  data/structures/pdb/8RQC_CRBN_ZF2.pdb         (extracted chains)
  data/structures/prepared/8RQC_CRBN_ZF2_prep.pdb  (PDBFixer output)
  data/docking/receptors/IKZF1_8RQC_CRBN_prep.pdbqt
"""

import os
import re
import logging
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
PDB_DIR  = PROJECT / "data/structures/pdb"
PREP_DIR = PROJECT / "data/structures/prepared"
RECV_DIR = PROJECT / "data/docking/receptors"
CONF_DIR = PROJECT / "data/docking/configs"

pdb_8rqc = PDB_DIR / "8RQC.pdb"
OBABEL = os.environ.get("OBABEL", "obabel")

if not pdb_8rqc.exists():
    log.error(f"8RQC.pdb not found at {pdb_8rqc}. Run 01_docking_prep.py first.")
    raise SystemExit(1)

# ── Parse COMPND to identify chain→molecule mapping ──────────────────────────
log.info("Parsing 8RQC chain assignments from COMPND records...")

lines = pdb_8rqc.read_text().splitlines()

# Parse COMPND block
compnd_block = "\n".join(l[10:].strip() for l in lines if l.startswith("COMPND"))
log.info(f"COMPND:\n{compnd_block}")

# Parse chain assignments: CHAIN: A, B; pattern
chain_mol = {}   # chain → molecule name
current_mol = None
for l in lines:
    if not l.startswith("COMPND"):
        continue
    content = l[10:].strip()
    if "MOL_ID:" in content:
        current_mol = None
    if "MOLECULE:" in content:
        current_mol = content.split("MOLECULE:")[-1].strip().rstrip(";")
    if "CHAIN:" in content and current_mol:
        chains_str = content.split("CHAIN:")[-1].strip().rstrip(";")
        for ch in [c.strip() for c in chains_str.split(",")]:
            chain_mol[ch] = current_mol

log.info("Chain → molecule mapping:")
for ch, mol in sorted(chain_mol.items()):
    log.info(f"  Chain {ch}: {mol}")

# Identify CRBN and IKZF1 chains
crbn_chains  = {ch for ch, mol in chain_mol.items() if "CEREBLON" in mol.upper() or "CRBN" in mol.upper()}
ikzf1_chains = {ch for ch, mol in chain_mol.items() if "IKAROS" in mol.upper() or "IKZF1" in mol.upper()}
ddb1_chains  = {ch for ch, mol in chain_mol.items() if "DDB1" in mol.upper() or "DAMAGE" in mol.upper()}

log.info(f"CRBN chains:  {crbn_chains}")
log.info(f"IKZF1 chains: {ikzf1_chains}")
log.info(f"DDB1 chains:  {ddb1_chains} (will be excluded)")

keep_chains = crbn_chains | ikzf1_chains
if not keep_chains:
    log.warning("Could not auto-detect chains — keeping all non-DDB1 chains")
    all_chains = set(chain_mol.keys())
    keep_chains = all_chains - ddb1_chains

# ── Extract CRBN + IKZF1_ZF2 chains; remove QFC; keep ZN ───────────────────
log.info(f"Keeping chains: {keep_chains}")
log.info("Removing: QFC ligand | DDB1 chains")
log.info("Keeping:  ZN ions (structural in ZF2)")

out_lines = []
skipped_qfc = 0
for line in lines:
    rec = line[:6].strip()

    if rec in ("ATOM", "HETATM"):
        chain = line[21]
        resname = line[17:20].strip()

        # Skip DDB1 chains
        if chain not in keep_chains:
            continue

        # Remove QFC (mezigdomide) — will be re-docked
        if rec == "HETATM" and resname == "QFC":
            skipped_qfc += 1
            continue

        # Keep ZN (structural zinc in ZF2)
        # Keep other small molecules? Remove water (HOH) for cleaner docking
        if rec == "HETATM" and resname == "HOH":
            continue

        out_lines.append(line)

    elif rec in ("TER", "END", "HEADER", "TITLE", "REMARK"):
        out_lines.append(line)

    # Skip SEQRES, CONECT, MASTER etc. for cleaner PDB

log.info(f"Removed {skipped_qfc} QFC atoms")

out_pdb = PDB_DIR / "8RQC_CRBN_ZF2.pdb"
out_pdb.write_text("\n".join(out_lines) + "\n")
log.info(f"Extracted receptor: {out_pdb.name}  ({out_pdb.stat().st_size // 1024} KB)")

# ── PDBFixer: add missing atoms + hydrogens ──────────────────────────────────
log.info("Running PDBFixer on extracted receptor...")

prep_pdb = PREP_DIR / "8RQC_CRBN_ZF2_prep.pdb"

if not HAS_PDBFIXER:
    log.warning("pdbfixer/openmm not available — using extracted PDB directly")
    shutil.copy(out_pdb, prep_pdb)
else:
    try:
        fixer = PDBFixer(filename=str(out_pdb))
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()

        # Skip terminal extensions
        chains = list(fixer.topology.chains())
        filtered = {}
        for key, residues in fixer.missingResidues.items():
            chain_idx, res_idx = key
            n_res = sum(1 for _ in chains[chain_idx].residues())
            if 0 < res_idx < n_res:
                filtered[key] = residues
        fixer.missingResidues = filtered

        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        with open(str(prep_pdb), "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        log.info(f"PDBFixer done: {prep_pdb.name}  ({prep_pdb.stat().st_size // 1024} KB)")

    except Exception as e:
        log.warning(f"PDBFixer failed ({e}) — using extracted PDB directly")
        shutil.copy(out_pdb, prep_pdb)

# ── Convert to PDBQT ─────────────────────────────────────────────────────────
pdbqt_out = RECV_DIR / "IKZF1_8RQC_CRBN_prep.pdbqt"
cmd = [OBABEL, str(prep_pdb), "-O", str(pdbqt_out), "-xr",
       "--partialcharge", "gasteiger", "-h"]
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode == 0 and pdbqt_out.exists():
    log.info(f"Receptor PDBQT: {pdbqt_out.name}  ({pdbqt_out.stat().st_size // 1024} KB)")
else:
    log.error(f"obabel failed: {result.stderr.strip()[:300]}")
    raise SystemExit(1)

# ── Update Vina config to point to new PDBQT ─────────────────────────────────
conf_path = CONF_DIR / "IKZF1_8RQC_CRBN.conf"
if conf_path.exists():
    conf_text = conf_path.read_text()
    conf_text = re.sub(r"^receptor = .*$", f"receptor = {pdbqt_out}", conf_text, flags=re.MULTILINE)
    conf_path.write_text(conf_text)
    log.info(f"Updated config: {conf_path.name}")

log.info("\n=== IKZF1 PROTAC receptor ready ===")
log.info(f"  Chains retained: {keep_chains}")
log.info(f"  QFC removed (re-dock site)")
log.info(f"  ZN ions kept (ZF2 structural)")
log.info(f"  Receptor PDBQT: {pdbqt_out}")
log.info("Next: 03_prep_ligands.py → 20_run_vina.slurm")
