"""
DAM-DRUG Phase 3 — Structure preparation with PDBFixer
========================================================
Mandatory binding-site refinement step before fpocket / docking.

For each trimmed structure:
  1. PDBFixer: add missing residues, add missing heavy atoms, add hydrogens
  2. OpenMM minimization (1000 steps, implicit solvent) to relieve clashes
  3. Write prepared PDB to data/structures/prepared/

Requires: pdbfixer, openmm
  conda install -n base -c conda-forge pdbfixer openmm

Run on compute node (CPU OpenMM):
  conda run -n base python code/phase3_structure/08_prepare_structures.py
"""

import os
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
TRIM_DIR = PROJECT / "data/structures/trimmed"
PREP_DIR = PROJECT / "data/structures/prepared"
PREP_DIR.mkdir(parents=True, exist_ok=True)

try:
    from pdbfixer import PDBFixer
    from openmm import app, unit, LangevinMiddleIntegrator
    from openmm.app import PDBFile, Modeller, ForceField, Simulation
    HAS_OPENMM = True
except ImportError:
    log.warning("openmm not available — will run PDBFixer only (no minimization)")
    HAS_OPENMM = False

try:
    from pdbfixer import PDBFixer
    HAS_PDBFIXER = True
except ImportError:
    log.error("pdbfixer not found — install with: conda install -c conda-forge pdbfixer")
    raise SystemExit(1)


def prepare(pdb_in: Path, pdb_out: Path, ph: float = 7.4):
    """Run PDBFixer + optional OpenMM minimization on a trimmed PDB."""
    if pdb_out.exists() and pdb_out.stat().st_size > 1000:
        log.info(f"  Already prepared: {pdb_out.name}")
        return

    log.info(f"  Preparing {pdb_in.name} ...")

    # ── PDBFixer ────────────────────────────────────────────────────────────────
    fixer = PDBFixer(filename=str(pdb_in))

    # Remove non-standard residues (ligands, ions, solvent) inherited from crystal
    fixer.removeHeterogens(keepWater=False)

    # Add missing residues (loops not resolved in crystal) within existing chains
    fixer.findMissingResidues()
    # Only fill gaps within the chain — skip terminal extensions that would add
    # residues outside the trimmed domain window
    missing = fixer.missingResidues
    chains = list(fixer.topology.chains())
    filtered = {}
    for key, residues in missing.items():
        chain_idx, res_idx = key
        n_chain_res = sum(1 for _ in chains[chain_idx].residues())
        # Skip terminal missing (res_idx == 0 or == n_chain_res)
        if 0 < res_idx < n_chain_res:
            filtered[key] = residues
    fixer.missingResidues = filtered

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    if not HAS_OPENMM:
        # Write PDBFixer output directly (no minimization)
        with open(str(pdb_out), "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        log.info(f"  Saved (no minimization): {pdb_out.name}")
        return

    # ── OpenMM minimization ─────────────────────────────────────────────────────
    # Implicit solvent (GBn2) — no water box needed, fast for small domains
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")

    modeller = Modeller(fixer.topology, fixer.positions)
    try:
        modeller.addHydrogens(ff, pH=ph)
    except Exception as e:
        log.warning(f"  addHydrogens failed ({e}), using PDBFixer hydrogens")

    try:
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )
    except Exception as e:
        log.warning(f"  Force field failed ({e}), writing without minimization")
        with open(str(pdb_out), "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        return

    integrator = LangevinMiddleIntegrator(300 * unit.kelvin,
                                          1 / unit.picosecond,
                                          0.002 * unit.picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    pre_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    log.info(f"    Pre-minimization energy: {pre_energy}")

    simulation.minimizeEnergy(maxIterations=200)

    post_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    log.info(f"    Post-minimization energy: {post_energy}")

    state = simulation.context.getState(getPositions=True)
    with open(str(pdb_out), "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    log.info(f"  Saved (minimized): {pdb_out.name}")


# ── Run on all trimmed structures ──────────────────────────────────────────────
pdb_files = sorted(TRIM_DIR.glob("*.pdb"))
log.info(f"Found {len(pdb_files)} trimmed structures")

failed = []
for pdb_in in pdb_files:
    out_name = pdb_in.stem + "_prep.pdb"
    pdb_out = PREP_DIR / out_name
    try:
        prepare(pdb_in, pdb_out)
    except Exception as e:
        log.error(f"  FAILED {pdb_in.name}: {e}")
        failed.append(pdb_in.name)

log.info("\n=== Prepared files ===")
for f in sorted(PREP_DIR.glob("*.pdb")):
    log.info(f"  {f.name}  {f.stat().st_size // 1024} KB")

if failed:
    log.error(f"Failed: {failed}")
    raise SystemExit(1)

log.info("Done. Ready for 09_run_fpocket.sh")
