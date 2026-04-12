"""
DAM-DRUG Phase 3 — Trim structures to druggable domains
=========================================================
Extracts relevant domain chains/residue ranges from PDB and AF2 structures.
fpocket gives cleaner results on isolated domains than full-length proteins.

Targets: PPARG, IKZF1, IRF8, BHLHE41, RUNX1, RUNX2, MAF (7 original)
         + SPI1, ACSL1, PIK3CA (3 added in step 3.2b)

Requires: biopython (available in scanpy-env)

Run on login node or compute node:
  apptainer exec containers/scenic.sif python code/phase3_structure/02_trim_domains.py
"""

import os
import logging
import numpy as np
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBIO, Select

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
STRUCT   = PROJECT / "data/structures"
AF2_DIR  = STRUCT / "af2"
PDB_DIR  = STRUCT / "pdb"
TRIM_DIR = STRUCT / "trimmed"
TRIM_DIR.mkdir(parents=True, exist_ok=True)

parser = PDB.PDBParser(QUIET=True)
io     = PDBIO()

# ── Selector classes ───────────────────────────────────────────────────────────
class DomainSelect(Select):
    """Keep ATOM records for specified chains + residue range. Strip HETATM."""
    def __init__(self, chains, res_min=None, res_max=None, min_plddt=None):
        self.chains   = set(chains)
        self.res_min  = res_min
        self.res_max  = res_max
        self.min_plddt = min_plddt  # for AF2: filter by B-factor (pLDDT)

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        if residue.id[0] != " ":   # skip HETATM residues
            return False
        rid = residue.id[1]
        if self.res_min is not None and rid < self.res_min:
            return False
        if self.res_max is not None and rid > self.res_max:
            return False
        return True

    def accept_atom(self, atom):
        if self.min_plddt is not None:
            return atom.bfactor >= self.min_plddt
        return True

def trim(source_pdb: Path, out_name: str, chains, res_min=None, res_max=None,
         min_plddt=None, notes=""):
    out_path = TRIM_DIR / f"{out_name}.pdb"
    if out_path.exists():
        log.info(f"  Already exists: {out_path.name}")
        return out_path

    if not source_pdb.exists():
        log.warning(f"  Source not found, skipping: {source_pdb}")
        return None

    structure = parser.get_structure(out_name, str(source_pdb))
    model = structure[0]

    sel = DomainSelect(chains, res_min, res_max, min_plddt)
    io.set_structure(model)
    io.save(str(out_path), sel)

    # Report residue count
    n_res = sum(1 for r in model.get_residues()
                if r.id[0] == " "
                and r.parent.id in set(chains)
                and (res_min is None or r.id[1] >= res_min)
                and (res_max is None or r.id[1] <= res_max))
    log.info(f"  {out_path.name}  ({n_res} residues){' — ' + notes if notes else ''}")
    return out_path

# ── Trim definitions ───────────────────────────────────────────────────────────
log.info("=== PPARG ===")
# LBD: residues 230-477, chain A only; strip ligands (HETATM excluded by selector)
for pdbid in ["2PRG", "1FM9", "4EMA"]:
    trim(PDB_DIR / f"{pdbid}.pdb", f"PPARG_{pdbid}_LBD", chains=["A"],
         res_min=230, res_max=477, notes="LBD only, chain A, no ligand")
# AF2: LBD window + full high-confidence
trim(AF2_DIR / "PPARG_AF2.pdb", "PPARG_AF2_LBD", chains=["A"],
     res_min=230, res_max=477, min_plddt=50, notes="AF2 LBD, pLDDT≥50")
trim(AF2_DIR / "PPARG_AF2.pdb", "PPARG_AF2_full", chains=["A"],
     min_plddt=50, notes="AF2 full-length, pLDDT≥50")

log.info("=== IKZF1 ===")
# 6H0F: IKZF1 ZF2 domain in CRBN/DDB1/pomalidomide complex (3.25 Å, human)
# Used for PROTAC/degrader track — ZF2 is the degron surface for cereblon
# Chain assignments: check COMPND; extract IKZF1 chain only, strip CRBN/DDB1/ligand
trim(PDB_DIR / "6H0F.pdb", "IKZF1_6H0F_ZF2", chains=["C"],
     notes="IKZF1 ZF2 in CRBN complex; PROTAC track; strip CRBN/DDB1 (keep chain C)")
# AF2: ZF1-4 DNA-binding domain (residues 162-320) — primary for DNA-binding pocket docking
trim(AF2_DIR / "IKZF1_AF2.pdb", "IKZF1_AF2_ZF14", chains=["A"],
     res_min=162, res_max=320, min_plddt=50, notes="ZF1-4 DBD, pLDDT≥50")

log.info("=== IRF8 ===")
# No IRF8 crystal structure in PDB — AF2 only
# AF2: DBD (residues 1-240)
trim(AF2_DIR / "IRF8_AF2.pdb", "IRF8_AF2_DBD", chains=["A"],
     res_min=1, res_max=240, min_plddt=50, notes="DBD, pLDDT≥50 — no PDB structure available")

log.info("=== BHLHE41 ===")
# AF2 only: bHLH domain (residues 1-200), strip low-confidence C-term
trim(AF2_DIR / "BHLHE41_AF2.pdb", "BHLHE41_AF2_bHLH", chains=["A"],
     res_min=1, res_max=200, min_plddt=50, notes="bHLH, pLDDT≥50")

log.info("=== RUNX1 ===")
# 4L0Z: RUNX1+ETS1+DNA (likely mouse, highly conserved with human; note in Methods)
# Chain A = RUNX1 Runt domain; strip ETS1 (chain B) and DNA (chains C/D)
trim(PDB_DIR / "4L0Z.pdb", "RUNX1_4L0Z_Runt", chains=["A"],
     notes="Runt domain, chain A; ETS1+DNA stripped; likely mouse — note in Methods")
# 1LJM: human RUNX1 Runt domain alone, 2.50 Å — ideal human structure
trim(PDB_DIR / "1LJM.pdb", "RUNX1_1LJM_Runt", chains=["A"],
     notes="Human Runt domain, 2.50 Å — preferred human structure")
# AF2: Runt domain (residues 51-188)
trim(AF2_DIR / "RUNX1_AF2.pdb", "RUNX1_AF2_Runt", chains=["A"],
     res_min=51, res_max=188, min_plddt=50, notes="Runt domain, pLDDT≥50")

log.info("=== RUNX2 ===")
# No RUNX2 PDB structure below 3.0 Å threshold — AF2 only
# (6VGD/6VGE exist at 4.2 Å but exceed plan's resolution cutoff)
# AF2: Runt domain (residues 80-220, similar to RUNX1)
trim(AF2_DIR / "RUNX2_AF2.pdb", "RUNX2_AF2_Runt", chains=["A"],
     res_min=80, res_max=220, min_plddt=50, notes="Runt domain, pLDDT≥50 — no PDB < 3.0 Å")

log.info("=== MAF ===")
# 4EOT: MafA homodimer (closest large-Maf PDB; no human c-Maf structure exists)
# Chain A = MafA bZIP; trim to bZIP core (residues 275-359 in MafA numbering)
trim(PDB_DIR / "4EOT.pdb", "MAF_4EOT_bZIP", chains=["A"],
     res_min=275, res_max=359, notes="MafA bZIP homolog (no human c-Maf PDB)")
# AF2: bZIP domain (residues 258-359, canonical bZIP window)
trim(AF2_DIR / "MAF_AF2.pdb", "MAF_AF2_bZIP", chains=["A"],
     res_min=258, res_max=359, min_plddt=50, notes="bZIP domain, pLDDT≥50")

log.info("=== SPI1 ===")
# 8EE9: human PU.1 ETS domain + DNA (1.22 Å, Cell Reports 2023, highest-res human structure)
# Chain F = ETS domain protein; other chains = DNA — strip DNA, keep chain F only
# ETS domain: UniProt P17947 feature ~165-270
trim(PDB_DIR / "8EE9.pdb", "SPI1_8EE9_ETS", chains=["F"],
     res_min=165, res_max=270,
     notes="ETS domain only; DNA chains stripped; chain F; 1.22 Å human structure")
# AF2: ETS domain window with pLDDT filter (SPI1 is well-modelled by AF2)
trim(AF2_DIR / "SPI1_AF2.pdb", "SPI1_AF2_ETS", chains=["A"],
     res_min=165, res_max=270, min_plddt=50,
     notes="ETS domain, pLDDT≥50; human P17947")

log.info("=== ACSL1 ===")
# No mammalian crystal structure in PDB for ACSL1 (P33121) — AF2 only
# AF2 human ACSL1: full-length AMP-binding domain; pLDDT filter removes disordered regions
# Active site is in the N-terminal large domain (~1-490); C-gate domain (~490-698)
trim(AF2_DIR / "ACSL1_AF2.pdb", "ACSL1_AF2_full", chains=["A"],
     min_plddt=50, notes="AF2 only — no mammalian PDB structure; full-length pLDDT≥50; human P33121")

log.info("=== PIK3CA ===")
# 4OVU: PIK3CA p110α + p85α-nSH2 + BEZ235 dual PI3K/mTOR inhibitor (2.93 Å)
# Chain A = p110α catalytic subunit; chain B = p85α-nSH2 regulatory fragment
# Kinase domain: ~726-1068 (canonical PI3K KD; ATP-binding site at K802)
trim(PDB_DIR / "4OVU.pdb", "PIK3CA_4OVU_KD", chains=["A"],
     res_min=726, res_max=1068,
     notes="Kinase domain only; p85-nSH2 (chain B) stripped; BEZ235 removed (HETATM); 2.93 Å")
# AF2 human PIK3CA: kinase domain window
trim(AF2_DIR / "PIK3CA_AF2.pdb", "PIK3CA_AF2_KD", chains=["A"],
     res_min=726, res_max=1068, min_plddt=50,
     notes="Kinase domain, pLDDT≥50; human P42336")

log.info("\n=== Trimmed files ===")
for f in sorted(TRIM_DIR.glob("*.pdb")):
    log.info(f"  {f.name}  {f.stat().st_size//1024} KB")
log.info("Done. Ready for 03_prepare_structures.py (PDBFixer).")
