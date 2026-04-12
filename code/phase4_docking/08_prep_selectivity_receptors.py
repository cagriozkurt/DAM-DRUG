"""
DAM-DRUG Phase 4 — Selectivity off-target receptor preparation
===============================================================
Downloads and prepares three CNS off-target receptors for selectivity docking
of the 6 validated MM-GBSA hits (CELECOXIB, TERIFLUNOMIDE, AMLEXANOX,
CAPMATINIB, LASMIDITAN, RISDIPLAM).

Off-targets:
  DRD2    — Dopamine D2 receptor        (PDB 6CM4, risperidone/8NU, 2.87 Å)
  HTR2A   — 5-HT2A serotonin receptor   (PDB 6A94, zotepine/ZOT, 2.9 Å)
  HERG    — hERG potassium channel      (PDB 7CN1, K+ cavity marker, 3.7 Å cryo-EM)

Run on TRUBA login node (no GPU needed):
    apptainer exec containers/scenic.sif python code/phase4_docking/08_prep_selectivity_receptors.py

Requires: obabel (available in container)
"""

import logging
import os
import subprocess
import urllib.request
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT   = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
PDB_DIR   = PROJECT / "data/structures/pdb"
RECV_DIR  = PROJECT / "data/docking/receptors"
CONF_DIR  = PROJECT / "data/docking/configs"
OBABEL    = os.environ.get("OBABEL", "obabel")

for d in (PDB_DIR, RECV_DIR, CONF_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Off-target receptor definitions
# Fields:
#   chains       — PDB chain IDs to extract as receptor
#   lig_hint     — HETATM residue name (CCD code) for binding-site detection
#   min_lig_atoms— minimum heavy atoms for a HETATM group to be used as box anchor
#                  (set to 1 for ion markers like K+)
#   box_z_offset — shift applied to centroid Z after detection (use for ion markers
#                  that sit above the true drug-binding cavity)
#   box          — Vina grid box side length (Å)
OFFTARGETS = [
    {
        "name":         "DRD2",
        "pdb_id":       "6CM4",
        "chains":       {"A"},          # DRD2 receptor (chain A)
        "lig_hint":     "8NU",          # risperidone CCD code
        "min_lig_atoms": 8,
        "box_z_offset":  0.0,
        "box":          22,
    },
    {
        "name":         "HTR2A",
        "pdb_id":       "6A94",
        "chains":       {"A"},          # 5-HT2A receptor (chain A); b562 fusion is chain B
        "lig_hint":     "ZOT",          # zotepine (atypical antipsychotic, antagonist)
        "min_lig_atoms": 8,
        "box_z_offset":  0.0,
        "box":          22,
    },
    {
        "name":         "HERG",
        "pdb_id":       "7CN1",
        "chains":       {"A", "B", "C", "D"},  # full homotetramer (drug site in central cavity)
        "lig_hint":     "K",            # K+ ions in selectivity filter mark the cavity axis
        "min_lig_atoms": 1,             # K+ is a single atom
        "box_z_offset": -5.0,           # drug cavity is ~5 Å below the selectivity filter
        "box":          26,
    },
]


def download_pdb(pdb_id: str, out_path: Path) -> bool:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log.info(f"  Downloading {pdb_id} from RCSB...")
    try:
        urllib.request.urlretrieve(url, out_path)
        log.info(f"  Downloaded: {out_path.name}  ({out_path.stat().st_size // 1024} KB)")
        return True
    except Exception as e:
        log.error(f"  Download failed: {e}")
        return False


def find_ligand_centroid(pdb_text: str, lig_hint: str, chains: set,
                         min_atoms: int = 8) -> tuple | None:
    """
    Find centroid of co-crystallized ligand in given chains.
    Tries lig_hint first; falls back to any HETATM with ≥min_atoms heavy atoms
    (excludes common buffer/solvent residues when min_atoms > 1).
    """
    SKIP = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "MPD",
            "PEG", "TRS", "MES", "HEX", "IMD", "ACT", "FMT"}

    # Group HETATM atoms by (chain, resseq, resname)
    groups: dict = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        chain   = line[21]
        resseq  = line[22:26].strip()
        if chain not in chains:
            continue
        # Skip common buffers/solvents only when looking for drug-sized molecules
        if min_atoms > 1 and resname in SKIP:
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except (ValueError, IndexError):
            continue
        key = (chain, resseq, resname)
        groups.setdefault(key, []).append((x, y, z))

    if not groups:
        return None

    # Prefer ligand matching hint
    hint = lig_hint[:3].upper()
    candidates = [(k, v) for k, v in groups.items() if k[2].upper().startswith(hint)]
    if not candidates:
        # Fall back to largest group meeting minimum atom count
        candidates = sorted(groups.items(), key=lambda kv: -len(kv[1]))
        candidates = [(k, v) for k, v in candidates if len(v) >= min_atoms]

    if not candidates:
        return None

    # Use all matching instances (e.g. 4 K+ ions across tetramer) for centroid
    all_coords = [c for _, v in candidates for c in v]
    n = len(all_coords)
    cx = sum(c[0] for c in all_coords) / n
    cy = sum(c[1] for c in all_coords) / n
    cz = sum(c[2] for c in all_coords) / n
    rname = candidates[0][0][2]
    log.info(f"  Ligand centroid ({rname}, {len(candidates)} instance(s), {n} atoms): "
             f"({cx:.2f}, {cy:.2f}, {cz:.2f})")
    return (cx, cy, cz)


def extract_receptor(pdb_text: str, chains: set) -> str:
    """Keep only ATOM records from specified chains + TER/END. Strip HETATM."""
    out_lines = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[21] in chains:
            out_lines.append(line)
        elif line.startswith(("TER", "END")):
            out_lines.append(line)
    return "\n".join(out_lines) + "\n"


def prepare_pdbqt(pdb_in: Path, pdbqt_out: Path) -> bool:
    cmd = [OBABEL, str(pdb_in), "-O", str(pdbqt_out),
           "-xr", "--partialcharge", "gasteiger", "-h"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not pdbqt_out.exists() or pdbqt_out.stat().st_size < 100:
        log.error(f"  obabel failed: {result.stderr.strip()[:300]}")
        return False
    log.info(f"  PDBQT: {pdbqt_out.name}  ({pdbqt_out.stat().st_size // 1024} KB)")
    return True


def write_vina_conf(conf_path: Path, pdbqt_path: Path,
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
    log.info(f"  Config: {conf_path.name}  (box={box}Å, centre=({cx:.1f},{cy:.1f},{cz:.1f}))")


def main():
    ok_count = 0
    for t in OFFTARGETS:
        name    = t["name"]
        pdb_id  = t["pdb_id"]
        chains  = t["chains"]
        log.info(f"\n=== {name} ({pdb_id}) ===")

        pdb_raw  = PDB_DIR  / f"{pdb_id}.pdb"
        pdb_rec  = PDB_DIR  / f"{name}_receptor.pdb"
        pdbqt    = RECV_DIR / f"{name}_prep.pdbqt"
        conf     = CONF_DIR / f"{name}.conf"

        # Download
        if not pdb_raw.exists():
            if not download_pdb(pdb_id, pdb_raw):
                log.error(f"  Skipping {name}")
                continue

        pdb_text = pdb_raw.read_text()

        # Find binding site centre from co-crystallized ligand/marker
        centre = find_ligand_centroid(pdb_text, t["lig_hint"], chains,
                                      min_atoms=t.get("min_lig_atoms", 8))
        if centre is None:
            log.warning(f"  Ligand/marker not found — check chain IDs and lig_hint")
            continue

        cx, cy, cz = centre
        z_off = t.get("box_z_offset", 0.0)
        if z_off:
            cz += z_off
            log.info(f"  Z offset {z_off:+.1f} Å applied → ({cx:.2f}, {cy:.2f}, {cz:.2f})")

        # Extract receptor chains (ATOM only, no HETATM)
        rec_text = extract_receptor(pdb_text, chains)
        n_atoms  = rec_text.count("\nATOM")
        if n_atoms < 100:
            log.error(f"  Too few ATOM records ({n_atoms}) for chains {chains} — check chain IDs")
            continue
        pdb_rec.write_text(rec_text)
        log.info(f"  Receptor PDB: {n_atoms} atoms, chains {chains}")

        # Convert to PDBQT
        if pdbqt.exists():
            log.info(f"  PDBQT exists, skipping obabel")
        else:
            if not prepare_pdbqt(pdb_rec, pdbqt):
                continue

        # Write Vina config
        write_vina_conf(conf, pdbqt, cx, cy, cz, t["box"])
        ok_count += 1

    log.info(f"\n{'='*60}")
    log.info(f"Prepared {ok_count}/{len(OFFTARGETS)} off-target receptors")
    log.info(f"Receptors: {RECV_DIR}/DRD2_prep.pdbqt, HTR2A_prep.pdbqt, HERG_prep.pdbqt")
    log.info(f"Configs:   {CONF_DIR}/DRD2.conf, HTR2A.conf, HERG.conf")
    log.info(f"Next: sbatch code/slurm/31_run_selectivity_docking.slurm")


if __name__ == "__main__":
    main()
