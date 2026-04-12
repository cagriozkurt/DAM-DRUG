"""
DAM-DRUG Phase 3 — Parse fpocket output
=========================================
Reads *_info.txt files from all fpocket output directories and extracts
top-pocket metrics into a tidy CSV.

Targets: PPARG, IKZF1, IRF8, BHLHE41, RUNX1, RUNX2, MAF (7 original)
         + SPI1, ACSL1, PIK3CA (3 added in step 3.2b)
Handles both fpocket 3.x and 4.x output formats.
"""

import os
import re
import logging
import pandas as pd
import numpy as np
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
FPOCKET = PROJECT / "results/phase3/fpocket"
OUT_CSV = PROJECT / "results/phase3/pocket_summary.csv"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

# Map prepared file stem → TF and source
# Keys match stems in data/structures/prepared/ (trimmed stem + _prep)
STEM_MAP = {
    "PPARG_2PRG_LBD_prep":      ("PPARG",   "pdb", "2PRG"),
    "PPARG_1FM9_LBD_prep":      ("PPARG",   "pdb", "1FM9"),
    "PPARG_4EMA_LBD_prep":      ("PPARG",   "pdb", "4EMA"),
    "PPARG_AF2_LBD_prep":       ("PPARG",   "af2", "AF2_LBD"),
    "PPARG_AF2_full_prep":      ("PPARG",   "af2", "AF2_full"),
    "IKZF1_6H0F_ZF2_prep":      ("IKZF1",   "pdb", "6H0F"),  # ZF2 in CRBN complex; PROTAC track
    "IKZF1_AF2_ZF14_prep":      ("IKZF1",   "af2", "AF2_ZF14"),
    "IRF8_AF2_DBD_prep":        ("IRF8",    "af2", "AF2_DBD"),  # no PDB structure available
    "BHLHE41_AF2_bHLH_prep":    ("BHLHE41", "af2", "AF2_bHLH"),
    "RUNX1_4L0Z_Runt_prep":     ("RUNX1",   "pdb", "4L0Z"),   # likely mouse; note in Methods
    "RUNX1_1LJM_Runt_prep":     ("RUNX1",   "pdb", "1LJM"),   # human, 2.50 Å — preferred
    "RUNX1_AF2_Runt_prep":      ("RUNX1",   "af2", "AF2_Runt"),
    # RUNX2_AF2_Runt omitted: AF2 Q13580 is a 120-residue isoform; Runt domain truncated at
    # residue 120 — only 41 of ~170 Runt domain residues present; insufficient for fpocket.
    # RUNX2 therapeutic direction is ACTIVATION not inhibition (Sun 2023) — deprioritize docking.
    "MAF_4EOT_bZIP_prep":       ("MAF",     "pdb", "4EOT"),  # MafA homolog; no human c-Maf PDB
    # MAF_AF2_bZIP omitted: AF2 O15525 is only 162 residues; bZIP domain (258-359) absent.
    # ── Step 3.2b additions ────────────────────────────────────────────────────
    "SPI1_8EE9_ETS_prep":       ("SPI1",    "pdb", "8EE9"),  # ETS domain, chain F, 1.22 Å human
    "SPI1_AF2_ETS_prep":        ("SPI1",    "af2", "AF2_ETS"),
    "ACSL1_AF2_full_prep":      ("ACSL1",   "af2", "AF2_full"),  # no mammalian PDB
    "PIK3CA_4OVU_KD_prep":      ("PIK3CA",  "pdb", "4OVU"),  # kinase domain, 2.93 Å
    "PIK3CA_AF2_KD_prep":       ("PIK3CA",  "af2", "AF2_KD"),
}

def parse_info_file(info_path: Path) -> dict:
    """Parse fpocket *_info.txt and return metrics for pocket #1."""
    text = info_path.read_text()

    # Split into pocket blocks
    blocks = re.split(r"Pocket\s+\d+\s*:", text)
    if len(blocks) < 2:
        return {}

    pocket1 = blocks[1]  # first pocket (highest ranked)
    n_pockets = len(blocks) - 1

    def extract(pattern, text, default=np.nan):
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                return default
        return default

    return {
        "n_pockets_total":   n_pockets,
        "pocket_score":      extract(r"\n\s*Score\s*:\s*([\d.]+)", pocket1),
        "drug_score":        extract(r"Druggability Score\s*:\s*([\d.]+)", pocket1),
        "n_alpha_spheres":   extract(r"Number of Alpha Spheres\s*:\s*([\d.]+)", pocket1),
        "volume_A3":         extract(r"Volume\s*:\s*([\d.]+)", pocket1),
        "hydrophobicity":    extract(r"Hydrophobicity score\s*:\s*([\d.-]+)", pocket1),
        "polarity":          extract(r"Polarity score\s*:\s*([\d.-]+)", pocket1),
        "mean_local_hydro":  extract(r"Mean local hydrophobic density\s*:\s*([\d.]+)", pocket1),
        "charge_score":      extract(r"Charge score\s*:\s*([\d.-]+)", pocket1),
    }

rows = []

for out_dir in sorted(FPOCKET.iterdir()):
    if not out_dir.is_dir() or not out_dir.name.endswith("_out"):
        continue

    stem = out_dir.name[:-4]  # strip _out
    if stem not in STEM_MAP:
        log.warning(f"  Unknown structure: {stem}")
        continue

    tf, source, pdb_id = STEM_MAP[stem]

    # Find *_info.txt
    info_files = list(out_dir.glob("*_info.txt"))
    if not info_files:
        log.warning(f"  No info file in {out_dir.name}")
        continue

    info_path = info_files[0]
    metrics = parse_info_file(info_path)

    if not metrics:
        log.warning(f"  Could not parse {info_path.name}")
        continue

    row = {"tf": tf, "structure_id": stem, "source": source, "pdb_id": pdb_id}
    row.update(metrics)

    # For AF2: compute mean pLDDT of top-pocket atoms from pocket1_atm.pdb
    if source == "af2":
        pocket_atm = out_dir / "pockets" / "pocket1_atm.pdb"
        if pocket_atm.exists():
            bfactors = []
            for line in pocket_atm.read_text().splitlines():
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        bfactors.append(float(line[60:66].strip()))
                    except (ValueError, IndexError):
                        pass
            row["mean_plddt_pocket"] = np.mean(bfactors) if bfactors else np.nan
        else:
            row["mean_plddt_pocket"] = np.nan
    else:
        row["mean_plddt_pocket"] = np.nan

    log.info(f"  {stem:35s}  score={row.get('pocket_score', 'N/A'):.3f}  "
             f"drug={row.get('drug_score', 'N/A'):.3f}  "
             f"vol={row.get('volume_A3', 'N/A'):.0f}Å³  "
             f"n_pockets={row.get('n_pockets_total', '?')}")
    rows.append(row)

df = pd.DataFrame(rows)

col_order = ["tf", "structure_id", "source", "pdb_id", "n_pockets_total",
             "pocket_score", "drug_score", "n_alpha_spheres", "volume_A3",
             "hydrophobicity", "polarity", "mean_local_hydro", "charge_score",
             "mean_plddt_pocket"]
df = df[[c for c in col_order if c in df.columns]]
df.to_csv(OUT_CSV, index=False)

log.info(f"\nSaved {len(df)} rows → {OUT_CSV}")
log.info(df[["tf", "structure_id", "pocket_score", "drug_score", "volume_A3"]].to_string(index=False))
