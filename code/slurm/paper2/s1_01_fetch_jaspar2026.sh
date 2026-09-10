#!/bin/bash
# WP1.1 — Fetch JASPAR 2026 CORE vertebrate motifs + build a motif2tf table.
# Run on a login node (no SLURM). Decisions D1 (CORE vertebrate only), D2 (genome-wide DB).
#
#   export DAM_DRUG_DIR=/arf/scratch/mozkurt/DAM-DRUG
#   bash code/slurm/paper2/s1_01_fetch_jaspar2026.sh
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
OUT=$PROJDIR/data/resources/jaspar2026
mkdir -p "$OUT/pfms" "$OUT/cb"

# ── JASPAR 2026 CORE vertebrate, non-redundant PFMs ──────────────────────────
# FIXME: confirm the 2026 release path at https://jaspar.elixir.no/downloads/
JASPAR_PFMS_URL="https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
JASPAR_META_URL="https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meta.txt"

if [ ! -s "$OUT/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt" ]; then
    echo "Downloading JASPAR 2026 CORE vertebrate PFMs..."
    wget -q -O "$OUT/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt" "$JASPAR_PFMS_URL"
fi
wget -q -O "$OUT/JASPAR2026_CORE_vertebrates_meta.txt" "$JASPAR_META_URL" || \
    echo "WARN: meta download failed — motif2tf will fall back to matrix names"

# ── Split the combined jaspar file into one PFM per motif ────────────────────
python - "$OUT" << 'PY'
import sys, re
from pathlib import Path
out = Path(sys.argv[1]); pfms = out / "pfms"
txt = (out / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt").read_text().splitlines()
cur, name, n = [], None, 0
def flush():
    global cur, name, n
    if name and cur:
        (pfms / f"{name}.pfm").write_text("\n".join(cur) + "\n"); n += 1
    cur = []
for line in txt:
    if line.startswith(">"):
        flush()
        # >MA0004.1  Arnt
        parts = line[1:].split()
        name = parts[0]
        cur = [line]
    else:
        cur.append(line)
flush()
print(f"wrote {n} PFM files to {pfms}")

# motif2tf table: motif_id \t gene_name  (RcisTarget/pyscenic --annotations_fname format)
meta = out / "JASPAR2026_CORE_vertebrates_meta.txt"
rows = ["#motif_id\tgene_name\tmotif_similarity_qvalue\torthologous_identity\tdescription\tannotation"]
if meta.exists():
    for ln in meta.read_text().splitlines()[1:]:
        f = ln.split("\t")
        if len(f) >= 2 and f[0]:
            for tf in re.split(r"[:;()]+", f[1]):
                tf = tf.strip().upper()
                if tf:
                    rows.append(f"{f[0]}\t{tf}\t0.0\t1.0\tjaspar2026\tgene is directly annotated")
else:
    for line in txt:
        if line.startswith(">"):
            p = line[1:].split()
            if len(p) >= 2:
                rows.append(f"{p[0]}\t{p[1].upper()}\t0.0\t1.0\tjaspar2026\tgene is directly annotated")
(out / "jaspar2026_motif2tf.tbl").write_text("\n".join(rows) + "\n")
print(f"wrote {len(rows)-1} motif->tf rows")
PY

echo "Done. Next: sbatch code/slurm/paper2/s1_02_build_cistarget_db.slurm"
