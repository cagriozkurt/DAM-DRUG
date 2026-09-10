#!/bin/bash
# WP1.1 — Fetch JASPAR 2026 CORE vertebrate motifs + build a motif2tf table.
# Run on a login node (no SLURM). Decisions D1 (CORE vertebrate, non-redundant),
# D2 (genome-wide cisTarget DB).
#
#   export DAM_DRUG_DIR=/arf/scratch/mozkurt/DAM-DRUG
#   bash code/slurm/paper2/s1_01_fetch_jaspar2026.sh
#
# URLs verified 2026-09-10: JASPAR2026 CORE vertebrates non-redundant =
#   1,019 matrices. transfac (PO matrices, best for cbust), jaspar (names),
#   MEME (not needed here).
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
OUT=$PROJDIR/data/resources/jaspar2026
mkdir -p "$OUT/transfac" "$OUT/cb"

BASE=https://jaspar.elixir.no/download/data/2026/CORE
TRANSFAC=$OUT/JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt
JASPAR=$OUT/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt
META=$OUT/ultimate_metadata_table_CORE.tsv

[ -s "$TRANSFAC" ] || { echo "Downloading transfac PFMs..."; \
  curl -fsSL -o "$TRANSFAC" "$BASE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt"; }
[ -s "$JASPAR" ]   || { echo "Downloading jaspar PFMs..."; \
  curl -fsSL -o "$JASPAR" "$BASE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"; }
# metadata is best-effort (mirror can be flaky); motif2tf falls back to the
# name embedded in the transfac ID line.
curl -fsSL --max-time 60 -o "$META" \
  "https://mencius.uio.no/JASPAR/JASPAR_metadata/2026/ultimate_metadata_table_CORE.tsv" \
  2>/dev/null || echo "  (metadata table unavailable — using transfac ID names)"

echo "transfac matrices: $(grep -c '^AC ' "$TRANSFAC")"

# ── split transfac into one file per motif (for cbust conversion in s1_02) ───
python3 - "$OUT" << 'PY'
import sys, re
from pathlib import Path
out = Path(sys.argv[1]); tdir = out / "transfac"
blocks = re.split(r'(?m)^//\s*$', (out / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt").read_text())
n = 0
for b in blocks:
    b = b.strip("\n ")
    if not b:
        continue
    ac = re.search(r'(?m)^AC\s+(\S+)', b)
    if not ac:
        continue
    mid = ac.group(1)
    (tdir / f"{mid}.transfac").write_text(b.rstrip() + "\n//\n")
    n += 1
print(f"wrote {n} per-motif transfac files -> {tdir}")

# ── motif2tf table (pyscenic ctx --annotations_fname format) ────────────────
rows = ["#motif_id\tgene_name\tmotif_similarity_qvalue\torthologous_identity\tdescription\tannotation"]
seen = set()
for line in (out / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt").read_text().splitlines():
    if not line.startswith(">"):
        continue
    parts = line[1:].split("\t") if "\t" in line else line[1:].split()
    if len(parts) < 2:
        continue
    mid, name = parts[0], parts[1]
    # split heterodimers (TFA::TFB) and fusions (TFA-TFB) into individual TFs
    for tf in re.split(r"::|(?<=[a-z])-(?=[A-Z])", name):
        tf = tf.strip().upper()
        if tf and (mid, tf) not in seen:
            seen.add((mid, tf))
            rows.append(f"{mid}\t{tf}\t0.0\t1.0\tjaspar2026_CORE_vertebrates\tgene is directly annotated")
(out / "jaspar2026_motif2tf.tbl").write_text("\n".join(rows) + "\n")
print(f"wrote {len(rows)-1} motif->TF rows -> {out / 'jaspar2026_motif2tf.tbl'}")
PY

echo "Done. Next: sbatch code/slurm/paper2/s1_02_build_cistarget_db.slurm"
