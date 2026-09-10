#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_cistarget_db
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=128G
# ─────────────────────────────────────────────────────────────────────────────
# WP1.2 — Build a genome-wide custom cisTarget rankings database from
#         JASPAR 2026 CORE vertebrate motifs (decision D2 = full genome-wide).
#
# Produces:
#   data/resources/jaspar2026/db/
#     jaspar2026_hg38_10kbp_up_10kbp_down_full_tx.genes_vs_motifs.rankings.feather
#     jaspar2026_hg38_10kbp_up_10kbp_down_full_tx.genes_vs_motifs.scores.feather
#
# Requires (fill FIXMEs):
#   - JASPAR PFMs from s1_01 in data/resources/jaspar2026/pfms/
#   - cbust (Cluster-Buster) binary
#   - hg38.fa
#   - the gene-region BED used by the v10 "10kbp_up_10kbp_down_full_tx" DB
#   - the aertslab create_cisTarget_databases repo
#
# Submit: sbatch code/slurm/paper2/s1_02_build_cistarget_db.slurm
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
JDIR=$PROJDIR/data/resources/jaspar2026
DB=$JDIR/db
TOOLS=$PROJDIR/tools
mkdir -p "$PROJDIR/logs" "$DB" "$TOOLS" "$JDIR/cb"

echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME — start $(date)"

# ── 0. Tooling ──────────────────────────────────────────────────────────────
# FIXME: point CBUST at a Cluster-Buster binary (module load, conda, or build)
CBUST=${CBUST:-$TOOLS/cbust}
[ -x "$CBUST" ] || { echo "ERROR: cbust not found at $CBUST — set \$CBUST"; exit 1; }

if [ ! -d "$TOOLS/create_cisTarget_databases" ]; then
    git clone https://github.com/aertslab/create_cisTarget_databases "$TOOLS/create_cisTarget_databases"
fi
CTDB=$TOOLS/create_cisTarget_databases

# ── 1. Reference FASTA + gene regions ───────────────────────────────────────
HG38=${HG38_FA:-$PROJDIR/data/resources/genome/hg38.fa}
if [ ! -s "$HG38" ]; then
    mkdir -p "$(dirname "$HG38")"
    echo "Downloading hg38.fa ..."
    wget -q -O "$HG38.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip "$HG38.gz"
fi

# FIXME: the exact BED for "10kbp upstream + 10kbp downstream + full transcript",
# hg38, matching the v10 database the accepted paper used. See the
# create_cisTarget_databases wiki (regions/regulatory-region definitions).
REGION_BED=${REGION_BED:-$PROJDIR/data/resources/genome/hg38_10kbp_up_10kbp_down_full_tx.bed}
[ -s "$REGION_BED" ] || { echo "ERROR: set \$REGION_BED (see FIXME)"; exit 1; }

REGION_FA=$JDIR/hg38_10kbp_up_10kbp_down_full_tx.regions.fa
if [ ! -s "$REGION_FA" ]; then
    module load bedtools 2>/dev/null || true
    bedtools getfasta -fi "$HG38" -bed "$REGION_BED" -name -fo "$REGION_FA"
fi

# ── 2. PFM -> Cluster-Buster singletons ─────────────────────────────────────
if [ -z "$(ls -A "$JDIR/cb" 2>/dev/null)" ]; then
    python "$CTDB/convert_motifs_to_cb_format/convert_jaspar_to_cb_format.py" \
        --jaspar "$JDIR/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
        --output_dir "$JDIR/cb" \
    || python - << 'PY'
# minimal fallback: JASPAR PFM -> cbust .cb
from pathlib import Path
import re
src = Path("data/resources/jaspar2026/pfms"); dst = Path("data/resources/jaspar2026/cb")
for p in src.glob("*.pfm"):
    lines = p.read_text().splitlines()
    hdr = lines[0].lstrip(">").split()
    mid = hdr[0]
    rows = {}
    for ln in lines[1:]:
        m = re.match(r"\s*([ACGT])\s*\[\s*(.*?)\s*\]", ln) or re.match(r"\s*([ACGT])\s+(.*)", ln)
        if m:
            rows[m.group(1)] = [float(x) for x in m.group(2).split()]
    if len(rows) != 4:
        continue
    n = len(next(iter(rows.values())))
    out = [f">{mid}"]
    for i in range(n):
        out.append("\t".join(f"{rows[b][i]:.4f}" for b in "ACGT"))
    (dst / f"{mid}.cb").write_text("\n".join(out) + "\n")
print("fallback cb conversion done")
PY
fi
ls "$JDIR/cb"/*.cb | sed 's:.*/::; s:\.cb$::' > "$JDIR/motifs.lst"
echo "$(wc -l < "$JDIR/motifs.lst") motifs in .cb format"

# ── 3. Build the rankings database ─────────────────────────────────────────
PREFIX=$DB/jaspar2026_hg38_10kbp_up_10kbp_down_full_tx
python "$CTDB/create_cistarget_motif_databases.py" \
    --fasta "$REGION_FA" \
    --motifs_dir "$JDIR/cb" \
    --motifs "$JDIR/motifs.lst" \
    --output "$PREFIX" \
    --bgpadding 1000 \
    --cbust "$CBUST" \
    --threads 20

echo "Done $(date)"
ls -lh "$DB"
