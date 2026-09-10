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
# Cluster-Buster — precompiled binary from the aertslab resources server.
CBUST=${CBUST:-$TOOLS/cbust}
if [ ! -x "$CBUST" ]; then
    echo "Downloading cbust..."
    wget -q -O "$CBUST" https://resources.aertslab.org/cistarget/programs/cbust
    chmod +x "$CBUST"
fi
"$CBUST" 2>&1 | head -1 || { echo "ERROR: cbust not runnable (glibc?)"; exit 1; }

if [ ! -d "$TOOLS/create_cisTarget_databases" ]; then
    git clone https://github.com/aertslab/create_cisTarget_databases "$TOOLS/create_cisTarget_databases"
fi
CTDB=$TOOLS/create_cisTarget_databases

module load bedtools 2>/dev/null || true

# ── 1. Reference genome + gene regulatory regions ──────────────────────────
HG38=${HG38_FA:-$PROJDIR/data/resources/genome/hg38.fa}
if [ ! -s "$HG38" ]; then
    mkdir -p "$(dirname "$HG38")"
    echo "Downloading hg38.fa (~1 GB)..."
    wget -q -O "$HG38.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip "$HG38.gz"
fi
[ -s "$HG38.fai" ] || samtools faidx "$HG38"

# Region set = each gene's span extended 10 kbp upstream + 10 kbp downstream,
# merged per gene symbol. Approximates the aertslab v10
# "10kbp_up_10kbp_down_full_tx" gene-based region set (v10 used RefSeq r80;
# here UCSC refGene — minor boundary differences, documented as a deviation).
GENOME_DIR=$(dirname "$HG38")
REFGENE=$GENOME_DIR/hg38.refGene.txt.gz
[ -s "$REFGENE" ] || wget -q -O "$REFGENE" \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

REGION_BED=$GENOME_DIR/hg38_10kbp_up_10kbp_down_full_tx.bed
if [ ! -s "$REGION_BED" ]; then
    HG38FAI="$HG38.fai" REFGENE="$REFGENE" OUT="$REGION_BED" python3 - << 'PY'
import gzip, os
from collections import defaultdict
chrom_len = {}
for ln in open(os.environ["HG38FAI"]):
    f = ln.split(); chrom_len[f[0]] = int(f[1])
PAD = 10000
spans = defaultdict(lambda: [None, None, None])   # gene -> [chrom, start, end]
with gzip.open(os.environ["REFGENE"], "rt") as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        chrom, txStart, txEnd, name2 = f[2], int(f[4]), int(f[5]), f[12]
        if chrom not in chrom_len or "_" in chrom:
            continue
        s = spans[name2]
        lo = max(0, txStart - PAD); hi = min(chrom_len[chrom], txEnd + PAD)
        if s[0] is None:
            spans[name2] = [chrom, lo, hi]
        elif s[0] == chrom:
            s[1] = min(s[1], lo); s[2] = max(s[2], hi)
with open(os.environ["OUT"], "w") as o:
    for g, (c, a, b) in sorted(spans.items()):
        if c:
            o.write(f"{c}\t{a}\t{b}\t{g}\n")
print(f"{len(spans)} gene regions -> {os.environ['OUT']}")
PY
fi

REGION_FA=$JDIR/hg38_10kbp_up_10kbp_down_full_tx.regions.fa
[ -s "$REGION_FA" ] || bedtools getfasta -fi "$HG38" -bed "$REGION_BED" -nameOnly -fo "$REGION_FA"
echo "region FASTA: $(grep -c '^>' "$REGION_FA") sequences"

# ── 2. transfac PFMs -> Cluster-Buster singletons ──────────────────────────
# transfac per-motif files written by s1_01 ($JDIR/transfac/*.transfac).
if [ -z "$(ls -A "$JDIR/cb" 2>/dev/null)" ]; then
    JDIR="$JDIR" python3 - << 'PY'
import os, re
from pathlib import Path
jd = Path(os.environ["JDIR"])
src, dst = jd / "transfac", jd / "cb"
dst.mkdir(exist_ok=True)
n = 0
for p in sorted(src.glob("*.transfac")):
    mid = None
    counts = []
    for ln in p.read_text().splitlines():
        if ln.startswith("AC "):
            mid = ln.split()[1]
        elif re.match(r"^\d+\s", ln):                       # PO matrix row: "01  A C G T"
            vals = ln.split()[1:5]
            if len(vals) == 4:
                counts.append([float(x) for x in vals])
    if mid and len(counts) >= 4:
        body = "\n".join("\t".join(f"{c:.4f}" for c in row) for row in counts)
        (dst / f"{mid}.cb").write_text(f">{mid}\n{body}\n")
        n += 1
print(f"{n} transfac -> Cluster-Buster .cb singletons")
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
