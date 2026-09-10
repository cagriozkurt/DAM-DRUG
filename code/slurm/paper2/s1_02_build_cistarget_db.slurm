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
# Inputs (all fetched/built by this script or s1_01 — no manual FIXMEs):
#   - per-motif transfac files from s1_01 in data/resources/jaspar2026/transfac/
#   - cbust: aertslab precompiled binary, with a source-build fallback
#   - hg38.fa: UCSC
#   - region BED: the CANONICAL aertslab v10 region file (verified URL below)
#   - aertslab create_cisTarget_databases repo (cloned)
#
# TRUBA env (verified 2026-09-10): scenic.sif has pyscenic + python3.10 in the
# `scenic` conda env but NOT cbust / bedtools / samtools; no host bedtools
# module. => cbust is downloaded/built here; FASTA extraction is pure-python.
#
# Submit: sbatch code/slurm/paper2/s1_02_build_cistarget_db.slurm
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif
JDIR=$PROJDIR/data/resources/jaspar2026
DB=$JDIR/db
TOOLS=$PROJDIR/tools
mkdir -p "$PROJDIR/logs" "$DB" "$TOOLS" "$JDIR/cb"
PY_SCENIC="apptainer exec --bind $PROJDIR:$PROJDIR $SIF conda run -n scenic python"

echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME — start $(date)"

# ── 0. Cluster-Buster ─────────────────────────────────────────────────────
CBUST=${CBUST:-$TOOLS/cbust}
if [ ! -x "$CBUST" ]; then
    echo "Fetching cbust (aertslab precompiled)..."
    wget -q -O "$CBUST" https://resources.aertslab.org/cistarget/programs/cbust && chmod +x "$CBUST" || true
fi
if ! "$CBUST" --help >/dev/null 2>&1 && ! "$CBUST" 2>&1 | grep -qi cluster; then
    echo "Precompiled cbust not runnable — building from source (ghuls fork)..."
    module load comp/gcc/12.3.0 2>/dev/null || module load gcc 2>/dev/null || true
    rm -rf "$TOOLS/cluster-buster"
    git clone https://github.com/ghuls/cluster-buster "$TOOLS/cluster-buster"
    make -C "$TOOLS/cluster-buster" cbust
    CBUST="$TOOLS/cluster-buster/cbust"
fi
echo "cbust: $CBUST"

if [ ! -d "$TOOLS/create_cisTarget_databases" ]; then
    git clone https://github.com/aertslab/create_cisTarget_databases "$TOOLS/create_cisTarget_databases"
fi
CTDB=$TOOLS/create_cisTarget_databases

# ── 1. Reference genome + canonical v10 region BED ────────────────────────
HG38=${HG38_FA:-$PROJDIR/data/resources/genome/hg38.fa}
if [ ! -s "$HG38" ]; then
    mkdir -p "$(dirname "$HG38")"
    echo "Downloading hg38.fa (~1 GB)..."
    wget -q -O "$HG38.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip "$HG38.gz"
fi
GENOME_DIR=$(dirname "$HG38")

# Canonical aertslab v10 gene-based region set: each gene split into
# named segments (GENE#1, GENE#2, ...), 10 kbp up + 10 kbp down + full
# transcript. This is the exact region definition behind
# hg38_10kbp_up_10kbp_down_full_tx_v10_clust.  (verified 2026-09-10: 200 OK,
# 92,636 lines; format = chrom start end GENE#N score strand)
REGION_BED=$GENOME_DIR/hg38-limited-upstream10000-tss-downstream10000-full-transcript.bed
[ -s "$REGION_BED" ] || wget -q -O "$REGION_BED" \
  https://resources.aertslab.org/cistarget/regions/hg38-limited-upstream10000-tss-downstream10000-full-transcript.bed

# ── FASTA extraction (pure-python; no bedtools/samtools on TRUBA) ──────────
REGION_FA=$JDIR/hg38_10kbp_up_10kbp_down_full_tx.regions.fa
if [ ! -s "$REGION_FA" ]; then
  HG38="$HG38" BED="$REGION_BED" FA="$REGION_FA" $PY_SCENIC - << 'PY'
import os
hg38, bed, out = os.environ["HG38"], os.environ["BED"], os.environ["FA"]

# 1. index the FASTA: chrom -> (data_offset, line_bases, line_width)
idx, chrom, off = {}, None, 0
with open(hg38, "rb") as f:
    while True:
        pos = f.tell()
        line = f.readline()
        if not line:
            break
        if line.startswith(b">"):
            chrom = line[1:].split()[0].decode()
            data_start = f.tell()
            nxt = f.readline()
            lb = len(nxt.rstrip(b"\n")); lw = len(nxt)
            idx[chrom] = [data_start, lb, lw]
            f.seek(data_start)

def fetch(c, s, e):
    if c not in idx:
        return None
    ds, lb, lw = idx[c]
    def byte(p):  # 0-based base position -> byte offset
        return ds + (p // lb) * lw + (p % lb)
    with open(hg38, "rb") as f:
        f.seek(byte(s))
        raw = f.read(byte(e) - byte(s))
    return raw.replace(b"\n", b"").upper().decode()

n = 0
with open(bed) as b, open(out, "w") as o:
    for ln in b:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        c, s, e, name = p[0], int(p[1]), int(p[2]), p[3]
        seq = fetch(c, s, e)
        if not seq or set(seq) <= {"N"}:
            continue
        o.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            o.write(seq[i:i + 60] + "\n")
        n += 1
print(f"{n} region sequences -> {out}")
PY
fi
echo "region FASTA: $(grep -c '^>' "$REGION_FA") sequences"

# ── 2. transfac PFMs -> Cluster-Buster singletons ──────────────────────────
# transfac per-motif files written by s1_01 ($JDIR/transfac/*.transfac).
if [ -z "$(ls -A "$JDIR/cb" 2>/dev/null)" ]; then
    JDIR="$JDIR" $PY_SCENIC - << 'PY'
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
# create_cistarget_motif_databases.py needs numpy/pandas/pyarrow/numba. The
# `scenic` env in scenic.sif has pyarrow but numba is unverified — if the run
# fails on ImportError, create the aertslab env once:
#   apptainer exec $SIF conda env create -f $CTDB/environment.yml -n ctdb
# and switch the two `conda run -n scenic` calls below to `-n ctdb`.
PREFIX=$DB/jaspar2026_hg38_10kbp_up_10kbp_down_full_tx
apptainer exec --bind "$PROJDIR:$PROJDIR" "$SIF" conda run -n scenic \
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
[ -s "$PREFIX.genes_vs_motifs.rankings.feather" ] || \
    { echo "ERROR: rankings feather not produced"; exit 1; }
