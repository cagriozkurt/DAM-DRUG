#!/bin/bash
# Account via $SBATCH_ACCOUNT
#SBATCH --job-name=p2_ctx_jaspar
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=180G
# ─────────────────────────────────────────────────────────────────────────────
# WP1.3 — pyscenic ctx + aucell using the JASPAR 2026 custom database.
# Re-prunes the SAME 5-seed GRNBoost2 consensus edge table used for the
# accepted paper, so any change in retained regulons is attributable to the
# motif database alone.
#
# Inputs:
#   results/phase2/GRN/adj_matrix_aggregated.tsv        (5-seed consensus, unchanged)
#   results/phase2/GRN/microglia_raw.loom               (100K-cell expression, unchanged)
#   data/resources/jaspar2026/db/*.rankings.feather     (from s1_02)
#   data/resources/jaspar2026/jaspar2026_motif2tf.tbl   (from s1_01)
# Outputs -> results/phase7/grn_jaspar2026/
#   regulons_jaspar2026.csv
#   scenic_auc_jaspar2026.loom
#   regulon_auc_by_state_jaspar2026.csv
#   TF_regulon_summary_jaspar2026.csv
# ─────────────────────────────────────────────────────────────────────────────
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
export DAM_DRUG_DIR=$PROJDIR
SIF=$PROJDIR/containers/scenic.sif                       # D3: :latest kept; pin @sha256 at deposit
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest

GRN=$PROJDIR/results/phase2/GRN
OUT=$PROJDIR/results/phase7/grn_jaspar2026
mkdir -p "$PROJDIR/logs" "$OUT"

ADJ=$GRN/adj_matrix_aggregated.tsv
LOOM=$GRN/microglia_raw.loom
RANKINGS=$(ls "$PROJDIR"/data/resources/jaspar2026/db/*.genes_vs_motifs.rankings.feather | head -1)
MOTIF_TBL=$PROJDIR/data/resources/jaspar2026/jaspar2026_motif2tf.tbl
REG=$OUT/regulons_jaspar2026.csv
AUC=$OUT/scenic_auc_jaspar2026.loom

for f in "$ADJ" "$LOOM" "$RANKINGS" "$MOTIF_TBL"; do
    [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done
echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME — start $(date)"
echo "rankings: $RANKINGS"

BIND="--bind $PROJDIR:$PROJDIR --env LD_LIBRARY_PATH=/opt/conda/envs/scenic/lib"

if [ ! -s "$REG" ]; then
    echo "=== pyscenic ctx (JASPAR 2026): $(date) ==="
    apptainer exec $BIND "$SIF" pyscenic ctx \
        "$ADJ" "$RANKINGS" \
        --annotations_fname "$MOTIF_TBL" \
        --expression_mtx_fname "$LOOM" \
        --output "$REG" \
        --num_workers 20 --mode dask_multiprocessing --mask_dropouts
fi

if [ ! -s "$AUC" ]; then
    echo "=== pyscenic aucell: $(date) ==="
    apptainer exec $BIND "$SIF" pyscenic aucell \
        "$LOOM" "$REG" --output "$AUC" --num_workers 20
fi

echo "=== summarise regulon activity by substate: $(date) ==="
apptainer exec $BIND "$SIF" conda run -n scenic python - << 'PY'
import os, loompy, pandas as pd
from pathlib import Path
OUT = Path(os.environ["DAM_DRUG_DIR"]) / "results/phase7/grn_jaspar2026"
with loompy.connect(str(OUT / "scenic_auc_jaspar2026.loom"), mode="r", validate=False) as ds:
    a = ds.ca["RegulonsAUC"]; regs = list(a.dtype.names)
    df = pd.DataFrame({r: a[r] for r in regs}, index=ds.ca["CellID"])
    df["state"] = list(ds.ca["state"])
sm = df.groupby("state").mean()
sm.to_csv(OUT / "regulon_auc_by_state_jaspar2026.csv")
rows = []
for tf in ["SPI1","RUNX1","IRF8","PPARG","CEBPB","IKZF1","IKZF2","IKZF3","RELB","BHLHE40","BHLHE41"]:
    for reg in [c for c in sm.columns if c.split("(")[0] == tf]:
        s = sm[reg]
        rows.append({"TF": tf, "regulon": reg, "peak_state": s.idxmax(),
                     "peak_AUC": float(s.max()), **{f"AUC_{k}": float(v) for k, v in s.items()}})
    if not any(c.split("(")[0] == tf for c in sm.columns):
        rows.append({"TF": tf, "regulon": "NOT_RETAINED"})
pd.DataFrame(rows).to_csv(OUT / "TF_regulon_summary_jaspar2026.csv", index=False)
print(f"{sm.shape[1]} regulons retained; TF summary written")
PY

echo "Done $(date). Next: sbatch code/slurm/paper2/s1_04_benchmark_nulls.slurm"
