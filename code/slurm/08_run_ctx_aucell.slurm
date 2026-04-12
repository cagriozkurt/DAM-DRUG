#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=ctx_aucell
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/ctx_aucell-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/ctx_aucell-%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=180G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — pyscenic ctx + aucell on aggregated adjacency matrix
#
# Input:  results/phase2/GRN/adj_matrix_aggregated.tsv  (5-seed consensus)
# Output: results/phase2/GRN/regulons_aggregated.csv
#         results/phase2/GRN/scenic_auc_aggregated.loom
#         results/phase2/GRN/regulon_auc_by_state_aggregated.csv
#         results/phase2/GRN/TF_regulon_summary_aggregated.csv
#
# Submit: sbatch code/slurm/26_run_ctx_aucell.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-/arf/scratch/mozkurt/DAM-DRUG}
export DAM_DRUG_DIR=$PROJDIR
SIF=$HOME/containers/dam-drug-scanpy.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy:latest
GRN=$PROJDIR/results/phase2/GRN

ADJ=$GRN/adj_matrix_aggregated.tsv
LOOM=$GRN/microglia_raw.loom
RANKINGS=$PROJDIR/data/resources/motifs/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
MOTIF_TBL=$PROJDIR/data/resources/motifs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

REG=$GRN/regulons_aggregated.csv
AUC_LOOM=$GRN/scenic_auc_aggregated.loom

echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"
mkdir -p $PROJDIR/logs

[ -f "$ADJ" ]   || { echo "ERROR: $ADJ not found"; exit 1; }
[ -f "$LOOM" ]  || { echo "ERROR: $LOOM not found"; exit 1; }

# ── Step 1: pyscenic ctx ──────────────────────────────────────────────────────
if [ -f "$REG" ] && [ -s "$REG" ]; then
    echo "Regulons already exist ($REG) — skipping ctx"
else
    echo "=== pyscenic ctx: $(date) ==="
    apptainer exec \
      --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
      "$SIF" \
      conda run -n scenic \
        pyscenic ctx \
          "$ADJ" \
          "$RANKINGS" \
          --annotations_fname "$MOTIF_TBL" \
          --expression_mtx_fname "$LOOM" \
          --output "$REG" \
          --num_workers 20 \
          --mode dask_multiprocessing \
          --mask_dropouts
    echo "ctx done: $(date)"
fi

# ── Step 2: pyscenic aucell ───────────────────────────────────────────────────
if [ -f "$AUC_LOOM" ] && [ -s "$AUC_LOOM" ]; then
    echo "AUC loom already exists ($AUC_LOOM) — skipping aucell"
else
    echo "=== pyscenic aucell: $(date) ==="
    apptainer exec \
      --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
      "$SIF" \
      conda run -n scenic \
        pyscenic aucell \
          "$LOOM" \
          "$REG" \
          --output "$AUC_LOOM" \
          --num_workers 20
    echo "aucell done: $(date)"
fi

# ── Step 3: summarize regulon activity by state ───────────────────────────────
echo "=== Summarizing regulon activity: $(date) ==="
apptainer exec \
  --bind /arf/scratch/mozkurt:/arf/scratch/mozkurt \
  "$SIF" \
  conda run -n scenic python - << 'PYEOF'
import os, loompy, pandas as pd
from pathlib import Path

GRN = Path(os.environ.get("DAM_DRUG_DIR", "/arf/scratch/mozkurt/DAM-DRUG")) / "results/phase2/GRN"
AUC_LOOM = GRN / "scenic_auc_aggregated.loom"
TF_TARGETS = ["SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1",
              "RELB", "BHLHE40", "BHLHE41"]

with loompy.connect(str(AUC_LOOM), mode="r", validate=False) as ds:
    auc_raw = ds.ca["RegulonsAUC"]
    regulon_names = list(auc_raw.dtype.names)
    auc_mtx = pd.DataFrame(
        {name: auc_raw[name] for name in regulon_names},
        index=ds.ca["CellID"],
    )
    cell_states = pd.Series(ds.ca["state"], index=ds.ca["CellID"])

auc_mtx["state"] = cell_states
state_mean = auc_mtx.groupby("state").mean()
state_mean.to_csv(GRN / "regulon_auc_by_state_aggregated.csv")
print(f"State AUC matrix: {state_mean.shape}")

summary_rows = []
for tf in TF_TARGETS:
    matches = [c for c in state_mean.columns if c.startswith(f"{tf}(")]
    if matches:
        for reg in matches:
            aucs = state_mean[reg]
            peak = aucs.idxmax()
            row = {"TF": tf, "regulon": reg,
                   "peak_state": peak, "peak_AUC": float(aucs.max())}
            for state in state_mean.index:
                row[f"AUC_{state}"] = float(aucs[state])
            summary_rows.append(row)
            print(f"  {reg:<40} peak={peak} ({aucs.max():.4f})")
    else:
        print(f"  {tf:<40} NO REGULON FOUND")
        summary_rows.append({"TF": tf, "regulon": "NOT_FOUND"})

pd.DataFrame(summary_rows).to_csv(
    GRN / "TF_regulon_summary_aggregated.csv", index=False
)
print(f"\nSaved: {GRN / 'TF_regulon_summary_aggregated.csv'}")
PYEOF

echo ""
echo "All done: $(date)"
echo "Outputs:"
echo "  $REG"
echo "  $AUC_LOOM"
echo "  $GRN/regulon_auc_by_state_aggregated.csv"
echo "  $GRN/TF_regulon_summary_aggregated.csv"
