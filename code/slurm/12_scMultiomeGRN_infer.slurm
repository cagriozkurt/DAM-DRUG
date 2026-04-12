#!/bin/bash
# Account is read from $SBATCH_ACCOUNT — set it in your shell profile:
#   export SBATCH_ACCOUNT=your_account_name
#SBATCH --job-name=scmgrn_infer
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=akya-cuda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --array=0-5
# ─────────────────────────────────────────────────────────────────────────────
# Phase 2.8 — scMultiomeGRN inference per microglial state (SLURM array)
#
# Array index → state:
#   0=Homeostatic  1=DAM  2=DAM-IRM  3=IRM  4=LateAD-DAM  5=LAM
#
# Submit: sbatch code/slurm/12_scMultiomeGRN_infer.slurm
# One state test: sbatch --array=1 code/slurm/12_scMultiomeGRN_infer.slurm
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
TOOLDIR=$PROJDIR/tools/scMultiomeGRN/extracted/ScmultiomeGRN-main
DATADIR=$PROJDIR/tools/scMultiomeGRN/data_resource
INPUTDIR=$PROJDIR/results/phase2/scMultiomeGRN/input
OUTDIR=$PROJDIR/results/phase2/scMultiomeGRN/results

STATES=("Homeostatic" "DAM" "DAM-IRM" "IRM" "LateAD-DAM" "LAM")
STATE=${STATES[$SLURM_ARRAY_TASK_ID]}
STATE_TAG=$(echo "$STATE" | tr '-' '_')

mkdir -p "$PROJDIR/logs" "$OUTDIR"

echo "=== scMultiomeGRN: $STATE (task $SLURM_ARRAY_TASK_ID) ==="
echo "Job: $SLURM_JOB_ID on $SLURMD_NODENAME"
echo "Start: $(date)"

SIF=$PROJDIR/containers/scmultiomegrn.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scmultiomegrn:latest

echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'not detected yet')"

# ── Check input files ─────────────────────────────────────────────────────────
CELL_DIR=$INPUTDIR/$STATE_TAG
if [ ! -f "$CELL_DIR/scrna/matrix.mtx" ] || [ ! -f "$CELL_DIR/atac/matrix.mtx" ]; then
    echo "ERROR: Input MTX files missing for $STATE at $CELL_DIR"
    echo "Run: sbatch code/slurm/obsolete_scMultiomeGRN_env.slurm"
    exit 1
fi

# ── Write inference script to /tmp ────────────────────────────────────────────
SCRIPT=/tmp/scmgrn_${STATE_TAG}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
cat > "$SCRIPT" << PYEOF
import os, sys

# scMultiomeGRN source on path
sys.path.insert(0, "$TOOLDIR")
os.chdir("$TOOLDIR")

from src.utils import (
    get_TFBS_from_ATAC, get_TFBS_from_promoter,
    build_Graph, build_cross_graph,
    extract_atac_feature, extract_scrna_feature2,
)
from src.trainer import Trainer
import torch
import pandas as pd, numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────────
cell_type   = "$STATE_TAG"
cell_dir    = "$CELL_DIR"
atac_file   = os.path.join(cell_dir, "atac",  "matrix.mtx")
scrna_file  = os.path.join(cell_dir, "scrna", "matrix.mtx")
save_dir    = os.path.join("$OUTDIR", cell_type, "graph")
out_dir     = os.path.join("$OUTDIR", cell_type)

motif_dir     = "$DATADIR/HOCOMOCOv11"
genome_file   = "$DATADIR/download_genome/hg38"
promoter_file = "$TOOLDIR/data_resource/gencode.v38.ProteinCoding_gene_promoter.txt"

os.makedirs(save_dir, exist_ok=True)
os.makedirs(out_dir,  exist_ok=True)

# ── Step 1: TFBS from ATAC peaks ──────────────────────────────────────────────
print(f"[{cell_type}] Step 1: TFBS from ATAC...", flush=True)
atac_region_file = get_TFBS_from_ATAC(
    atac_file, motif_dir, genome_file, overwrite=False, save_dir=save_dir)
print(f"  -> {atac_region_file}", flush=True)

# ── Step 2: TFBS from promoters ───────────────────────────────────────────────
print(f"[{cell_type}] Step 2: TFBS from promoters...", flush=True)
promoter_region_file = get_TFBS_from_promoter(
    promoter_file, motif_dir, overwrite=False, save_dir=save_dir)
print(f"  -> {promoter_region_file}", flush=True)

# ── Step 3: Build TF-gene graph ───────────────────────────────────────────────
print(f"[{cell_type}] Step 3: Build graph...", flush=True)
graph_file, node_file = build_Graph(
    atac_region_file, promoter_region_file, overwrite=False, save_dir=save_dir)
print(f"  graph : {graph_file}", flush=True)
print(f"  nodes : {node_file}", flush=True)

# ── Step 4: Extract features ──────────────────────────────────────────────────
print(f"[{cell_type}] Step 4: Extract features...", flush=True)
atac_feature_file  = extract_atac_feature(
    atac_file, node_file, cell_num=None, seed=666, overwrite=False, save_dir=save_dir)
scrna_feature_file = extract_scrna_feature2(
    scrna_file, node_file, cell_num=None, seed=666, overwrite=False,
    just_node=False, save_dir=save_dir)
print(f"  atac_feat : {atac_feature_file}", flush=True)
print(f"  scrna_feat: {scrna_feature_file}", flush=True)

# ── Step 5: Cross-modal graph ─────────────────────────────────────────────────
print(f"[{cell_type}] Step 5: Cross-modal graph...", flush=True)
scrna_tmp_file = os.path.join(save_dir, "tmp", "scrna.csv")
cross_graph_file, cross_node_file, cross_atac_file, cross_scrna_file = build_cross_graph(
    graph_file, atac_feature_file, scrna_feature_file,
    overwrite=False, save_dir=save_dir)
print(f"  cross_graph: {cross_graph_file}", flush=True)

# ── Step 6: Train scMultiomeGRN ───────────────────────────────────────────────
print(f"[{cell_type}] Step 6: Training scMultiomeGRN (GPU)...", flush=True)
print(f"  torch: {torch.__version__}, CUDA: {torch.cuda.is_available()}", flush=True)

parser = Trainer.add_argparse_args(
    name=cell_type,
    adj_file=cross_graph_file,
    node_file=cross_node_file,
    scrna_file=scrna_tmp_file,
    atac_feature_file=cross_atac_file,
    scrna_feature_file=cross_scrna_file,
    model_cls="scMultiomeGRN",
    comment="DAM-DRUG",
    mask_rate=0.0,
    data_dir=os.path.join(out_dir, "tmp"),
    max_epochs=2000,
    patience=100,
    n_splits=10,
)
config = parser.parse_args([])
config.accelerator = "auto"

log_dir = Trainer.run_single_fold(**vars(config))
print(f"[{cell_type}] Training done. log_dir: {log_dir}", flush=True)

# ── Step 7: Hub TF summary ────────────────────────────────────────────────────
best_adj = pd.read_csv(os.path.join(log_dir, "best_all_adj.csv"), index_col=0)
node_names = best_adj.columns.tolist()
hub_score = best_adj.sum(axis=1).values - np.diag(best_adj.values)
hub_df = pd.DataFrame({"TF": node_names, "n_targets": hub_score})
hub_df = hub_df.sort_values("n_targets", ascending=False).reset_index(drop=True)
hub_df["state"] = "$STATE"
hub_df.to_csv(os.path.join(out_dir, "hub_scores.csv"), index=False)
print(f"Top 10 hub TFs ({cell_type}):")
print(hub_df.head(10).to_string(index=False))

# Highlight target TFs
target_tfs = {"SPI1","RUNX1","IRF8","PPARG","CEBPB","IKZF1","RELB","BHLHE40","BHLHE41"}
found = hub_df[hub_df["TF"].isin(target_tfs)]
print(f"\nTarget TF ranks ({cell_type}):")
print(found.to_string(index=False))
PYEOF

# ── Run inside scMultiomeGRN container ───────────────────────────────────────
apptainer exec --nv \
  --bind "$PROJDIR:$PROJDIR" \
  "$SIF" \
  python "$SCRIPT"

echo ""
echo "=== $STATE done: $(date) ==="
