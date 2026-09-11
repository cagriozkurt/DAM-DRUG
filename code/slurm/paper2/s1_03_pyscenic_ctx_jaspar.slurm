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
# Scores against BOTH cisTarget windows built by s1_02 (wide 10kbp-up/10kbp-
# down + narrow 500bp-up/100bp-down), passed together to `pyscenic ctx` (it
# takes multiple db_fnames and scores each motif against all of them, taking
# the max). FINDING 2026-09-11 (job 6346282): the wide window ALONE retained
# only 9 regulons and NONE of the 11 target TFs (IKZF1 included) despite every
# target TF having a present motif and a large GRNBoost2 candidate-target set
# -- diagnosed as enrichment dilution from using one wide window in isolation,
# not a real absence of TF-target relationships (see s1_02's header). Rerun
# with the matched dual-window pair aertslab's own production databases use.
#
# Inputs:
#   results/phase2/GRN/adj_matrix_aggregated.tsv        (5-seed consensus, unchanged)
#   results/phase2/GRN/microglia_raw.loom               (100K-cell expression, unchanged)
#   data/resources/jaspar2026/db/*.rankings.feather     (both windows, from s1_02)
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
SIF=$PROJDIR/containers/scenic.sif
[ -f "$SIF" ] || apptainer pull "$SIF" docker://ghcr.io/cagriozkurt/dam-drug-scanpy@sha256:a16bd5bee57c5fb1918f2e4206368d342861392126c3a4159a462ddd8910efc3

GRN=$PROJDIR/results/phase2/GRN
OUT=$PROJDIR/results/phase7/grn_jaspar2026
mkdir -p "$PROJDIR/logs" "$OUT"

ADJ=$GRN/adj_matrix_aggregated.tsv
LOOM=$GRN/microglia_raw.loom
RANKINGS=($(ls "$PROJDIR"/data/resources/jaspar2026/db/*.genes_vs_motifs.rankings.feather))
MOTIF_TBL=$PROJDIR/data/resources/jaspar2026/jaspar2026_motif2tf.tbl
REG=$OUT/regulons_jaspar2026.csv
AUC=$OUT/scenic_auc_jaspar2026.loom

[ "${#RANKINGS[@]}" -ge 2 ] || { echo "ERROR: expected 2 rankings feathers (wide+narrow) from s1_02, found ${#RANKINGS[@]}: ${RANKINGS[*]}"; exit 1; }
for f in "$ADJ" "$LOOM" "$MOTIF_TBL" "${RANKINGS[@]}"; do
    [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done
echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME — start $(date)"
echo "rankings (${#RANKINGS[@]} dbs): ${RANKINGS[*]}"

# Three env fixes discovered running this 2026-09-11 (the scenic env has
# drifted since the accepted-paper runs used it):
#  1. bare `apptainer exec $SIF pyscenic ...` fails ("executable file not
#     found in $PATH") -- pyscenic lives in the `scenic` conda env, never on
#     the container's base PATH. Must go through `conda run -n scenic`.
#  2. ctxcore imports `pkg_resources`, which setuptools >=81 (installed in
#     the image, 82.0.1) no longer ships. Pin an older setuptools to
#     $TOOLS/pylibs (writable target dir, not the read-only image) + PYTHONPATH
#     (same mechanism s1_02 used for flatbuffers/loguru).
#  3. pyscenic's transform.py still uses `np.object`, removed in current
#     numpy -> AttributeError. Fixed via a sitecustomize.py in the same
#     pylibs dir (auto-imported by every `python`/`conda run` invocation that
#     has it on PYTHONPATH): restores object/bool/int/float/str aliases on
#     the numpy module if missing. Purely additive, no numpy reinstall (which
#     would risk an ABI mismatch with already-compiled pandas/pyarrow).
#  4. loompy -> sqlite3 -> libicui18n.so.78 needs CXXABI_1.3.15 from a newer
#     libstdc++ than the SYSTEM one apptainer's auto library binds put ahead
#     of the conda env's own (newer) libstdc++.so.6, even with
#     LD_LIBRARY_PATH set correctly -- confirmed via `echo $LD_LIBRARY_PATH`
#     inside the container. LD_PRELOAD (higher priority than RPATH) forces
#     the correct one.
#  5. pyscenic's own aucell writer (cli/utils.py append_auc_mtx) calls
#     pandas Series.iteritems(), removed in pandas >=2.0 -> AttributeError,
#     surfaces only at the very end of aucell after ~10 min of enrichment
#     calc (job 6346130, 2026-09-11). Same sitecustomize.py mechanism,
#     additive alias (Series.items already does what iteritems did).
TOOLS=$PROJDIR/tools
PYLIBS=$TOOLS/pylibs
mkdir -p "$PYLIBS"
apptainer exec --bind "$PROJDIR:$PROJDIR" "$SIF" conda run -n scenic \
    pip install --no-input --target "$PYLIBS" "setuptools<81" >/dev/null
cat > "$PYLIBS/sitecustomize.py" << 'PY'
try:
    import numpy as _np
    for _name, _repl in (("object", object), ("bool", bool), ("int", int),
                         ("float", float), ("str", str)):
        if not hasattr(_np, _name):
            setattr(_np, _name, _repl)
except Exception:
    pass

try:
    import pandas as _pd
    if not hasattr(_pd.Series, "iteritems"):
        _pd.Series.iteritems = _pd.Series.items
except Exception:
    pass
PY

BIND="--bind $PROJDIR:$PROJDIR --env LD_LIBRARY_PATH=/opt/conda/envs/scenic/lib --env LD_PRELOAD=/opt/conda/envs/scenic/lib/libstdc++.so.6 --env PYTHONPATH=$PYLIBS"

if [ ! -s "$REG" ]; then
    echo "=== pyscenic ctx (JASPAR 2026, ${#RANKINGS[@]} dbs): $(date) ==="
    apptainer exec $BIND "$SIF" conda run -n scenic pyscenic ctx \
        "$ADJ" "${RANKINGS[@]}" \
        --annotations_fname "$MOTIF_TBL" \
        --expression_mtx_fname "$LOOM" \
        --output "$REG" \
        --num_workers 20 --mode dask_multiprocessing --mask_dropouts
fi

if [ ! -s "$AUC" ]; then
    echo "=== pyscenic aucell: $(date) ==="
    apptainer exec $BIND "$SIF" conda run -n scenic pyscenic aucell \
        "$LOOM" "$REG" --output "$AUC" --num_workers 20
fi

echo "=== summarise regulon activity by substate: $(date) ==="
cat > "$OUT/_summarise.py" << 'PY'
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
apptainer exec $BIND "$SIF" conda run -n scenic python "$OUT/_summarise.py"

echo "Done $(date). Next: sbatch code/slurm/paper2/s1_04_benchmark_nulls.slurm"
