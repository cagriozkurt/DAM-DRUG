"""
WP1.5 — Direct NES diagnostic: does IKZF1's own JASPAR2026 motif recover its
own GRNBoost2-predicted targets at all?
==============================================================================
`pyscenic ctx` builds several derived TF-modules per TF internally (top-N
targets, importance-quantile thresholds, etc.) and prunes each against every
motif in the database; it doesn't expose per-motif NES scores for modules
that fail to pass, so a "not retained" regulon gives no diagnostic detail on
its own. This script reproduces the recovery-curve / NES calculation
directly via ctxcore's public API (the same one pyscenic ctx uses
internally) for one concrete, best-case module — IKZF1's top-50 targets by
GRNBoost2 importance — against both cisTarget windows, to get an actual NES
number rather than a binary pass/fail.

Result (2026-09-11): IKZF1's own motif (MA1508.2) scores NES = -0.02 (wide
window) / +0.57 (narrow window) for its own top-50 target module — i.e.
indistinguishable from the null (NES~0), nowhere near the ctx default
nes_threshold of 3.0. This rules out both prior hypotheses (window-width
dilution, motif-universe size) as insufficient explanations: even in the
single most favourable case (tightest, highest-confidence target set,
either window), there is no recovery signal for IKZF1's own motif. Other
(non-IKZF1-annotated) motifs DO reach NES>3 for this same target set,
confirming the targets share SOME regulatory signature -- just not one
matching IKZF1's JASPAR2026 PWM.

Run: apptainer exec --env PYTHONPATH=$PROJDIR/tools/pylibs
       --env LD_LIBRARY_PATH=/opt/conda/envs/scenic/lib
       --env LD_PRELOAD=/opt/conda/envs/scenic/lib/libstdc++.so.6
       containers/scenic.sif conda run -n scenic python s1_05_diagnose_nes.py
"""
import pandas as pd, numpy as np
from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.recovery import recovery
from ctxcore.genesig import GeneSignature

WIDE = 'data/resources/jaspar2026/db/jaspar2026_hg38_10kbp_up_10kbp_down_full_tx.genes_vs_motifs.rankings.feather'
NARROW = 'data/resources/jaspar2026/db/jaspar2026_hg38_500bp_up_100bp_down_full_tx.genes_vs_motifs.rankings.feather'

RANK_THRESHOLD = 5000   # ctx default
AUC_THRESHOLD = 0.05    # ctx default
TARGET_TF = 'IKZF1'
TARGET_MOTIF_SUBSTR = 'MA1508'  # IKZF1's motif ID in JASPAR2026

adj = pd.read_csv('results/phase2/GRN/adj_matrix_aggregated.tsv', sep='\t')
tf_rows = adj[adj['TF'] == TARGET_TF].sort_values('importance', ascending=False)
top50 = {g: 1.0 for g in tf_rows.head(50)['target']}
print(f'{TARGET_TF} total targets: {len(tf_rows)}; top50 genes: {len(top50)}')

for label, path in [('WIDE', WIDE), ('NARROW', NARROW)]:
    db = FeatherRankingDatabase(path, name=label)
    gs = GeneSignature(name=f'{TARGET_TF}_top50', gene2weight=top50)
    print(label, 'genes found in db:', len(set(gs.genes) & set(db.genes)), '/', len(top50))
    rnk = db.load(gs)
    weights = np.array([top50.get(g, 0.0) for g in rnk.columns.values])
    _rccs, aucs = recovery(rnk, db.total_genes, weights, RANK_THRESHOLD, AUC_THRESHOLD)
    aucs = pd.Series(aucs, index=rnk.index)
    mean_auc, std_auc = aucs.mean(), aucs.std()
    nes = (aucs - mean_auc) / std_auc
    own_motif_nes = nes[nes.index.str.contains(TARGET_MOTIF_SUBSTR)]
    print(label, 'n_motifs', len(aucs), 'null mean', round(mean_auc, 5), 'std', round(std_auc, 5))
    print(label, f'{TARGET_TF} motif ({TARGET_MOTIF_SUBSTR}) NES:', own_motif_nes.values)
    print(label, 'top5 NES overall (any motif):')
    print(nes.sort_values(ascending=False).head(5))
    print(label, 'n motifs NES>3 (ctx pass threshold):', int((nes > 3).sum()))
