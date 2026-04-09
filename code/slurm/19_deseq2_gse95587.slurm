#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=deseq2_gse95587
#SBATCH --output=/arf/scratch/mozkurt/DAM-DRUG/logs/deseq2_gse95587-%j.out
#SBATCH --error=/arf/scratch/mozkurt/DAM-DRUG/logs/deseq2_gse95587-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
# ─────────────────────────────────────────────────────────────────────────────
# Phase 5 — Bulk RNA-seq validation: GSE95587 (Friedman 2018, Genentech)
#   Fusiform gyrus, 117 samples (AD vs neurologically-normal control)
#   Fields available: diagnosis, sex, age at death, Braak I-VI, RIN, PMI
#
# Steps:
#   1. Parse series matrix → metadata.csv
#   2. Merge per-sample TSVs  → raw_counts.tsv  (Entrez ID × GSM)
#   3. pyDESeq2              → deseq2_all_genes.csv + deseq2_target_tfs.csv
#   4. Forest plot           → forest_plot_target_tfs.pdf + .png
#
# Packages required (pip install): pydeseq2 pandas matplotlib adjustText
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

DATA_DIR=/arf/scratch/mozkurt/DAM-DRUG/data/raw/Mayo/suppl
OUT_DIR=/arf/scratch/mozkurt/DAM-DRUG/results/phase5/gse95587
mkdir -p "$OUT_DIR"

echo "Start: $(date)"
echo "Node:  $SLURMD_NODENAME"

# ── Environment ──────────────────────────────────────────────────────────────
source /etc/profile.d/modules.sh
module load miniconda3
source /arf/sw/comp/python/miniconda3/etc/profile.d/conda.sh
conda activate mmpbsa || { echo "ERROR: cannot activate mmpbsa env"; exit 1; }
export PATH="$HOME/.local/bin:$PATH"

# Install packages if needed (no-op if already installed)
pip install -q pydeseq2==0.5.4 adjustText 2>/dev/null || true

# ── Main Python pipeline ─────────────────────────────────────────────────────
python3 << 'PYEOF'
import gzip, sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings('ignore')

DATA_DIR = Path('/arf/scratch/mozkurt/DAM-DRUG/data/raw/Mayo/suppl')
OUT_DIR  = Path('/arf/scratch/mozkurt/DAM-DRUG/results/phase5/gse95587')

# ─────────────────────────────────────────────────────────────────────────────
# Target TF Entrez gene IDs (Homo sapiens)
# ─────────────────────────────────────────────────────────────────────────────
TF_ENTREZ = {
    'IRF8':    '3394',
    'PPARG':   '5468',
    'RUNX1':   '861',
    'BHLHE41': '79365',
    'BHLHE40': '8553',
    'SPI1':    '6688',
    'IKZF1':   '10320',
    'CEBPB':   '1051',
    'MAF':     '4094',
    'RELB':    '5971',
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: parse series matrix
# ─────────────────────────────────────────────────────────────────────────────
print("=== Step 1: parse series matrix ===")

series_file = DATA_DIR / 'GSE95587_series_matrix.txt.gz'
fields = {}       # field_name → list of values (one per sample)
char_lines = []   # all !Sample_characteristics_ch1 rows

with gzip.open(series_file, 'rt') as f:
    for line in f:
        line = line.rstrip('\n')
        if not line.startswith('!Sample_'): continue
        parts = line.split('\t')
        key   = parts[0].lstrip('!')
        vals  = [v.strip('"') for v in parts[1:]]
        if key == 'Sample_characteristics_ch1':
            char_lines.append(vals)
        else:
            fields[key] = vals

gsm_ids = fields['Sample_geo_accession']
n = len(gsm_ids)
print(f"  {n} samples: {gsm_ids[0]} … {gsm_ids[-1]}")

# Parse characteristics (each line is "key: value" for all samples)
chars = {}
for row in char_lines:
    if not row: continue
    first = row[0]
    if ': ' not in first: continue
    char_key = first.split(': ', 1)[0]
    chars[char_key] = [v.split(': ', 1)[1] if ': ' in v else '' for v in row]

print(f"  Characteristic fields: {list(chars.keys())}")

roman_map = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6}
meta = pd.DataFrame({
    'gsm':       gsm_ids,
    'sample_id': chars.get('sample_id', [''] * n),
    'diagnosis': chars.get('diagnosis', [''] * n),
    'sex':       chars.get('Sex', [''] * n),
    'age':       pd.to_numeric(chars.get('age at death', [None]*n), errors='coerce'),
    'braak_str': chars.get('braak score', [''] * n),
    'rin':       pd.to_numeric(chars.get('rin_score', [None]*n), errors='coerce'),
    'pmi':       pd.to_numeric(chars.get('post_mortem_interval', [None]*n), errors='coerce'),
})
meta['condition']  = meta['diagnosis'].apply(
    lambda x: 'AD' if "Alzheimer" in x else ('CON' if 'control' in x.lower() else 'unknown'))
meta['braak_num']  = meta['braak_str'].map(roman_map)

meta.to_csv(OUT_DIR / 'metadata.csv', index=False)
print(meta['condition'].value_counts().to_string())
print(f"  Written: metadata.csv")

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: build raw count matrix
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Step 2: build count matrix ===")

# GSM id → file path
gsm_to_file = {}
for f in sorted(DATA_DIR.glob('GSM*.tsv.gz')):
    gsm = f.name.split('_')[0]
    gsm_to_file[gsm] = f

print(f"  Count files found: {len(gsm_to_file)}")
missing = [g for g in gsm_ids if g not in gsm_to_file]
if missing:
    print(f"  WARN: {len(missing)} GSMs have no count file — will be 0-filled")

# Build sparse-style dict first to avoid huge memory
all_genes = set()
sample_counts = {}  # gsm → {gene: count}

for gsm in gsm_ids:
    if gsm not in gsm_to_file:
        sample_counts[gsm] = {}
        continue
    sc = {}
    with gzip.open(gsm_to_file[gsm], 'rt') as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('ID_REF'): continue
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            gene = parts[0]
            try:
                cnt = int(float(parts[2]))
            except ValueError:
                cnt = 0
            sc[gene] = cnt
            all_genes.add(gene)
    sample_counts[gsm] = sc

print(f"  Genes detected: {len(all_genes)}")

# Assemble DataFrame (genes × samples)
gene_list  = sorted(all_genes, key=lambda x: int(x) if x.isdigit() else 0)
counts_mat = pd.DataFrame(
    {gsm: [sample_counts[gsm].get(g, 0) for g in gene_list] for gsm in gsm_ids},
    index=gene_list
)
counts_mat.index.name = 'entrez_id'
counts_mat.to_csv(OUT_DIR / 'raw_counts.tsv', sep='\t')
print(f"  Count matrix: {counts_mat.shape[0]} genes × {counts_mat.shape[1]} samples")
print(f"  Written: raw_counts.tsv")

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: pyDESeq2
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Step 3: pyDESeq2 ===")

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Keep only AD and CON, drop samples with NA covariates
meta_deseq = meta[meta['condition'].isin(['AD', 'CON'])].copy()
meta_deseq = meta_deseq.dropna(subset=['age', 'rin', 'sex'])
meta_deseq = meta_deseq.set_index('gsm')

# Recode sex as 0/1
meta_deseq['sex_bin'] = (meta_deseq['sex'] == 'M').astype(int)

counts_deseq = counts_mat[meta_deseq.index].T  # samples × genes
print(f"  Samples in DESeq2: {len(meta_deseq)}  (AD={sum(meta_deseq.condition=='AD')}, CON={sum(meta_deseq.condition=='CON')})")

# Filter low-count genes: ≥ 10 reads in at least 10 samples
keep = (counts_deseq >= 10).sum(axis=0) >= 10
counts_deseq = counts_deseq.loc[:, keep]
print(f"  Genes after count filter: {counts_deseq.shape[1]}")

# Check target TFs are present
for gene_name, eid in TF_ENTREZ.items():
    if eid not in counts_deseq.columns:
        print(f"  WARN: {gene_name} (Entrez {eid}) filtered out (low count)")

# Build DESeq dataset
# Continuous covariates: age (z-scored), rin (z-scored); categorical: condition, sex
meta_deseq['age_z'] = (meta_deseq['age'] - meta_deseq['age'].mean()) / meta_deseq['age'].std()
meta_deseq['rin_z']  = (meta_deseq['rin'] - meta_deseq['rin'].mean()) / meta_deseq['rin'].std()

dds = DeseqDataSet(
    counts   = counts_deseq,
    metadata = meta_deseq[['condition', 'sex_bin', 'age_z', 'rin_z']],
    design_factors = ['condition', 'sex_bin', 'age_z', 'rin_z'],
    ref_level = [['condition', 'CON']],
    n_cpus   = 20,
    quiet    = False,
)

dds.deseq2()

# Extract results: AD vs CON
stat_res = DeseqStats(dds, contrast=['condition', 'AD', 'CON'], alpha=0.05, quiet=False)
stat_res.summary()
res_df = stat_res.results_df.copy()
res_df.index.name = 'entrez_id'
res_df = res_df.sort_values('pvalue', na_position='last')
res_df.to_csv(OUT_DIR / 'deseq2_all_genes.csv')
n_sig = (res_df['padj'] < 0.05).sum()
print(f"\n  Significant genes (padj<0.05): {n_sig} / {len(res_df)}")

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: target TF results table
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Step 4: target TF results ===")

rows = []
for gene_name, eid in TF_ENTREZ.items():
    if eid in res_df.index:
        r = res_df.loc[eid]
        rows.append({
            'gene':      gene_name,
            'entrez_id': eid,
            'baseMean':  round(r['baseMean'], 2),
            'log2FC':    round(r['log2FoldChange'], 3),
            'lfcSE':     round(r['lfcSE'], 3),
            'stat':      round(r['stat'], 3),
            'pvalue':    r['pvalue'],
            'padj':      r['padj'],
        })
    else:
        rows.append({
            'gene': gene_name, 'entrez_id': eid,
            'baseMean': None, 'log2FC': None, 'lfcSE': None,
            'stat': None, 'pvalue': None, 'padj': None,
        })

tf_df = pd.DataFrame(rows)
tf_df.to_csv(OUT_DIR / 'deseq2_target_tfs.csv', index=False)

print(tf_df[['gene', 'log2FC', 'pvalue', 'padj']].to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: forest plot
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Step 5: forest plot ===")

plot_df = tf_df.dropna(subset=['log2FC', 'lfcSE']).copy()
plot_df = plot_df.sort_values('log2FC')
plot_df['sig'] = plot_df['padj'].apply(
    lambda p: '**' if (p is not None and not np.isnan(p) and p < 0.01)
              else ('*' if (p is not None and not np.isnan(p) and p < 0.05) else ''))
plot_df['color'] = plot_df.apply(
    lambda r: '#d62728' if r['log2FC'] > 0 and r['sig'] != ''
              else ('#1f77b4' if r['log2FC'] < 0 and r['sig'] != '' else '#aaaaaa'),
    axis=1)

fig, ax = plt.subplots(figsize=(7, 0.5 * len(plot_df) + 2))
y_pos = range(len(plot_df))

for i, (_, row) in enumerate(plot_df.iterrows()):
    ci = 1.96 * row['lfcSE']
    ax.plot([row['log2FC'] - ci, row['log2FC'] + ci], [i, i],
            color=row['color'], linewidth=2, zorder=2)
    ax.scatter(row['log2FC'], i, color=row['color'], s=60, zorder=3)
    if row['sig']:
        ax.text(row['log2FC'] + ci + 0.02, i, row['sig'],
                va='center', fontsize=9, color=row['color'])

ax.axvline(0, color='grey', linestyle='--', linewidth=0.8, zorder=1)
ax.set_yticks(list(y_pos))
ax.set_yticklabels(plot_df['gene'].tolist(), fontsize=11)
ax.set_xlabel('log$_2$ fold change  (AD vs Control)', fontsize=11)
ax.set_title('Target TF expression — GSE95587\nFusiform gyrus, AD vs Control\n'
             'Design: condition + sex + age + RIN', fontsize=10)

patch_up  = mpatches.Patch(color='#d62728', label='Up (padj<0.05)')
patch_dn  = mpatches.Patch(color='#1f77b4', label='Down (padj<0.05)')
patch_ns  = mpatches.Patch(color='#aaaaaa', label='ns')
ax.legend(handles=[patch_up, patch_dn, patch_ns], fontsize=9,
          loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

for ext in ('pdf', 'png'):
    p = OUT_DIR / f'forest_plot_target_tfs.{ext}'
    plt.savefig(p, dpi=150 if ext == 'png' else None)
    print(f"  Written: {p.name}")

plt.close()
print("\nAll done.")
PYEOF

echo "Finished: $(date)"
echo "Results in: $OUT_DIR"
