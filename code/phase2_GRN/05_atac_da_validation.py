"""
DAM-DRUG Item 2.8 (revised): ATAC Differential Accessibility Validation
==========================================================================
Strategy:
  1. Filter MTG ATAC h5ad to microglia (Class == 'Microglia')
  2. Pseudobulk by donor: sum peak counts, split Severely Affected vs Neurotypical
  3. Wilcoxon test per peak → DAPs (BH FDR < 0.05, |log2FC| > 1)
  4. Extract DAP sequences from hg38 genome (via GetSequence.pl)
  5. FIMO scan with 9 target TF motifs from HOCOMOCOv11
  6. Report: enrichment of target TF motifs in AD-upregulated peaks

Requires scMultiomeGRN container (has FIMO, bedops, scipy, pandas):
  apptainer exec containers/scmultiomegrn.sif python code/phase2_GRN/05_atac_da_validation.py

Or submit via: sbatch code/slurm/13_atac_da_validation.slurm
"""

import os
import sys
import re
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as stats
from pathlib import Path
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
ATAC_H5AD = PROJECT / "data/raw/SEA-AD/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
TOOLDIR  = PROJECT / "tools/scMultiomeGRN/extracted/ScmultiomeGRN-main"
DATADIR  = PROJECT / "tools/scMultiomeGRN/data_resource"
OUT      = PROJECT / "results/phase2/ATAC_validation"
CKPT     = OUT / ".checkpoints"

MOTIF_DIR    = DATADIR / "HOCOMOCOv11"
GENOME_DIR   = DATADIR / "download_genome/hg38"
GET_SEQ_PERL = TOOLDIR / "src/GetSequence.pl"
PROMOTER_FILE = TOOLDIR / "data_resource/gencode.v38.ProteinCoding_gene_promoter.txt"

# Target TFs from DAM-DRUG project
TARGET_TFS = ["SPI1", "RUNX1", "IRF8", "PPARG", "CEBPB", "IKZF1", "RELB", "BHLHE40", "BHLHE41"]

OUT.mkdir(parents=True, exist_ok=True)
CKPT.mkdir(exist_ok=True)


def run(cmd, verbose=False):
    if verbose:
        print(f"  $ {cmd}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\n{r.stderr}")
    return r.stdout.strip()


# ── Step 1: Filter ATAC to microglia, pseudobulk ──────────────────────────────
PSEUDO_FILE = OUT / "pseudobulk_microglia.csv"

if PSEUDO_FILE.exists():
    print(f"[Step 1] SKIP — pseudobulk already exists ({PSEUDO_FILE})")
    pseudo = pd.read_csv(PSEUDO_FILE, index_col=0)
    meta   = pd.read_csv(OUT / "pseudobulk_meta.csv", index_col=0)
else:
    print("[Step 1] Loading ATAC h5ad (backed)...", flush=True)
    atac = sc.read_h5ad(ATAC_H5AD, backed='r')
    print(f"  Total: {atac.n_obs:,} cells × {atac.n_vars:,} peaks")

    # Print cell type distribution
    print(f"  Class values: {atac.obs['Class'].value_counts().head(10).to_dict()}")

    # Filter to microglia (Class="Non-neuronal and Non-neural", Subclass="Microglia-PVM")
    mg_mask = atac.obs["Subclass"] == "Microglia-PVM"
    mg_barcodes = atac.obs_names[mg_mask]
    print(f"  Microglia: {len(mg_barcodes):,} cells")

    # Print subtypes
    mg_obs = atac.obs.loc[mg_barcodes]
    print(f"  Subclass: {mg_obs['Subclass'].value_counts().to_dict()}")
    print(f"  Supertype unique: {mg_obs['Supertype'].nunique()} supertypes")
    print(f"  Neurotypical reference: {mg_obs['Neurotypical reference'].value_counts().to_dict()}")
    print(f"  Severely Affected: {mg_obs['Severely Affected Donor'].value_counts().to_dict()}")

    # Load microglia ATAC matrix into RAM
    print("  Loading microglia ATAC into RAM...", flush=True)
    mg = atac[mg_barcodes].to_memory()
    atac.file.close()

    # Pseudobulk: sum peaks per donor
    print("  Pseudobulking by donor...", flush=True)
    donors    = mg_obs["Donor ID"].values
    neuro_ref = mg_obs["Neurotypical reference"].values
    braak     = mg_obs["Braak"].values
    sev       = mg_obs["Severely Affected Donor"].values

    unique_donors = mg_obs["Donor ID"].unique()
    pb_rows = []
    pb_meta = []
    X = mg.X if sp.issparse(mg.X) else sp.csr_matrix(mg.X)

    for donor in tqdm(unique_donors, desc="Pseudobulk"):
        dmask = mg_obs["Donor ID"] == donor
        dX = X[dmask.values]
        pb_rows.append(np.asarray(dX.sum(axis=0)).ravel())
        row = mg_obs.loc[dmask].iloc[0]
        pb_meta.append({
            "donor_id":      donor,
            "n_cells":       int(dmask.sum()),
            "neurotypical":  row["Neurotypical reference"],
            "braak":         row["Braak"],
            "severely_affected": row["Severely Affected Donor"],
        })

    pseudo = pd.DataFrame(np.array(pb_rows), index=unique_donors, columns=mg.var_names)
    meta   = pd.DataFrame(pb_meta, index=unique_donors)

    pseudo.to_csv(PSEUDO_FILE)
    meta.to_csv(OUT / "pseudobulk_meta.csv")
    print(f"  Pseudobulk: {pseudo.shape[0]} donors × {pseudo.shape[1]} peaks")


# ── Step 2: Differential accessibility test ────────────────────────────────────
DA_FILE = OUT / "da_peaks.csv"

if DA_FILE.exists():
    print(f"\n[Step 2] SKIP — DA results already exist ({DA_FILE})")
    da = pd.read_csv(DA_FILE)
else:
    print("\n[Step 2] Differential accessibility (Wilcoxon, pseudobulk)...", flush=True)

    # Group: Severely Affected (Y) vs Not Severely Affected (N)
    # Only 2 Neurotypical donors in this ATAC dataset; severity is better powered
    sev_val = meta["severely_affected"].astype(str).str.strip()
    ad_mask   = (sev_val == "Y").values
    ctrl_mask = (sev_val == "N").values

    print(f"  Severely Affected (AD) donors: {ad_mask.sum()}")
    print(f"  Not Severely Affected donors : {ctrl_mask.sum()}")

    if ctrl_mask.sum() < 3 or ad_mask.sum() < 3:
        print("  WARNING: too few donors in one group.")
        print(f"  severely_affected unique values: {meta['severely_affected'].unique()}")
        sys.exit(1)

    # Normalize: log1p(CPM)
    total = pseudo.sum(axis=1).values
    normed = np.log1p(pseudo.values / total[:, None] * 1e6)

    ctrl_mean = normed[ctrl_mask].mean(axis=0)
    ad_mean   = normed[ad_mask].mean(axis=0)
    log2fc    = ad_mean - ctrl_mean  # log1p scale, interpretable as ~log2FC

    # Wilcoxon per peak
    pvals = np.zeros(pseudo.shape[1])
    for j in tqdm(range(pseudo.shape[1]), desc="Wilcoxon", mininterval=10):
        try:
            _, pvals[j] = stats.mannwhitneyu(
                normed[ad_mask, j], normed[ctrl_mask, j],
                alternative="two-sided"
            )
        except Exception:
            pvals[j] = 1.0

    # BH FDR
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")

    da = pd.DataFrame({
        "peak":    pseudo.columns,
        "log2fc":  log2fc,
        "pval":    pvals,
        "padj":    padj,
        "ctrl_mean": ctrl_mean,
        "ad_mean":   ad_mean,
    })
    da.to_csv(DA_FILE, index=False)

    n_sig = ((da["padj"] < 0.05) & (da["log2fc"].abs() > 1)).sum()
    n_up  = ((da["padj"] < 0.05) & (da["log2fc"] > 1)).sum()
    n_dn  = ((da["padj"] < 0.05) & (da["log2fc"] < -1)).sum()
    print(f"  Significant DAPs (FDR<0.05, |log2FC|>1): {n_sig}  (up={n_up}, down={n_dn})")


# ── Step 3: BED file of AD-upregulated DAPs ───────────────────────────────────
UP_BED = OUT / "daps_ad_up.bed"

if UP_BED.exists():
    print(f"\n[Step 3] SKIP — BED file exists ({UP_BED})")
else:
    print("\n[Step 3] Writing BED of AD-upregulated peaks...", flush=True)
    da = pd.read_csv(DA_FILE)
    dap_up = da[(da["padj"] < 0.05) & (da["log2fc"] > 1)].copy()
    print(f"  AD-upregulated DAPs: {len(dap_up)}")

    if len(dap_up) == 0:
        print("  WARNING: No significant AD-upregulated peaks. Relaxing threshold to top 1000 by log2FC.")
        dap_up = da.sort_values("log2fc", ascending=False).head(1000)

    def peak_to_bed(peak_name):
        """chr1:12345-67890 or chr1_12345_67890 → (chr, start, end)"""
        clean = re.sub(r"[:\-]", "_", peak_name)
        parts = clean.split("_")
        if len(parts) >= 3:
            return parts[0], parts[-2], parts[-1]
        return None, None, None

    bed_rows = []
    for peak in dap_up["peak"]:
        chrom, start, end = peak_to_bed(peak)
        if chrom:
            bed_rows.append(f"{chrom}\t{start}\t{end}\t{peak}")

    with open(UP_BED, "w") as f:
        f.write("\n".join(bed_rows) + "\n")
    print(f"  Written: {UP_BED} ({len(bed_rows)} peaks)")

    # Also background BED (all tested peaks, FDR>0.2)
    bg = da[da["padj"] > 0.2].sample(min(5000, (da["padj"]>0.2).sum()), random_state=42)
    bg_rows = []
    for peak in bg["peak"]:
        chrom, start, end = peak_to_bed(peak)
        if chrom:
            bg_rows.append(f"{chrom}\t{start}\t{end}\t{peak}")
    with open(OUT / "background_peaks.bed", "w") as f:
        f.write("\n".join(bg_rows) + "\n")


# ── Step 4: Extract FASTA sequences ───────────────────────────────────────────
UP_FASTA  = OUT / "daps_ad_up.fasta"
BG_FASTA  = OUT / "background_peaks.fasta"

def bed_to_fasta(bed_file, fasta_file):
    if fasta_file.exists():
        print(f"  SKIP — {fasta_file.name} exists")
        return
    perl_loc = run("which perl")
    cmd = (f"{perl_loc} {GET_SEQ_PERL} "
           f"{bed_file.resolve()} {fasta_file.resolve()} {GENOME_DIR.resolve()}")
    run(cmd, verbose=True)

print("\n[Step 4] Extracting peak sequences...", flush=True)
bed_to_fasta(UP_BED,               UP_FASTA)
bed_to_fasta(OUT / "background_peaks.bed", BG_FASTA)


# ── Step 5: FIMO scan for target TF motifs ────────────────────────────────────
FIMO_DIR = OUT / "fimo_results"
FIMO_DIR.mkdir(exist_ok=True)

print("\n[Step 5] FIMO motif scanning...", flush=True)
FIMO_LOC = run("which fimo")

# Find motif files for target TFs
motif_files = list(MOTIF_DIR.glob("*.meme"))
target_motifs = {}
for mf in motif_files:
    tf_name = mf.stem.split("_HUMAN")[0].split(".")[0]
    if tf_name in TARGET_TFS:
        target_motifs.setdefault(tf_name, []).append(mf)

print(f"  Found motif files for: {sorted(target_motifs.keys())}")

results = []
for tf, mfiles in sorted(target_motifs.items()):
    for mf in mfiles:
        motif_tag = mf.stem
        for label, fasta in [("ad_up", UP_FASTA), ("background", BG_FASTA)]:
            out_tsv = FIMO_DIR / f"{motif_tag}_{label}.tsv"
            if out_tsv.exists():
                continue
            tmp_dir = FIMO_DIR / f"tmp_{motif_tag}_{label}"
            tmp_dir.mkdir(exist_ok=True)
            try:
                run(f"{FIMO_LOC} --oc {tmp_dir} --thresh 1e-4 --no-qvalue "
                    f"{mf} {fasta}", verbose=False)
                fimo_tsv = tmp_dir / "fimo.tsv"
                if fimo_tsv.exists():
                    fimo_tsv.rename(out_tsv)
                run(f"rm -rf {tmp_dir}")
            except Exception as e:
                print(f"  FIMO error {motif_tag} {label}: {e}")

# Count hits per TF per set
print("\n[Step 6] Summarizing motif enrichment...", flush=True)
summary_rows = []
for tf in TARGET_TFS:
    for mf in target_motifs.get(tf, []):
        motif_tag = mf.stem
        def count_hits(label):
            tsv = FIMO_DIR / f"{motif_tag}_{label}.tsv"
            if not tsv.exists():
                return 0
            try:
                df = pd.read_csv(tsv, sep="\t", comment="#")
                return len(df[df["q-value"] < 0.05]) if "q-value" in df.columns else len(df)
            except Exception:
                return 0

        n_up = count_hits("ad_up")
        n_bg = count_hits("background")

        # Enrichment ratio (hits per peak)
        with open(UP_BED) as _f:
            n_up_peaks = sum(1 for _ in _f)
        with open(OUT / "background_peaks.bed") as _f:
            n_bg_peaks = sum(1 for _ in _f)
        rate_up = n_up / max(n_up_peaks, 1)
        rate_bg = n_bg / max(n_bg_peaks, 1)
        enrichment = rate_up / max(rate_bg, 1e-6)

        summary_rows.append({
            "TF": tf, "motif": motif_tag,
            "hits_ad_up": n_up, "hits_background": n_bg,
            "rate_up": rate_up, "rate_bg": rate_bg,
            "enrichment_ratio": enrichment,
        })

summary = pd.DataFrame(summary_rows)
if not summary.empty:
    summary = summary.sort_values("enrichment_ratio", ascending=False)
    summary.to_csv(OUT / "motif_enrichment_summary.csv", index=False)
    print("\nMotif enrichment (AD-upregulated peaks vs background):")
    print(summary[["TF", "hits_ad_up", "hits_background", "enrichment_ratio"]].to_string(index=False))

print(f"\n=== ATAC DA validation complete ===")
print(f"Outputs in: {OUT}")
print(f"  da_peaks.csv             — {(OUT/'da_peaks.csv').stat().st_size//1024} KB")
print(f"  motif_enrichment_summary.csv")
