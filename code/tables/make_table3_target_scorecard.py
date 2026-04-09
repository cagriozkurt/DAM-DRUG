"""
DAM-DRUG — Table 3: Target Scorecard
======================================
Compiles multi-evidence support for each of the 10 candidate TF/enzyme targets.

Evidence columns:
  pySCENIC_regulon  — in the 46-regulon pySCENIC GRN set (Y/N)
  AUC_LateAD_DAM    — mean AUCell score in LateAD-DAM state (SCENIC regulon members only)
  pseudotime_rho    — Spearman ρ with diffusion pseudotime (regulon members only)
  scRNA_DGE         — best disease-vs-homeostatic log2FC from pseudobulk DESeq2
  bulk_log2FC/padj  — GSE95587 fusiform gyrus replication (n=117)
  ATAC_enrichment   — fold-enrichment of TF motif in AD-upregulated ATAC peaks
  CellOracle_KO_mag — mean perturbation magnitude in LateAD-DAM on TF KO
  fpocket_score     — fpocket druggability drug_score (best structure)
  strategy          — proposed therapeutic modality

Run locally:
  python code/tables/make_table3_target_scorecard.py
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", "/Volumes/PortableSSD/untitled folder/DAM-DRUG"))
OUT_CSV = PROJECT / "results/tables/table3_target_scorecard.csv"
OUT_MD  = PROJECT / "results/tables/table3_target_scorecard.md"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────────
drugg    = pd.read_csv(PROJECT / "results/phase3/druggability_ranking.csv", index_col=0)
auc      = pd.read_csv(PROJECT / "results/phase2/GRN/regulon_auc_by_state_aggregated.csv",
                        index_col=0)
pseudotime = pd.read_csv(PROJECT / "results/phase2/GRN/regulon_pseudotime_corr.csv")
gse      = pd.read_csv(PROJECT / "results/phase5/gse95587/deseq2_target_tfs.csv")
dge_pb   = pd.read_csv(PROJECT / "results/phase1/DGE_pseudobulk/TF_pseudobulk_summary.csv")

# ── Hardcoded values not available as local CSVs ───────────────────────────────
# ATAC enrichment: TF motif fold-enrichment in 56 AD-upregulated DAPs (step 2.8)
atac = {"PPARG": 4.5, "IKZF1": 3.9}

# CellOracle KO magnitude in LateAD-DAM state (step 5.7; three TFs tested)
ko_mag = {"IKZF1": 0.047, "IRF8": 0.026, "SPI1": 0.013}

# Literature support (brief tag for strategy column)
strategy = {
    "PPARG":   "Agonist (LBD); LDAM lipid regulator",
    "IKZF1":   "PROTAC/degrader (ZF2 degron); no small-mol pocket",
    "IRF8":    "Small-molecule (DBD pocket); DAM master regulator",
    "BHLHE41": "Small-molecule (bHLH); repressor of DAM transition",
    "SPI1":    "Small-molecule (ETS domain); causal DAM driver (Ayata 2025)",
    "RUNX1":   "Small-molecule (Runt domain); DAM co-activator",
    "RUNX2":   "Activator strategy (Runt); low druggability — deprioritize",
    "MAF":     "Small-molecule (bZIP); MafA homolog pocket",
    "ACSL1":   "Metabolic inhibitor (AMP-binding); LDAM fatty-acid enzyme (Haney 2024)",
    "PIK3CA":  "Kinase inhibitor (ATP site); PI3K-α LDAM lipid metabolism (Haney 2024)",
}

TARGETS = ["PPARG", "IKZF1", "IRF8", "BHLHE41", "SPI1",
           "RUNX1", "RUNX2", "MAF", "ACSL1", "PIK3CA"]

# ── Best DGE contrast per TF ────────────────────────────────────────────────────
# Use the disease-vs-homeostatic contrast where TF is most significantly upregulated.
# Priority: DAM > IRM > LAM (in disease-progression order).
DGE_CONTRASTS = ["DAM_vs_Homeostatic", "IRM_vs_Homeostatic", "LAM_vs_Homeostatic"]

def best_dge(tf):
    """Return (contrast, log2FC, padj) for most significant upregulated contrast."""
    best = None
    for contrast in DGE_CONTRASTS:
        row = dge_pb[(dge_pb["names"] == tf) & (dge_pb["contrast"] == contrast)]
        if row.empty:
            continue
        lfc  = row.iloc[0]["log2FoldChange"]
        padj = row.iloc[0]["padj"]
        if lfc > 0 and (best is None or padj < best[2]):
            best = (contrast, lfc, padj)
    return best  # None if not upregulated in any contrast


# ── Build table ────────────────────────────────────────────────────────────────
rows = []
for tf in TARGETS:
    r = {"TF": tf}

    # pySCENIC regulon membership
    col = f"{tf}(+)"
    in_grn = col in auc.columns
    r["pySCENIC_regulon"] = "Y" if in_grn else "N"

    # AUCell score in LateAD-DAM (SCENIC members only)
    if in_grn and "LateAD-DAM" in auc.index:
        r["AUC_LateAD_DAM"] = round(auc.loc["LateAD-DAM", col], 4)
    else:
        r["AUC_LateAD_DAM"] = None

    # Pseudotime ρ (SCENIC members only)
    pt = pseudotime[pseudotime["TF"] == tf]
    if not pt.empty:
        r["pseudotime_rho"]  = round(pt.iloc[0]["rho"], 3)
        r["pseudotime_padj"] = float(f"{pt.iloc[0]['padj']:.2e}")
    else:
        r["pseudotime_rho"]  = None
        r["pseudotime_padj"] = None

    # Internal scRNA pseudobulk DGE
    dge_result = best_dge(tf)
    if dge_result:
        contrast, lfc, padj = dge_result
        r["scRNA_contrast"]  = contrast.replace("_vs_Homeostatic", "")
        r["scRNA_log2FC"]    = round(lfc, 2)
        r["scRNA_padj"]      = float(f"{padj:.2e}")
    else:
        # Check if any contrast exists at all (might be repressed in all)
        any_row = dge_pb[dge_pb["names"] == tf]
        r["scRNA_contrast"]  = "—"
        r["scRNA_log2FC"]    = None
        r["scRNA_padj"]      = None

    # GSE95587 bulk replication
    g = gse[gse["gene"] == tf]
    if not g.empty:
        r["bulk_log2FC"] = round(g.iloc[0]["log2FC"], 3)
        r["bulk_padj"]   = round(g.iloc[0]["padj"], 4)
    else:
        r["bulk_log2FC"] = None
        r["bulk_padj"]   = None

    # ATAC enrichment (hardcoded from step 2.8)
    r["ATAC_enrichment"] = atac.get(tf, None)

    # CellOracle KO magnitude (hardcoded from step 5.7)
    r["CellOracle_KO_LateAD_DAM"] = ko_mag.get(tf, None)

    # fpocket druggability
    dr = drugg[drugg["tf"] == tf]
    if not dr.empty:
        r["fpocket_drug_score"] = round(dr.iloc[0]["drug_score"], 3)
        r["fpocket_structure"]  = dr.iloc[0]["structure_id"].replace("_prep", "")
        r["fpocket_rank"]       = int(dr.index[0])
    else:
        r["fpocket_drug_score"] = None
        r["fpocket_structure"]  = None
        r["fpocket_rank"]       = None

    r["strategy"] = strategy.get(tf, "—")
    rows.append(r)

table = pd.DataFrame(rows)
table.to_csv(OUT_CSV, index=False)
print(f"Saved CSV → {OUT_CSV}")
print()
print(table[["TF", "pySCENIC_regulon", "AUC_LateAD_DAM", "pseudotime_rho",
             "scRNA_contrast", "scRNA_log2FC", "scRNA_padj",
             "bulk_log2FC", "bulk_padj",
             "ATAC_enrichment", "CellOracle_KO_LateAD_DAM",
             "fpocket_drug_score"]].to_string(index=False))

# ── Markdown version (manuscript-ready) ──────────────────────────────────────
def fmt(val, fmt_str=None, suffix=""):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return "—"
    if fmt_str:
        return format(val, fmt_str) + suffix
    return str(val) + suffix

md_rows = []
for _, r in table.iterrows():
    # Significance stars for padj
    def sig(padj):
        if padj is None: return ""
        if padj < 0.001: return "***"
        if padj < 0.01:  return "**"
        if padj < 0.05:  return "*"
        return " (ns)"

    _no_dge = r["scRNA_log2FC"] is None or (isinstance(r["scRNA_log2FC"], float) and np.isnan(r["scRNA_log2FC"]))
    scRNA_cell = ("—" if _no_dge else
                  f"{fmt(r['scRNA_log2FC'], '+.2f')} [{r['scRNA_contrast']}]{sig(r['scRNA_padj'])}")
    _no_bulk = r["bulk_log2FC"] is None or (isinstance(r["bulk_log2FC"], float) and np.isnan(r["bulk_log2FC"]))
    bulk_cell  = ("—" if _no_bulk else
                  f"{fmt(r['bulk_log2FC'], '+.3f')}{sig(r['bulk_padj'])}")
    atac_cell  = fmt(r["ATAC_enrichment"], ".1f", "×")
    ko_cell    = fmt(r["CellOracle_KO_LateAD_DAM"], ".3f")
    auc_cell   = fmt(r["AUC_LateAD_DAM"], ".4f")
    pt_cell    = fmt(r["pseudotime_rho"], "+.3f")
    grn_cell   = r["pySCENIC_regulon"]
    # RUNX2: fpocket omitted (AF2 isoform truncated; activation not inhibition strategy)
    if r["TF"] == "RUNX2":
        drug_cell = "n.d.†"
    else:
        drug_cell  = fmt(r["fpocket_drug_score"], ".3f")

    md_rows.append(f"| {r['TF']} | {grn_cell} | {auc_cell} | {pt_cell} | "
                   f"{scRNA_cell} | {bulk_cell} | {atac_cell} | {ko_cell} | "
                   f"{drug_cell} |")

header = ("| TF | GRN¹ | AUC² | PT ρ³ | scRNA DGE⁴ | Bulk rep.⁵ | "
          "ATAC⁶ | CellOracle⁷ | Drug score⁸ |")
sep    = "|---|:---:|:---:|:---:|---|---|:---:|:---:|:---:|"

footnotes = """
**Footnotes**
¹ pySCENIC 46-regulon GRN set (Y = included as master regulator).
² Mean AUCell score in LateAD-DAM state (SCENIC regulon members only).
³ Spearman ρ with diffusion pseudotime (HM→IRM→DAM→LAM→LateAD-DAM); regulon members only.
⁴ Best upregulated disease-vs-homeostatic pseudobulk DESeq2 contrast; [state] shown; * p<0.05, ** p<0.01, *** p<0.001.
⁵ GSE95587 fusiform gyrus bulk RNA-seq replication (n=117); * p<0.05, ** p<0.01, *** p<0.001.
⁶ Motif fold-enrichment in 56 AD-upregulated ATAC-seq DAPs (step 2.8).
⁷ Mean UMAP displacement magnitude in LateAD-DAM state under TF KO (CellOracle; step 5.7).
⁸ fpocket druggability drug_score for best-ranked structure (0–1; >0.5 = well-druggable).
† RUNX2: fpocket not determined — AF2 Q13580 is a 120-residue isoform with truncated Runt domain; therapeutic strategy is activation (Sun 2023), not inhibition.
"""

md_content = f"# Table 3 — Target Scorecard\n\n{header}\n{sep}\n" + "\n".join(md_rows) + "\n" + footnotes
OUT_MD.write_text(md_content)
print(f"\nSaved Markdown → {OUT_MD}")
