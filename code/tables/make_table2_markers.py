"""
DAM-DRUG Table 2 — Top Marker Genes per Microglial State
=========================================================
Source: results/phase1/DGE_pseudobulk/dge_pb_*_vs_Homeostatic.csv
        (DESeq2 pseudobulk, one sample = one donor)

Outputs:
  results/tables/table2_marker_genes.csv
  results/tables/table2_marker_genes.md

Runs locally — no h5ad needed.
"""

import os
import pandas as pd
from pathlib import Path

PROJECT = Path(os.environ.get("DAM_DRUG_DIR",
               str(Path(__file__).resolve().parents[2])))
DGE_DIR = PROJECT / "results/phase1/DGE_pseudobulk"
OUT     = PROJECT / "results/tables"
OUT.mkdir(parents=True, exist_ok=True)

TOP_N = 5

# lncRNA / uncharacterised transcript name patterns (GENCODE/HGNC conventions)
import re
NONCODING_RE = re.compile(
    r'^(LINC\d|AC\d{6}\.|AL\d{6}\.|AP\d{6}\.|BX\d{6}\.|AJ\d{6}\.|'
    r'XACT$|NEAT1$|MALAT1$|MIR\d|SNHG\d|SNORD\d|RNU\d)'
)

def is_protein_coding(gene: str) -> bool:
    return NONCODING_RE.match(gene) is None

# ── State label map ───────────────────────────────────────────────────────────
CONTRASTS = {
    "dge_pb_IRM_vs_Homeostatic.csv":  "IRM",
    "dge_pb_DAM_vs_Homeostatic.csv":  "DAM",
    "dge_pb_LAM_vs_Homeostatic.csv":  "LAM",
}
# Display order (trajectory order)
STATE_ORDER = ["IRM", "DAM", "LAM"]

# ── Load, filter, rank ────────────────────────────────────────────────────────
frames = []
for fname, state in CONTRASTS.items():
    path = DGE_DIR / fname
    df = pd.read_csv(path)
    # Filter: upregulated, significant, protein-coding
    df = df[(df["log2FoldChange"] > 0) & (df["padj"] < 0.05)].copy()
    df = df[df["names"].apply(is_protein_coding)]
    df = df.sort_values("log2FoldChange", ascending=False).head(TOP_N)
    df["State"] = state
    frames.append(df)

all_markers = pd.concat(frames, ignore_index=True)

# ── Format output columns ─────────────────────────────────────────────────────
result = pd.DataFrame({
    "State":            all_markers["State"],
    "Gene":             all_markers["names"],
    "log2FC":           all_markers["log2FoldChange"].round(2),
    "padj":             all_markers["padj"].apply(
                            lambda p: "<0.001" if p < 0.001 else f"{p:.3f}"),
    "% expressed":      (all_markers["pts"] * 100).round(0).astype(int).astype(str) + "%",
    "baseMean":         all_markers["baseMean"].round(1),
})

# Sort by trajectory order
result["_ord"] = result["State"].map({s: i for i, s in enumerate(STATE_ORDER)})
result = result.sort_values(["_ord", "log2FC"], ascending=[True, False]).drop(columns="_ord")
result = result.reset_index(drop=True)

# ── Save CSV ──────────────────────────────────────────────────────────────────
csv_path = OUT / "table2_marker_genes.csv"
result.to_csv(csv_path, index=False)
print(f"Wrote {csv_path}")

# ── Save Markdown ─────────────────────────────────────────────────────────────
cols = list(result.columns)
header = "| " + " | ".join(cols) + " |"
sep    = "| " + " | ".join(["---"] * len(cols)) + " |"
body   = "\n".join(
    "| " + " | ".join(str(result.iloc[i][c]) for c in cols) + " |"
    for i in range(len(result))
)
md = (
    "# Table 2 — Top Marker Genes per Microglial State\n\n"
    + header + "\n" + sep + "\n" + body
    + "\n\nTop 5 upregulated genes per state vs Homeostatic microglia "
    "(DESeq2 pseudobulk; one pseudobulk sample per donor). "
    "Sorted by log2 fold-change. padj: Benjamini–Hochberg FDR. "
    "% expressed: fraction of nuclei in the state expressing the gene."
)
md_path = OUT / "table2_marker_genes.md"
md_path.write_text(md)
print(f"Wrote {md_path}")

print("\n--- Table preview ---")
print(result.to_string(index=False))
