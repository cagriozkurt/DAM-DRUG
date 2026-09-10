"""
Section 2D.1 — Epistemic evidence-summary table for IKZF1
========================================================
Regenerates the Figure 3D / Figure 4 evidence table so that internal SEA-AD
modalities render as ONE non-independent block and external / robustness lines
are separated. Incorporates the Section 2 robustness results.

Output:
  results/phase7/epistemic/evidence_summary.csv
  results/phase7/epistemic/evidence_summary.png
"""
import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
OUT = PROJECT / "results/phase7/epistemic"
OUT.mkdir(parents=True, exist_ok=True)

ROWS = [
    # block, line, modality, source, result, epistemic_status
    ("Internal SEA-AD RNA (non-independent; same 84 donors)", 1,
     "Pseudobulk DGE", "SEA-AD snRNA", "IKZF1 up in DAM & IRM (padj<0.05)", "exploratory"),
    ("Internal SEA-AD RNA (non-independent; same 84 donors)", 2,
     "RcisTarget regulon", "SEA-AD snRNA", "IKZF1(+) peaks in LateAD-DAM (AUCell 0.153)", "exploratory"),
    ("Internal SEA-AD RNA (non-independent; same 84 donors)", 3,
     "Diffusion pseudotime", "SEA-AD snRNA", "only positive rho among 46 regulons (+0.309)", "exploratory"),
    ("Internal SEA-AD RNA (non-independent; same 84 donors)", 4,
     "CellOracle KO", "SEA-AD snRNA", "largest LateAD-DAM UMAP shift (ordinal)", "exploratory"),
    ("Internal SEA-AD chromatin (same donors, different assay)", 5,
     "ATAC motif enrichment", "SEA-AD ATAC", "IKZF1 motif 3.90x in AD-up peaks (descriptive)", "orthogonal-internal"),
    ("Robustness of the internal signal (Section 2)", 6,
     "Regional LMM", "SEA-AD snRNA, 10 regions", "DAM +0.57 / LateAD-DAM +0.55 log1p-CPM survive (1|region); 11/11 leave-one-region-out", "confirmatory-internal"),
    ("Robustness of the internal signal (Section 2)", 7,
     "Donor LOO + downsampling", "SEA-AD snRNA", "LateAD-DAM peak & rho stable to any single donor; gap>0 in 1000/1000 draws at 1 cell/donor", "confirmatory-internal"),
    ("External replication", 8,
     "Bulk RNA-seq", "GSE95587 (n=117)", "IKZF1 up in AD fusiform gyrus (log2FC +0.64, padj .004)", "external"),
    ("External replication", 9,
     "snRNA signature scoring", "Grubman 2019 GSE138852 (EC)", "IKZF1(+) signature reactive vs homeostatic microglia (delta +0.72, p 4e-17, perm p .001); NOT AD-vs-ct", "external-state-level"),
    ("Literature triangulation", 10,
     "Western blot", "Ballasch 2023 [30]", "IKZF1 protein up in AD hippocampus", "literature"),
]


def main():
    df = pd.DataFrame(ROWS, columns=["block", "line", "modality", "source",
                                     "result", "epistemic_status"])
    df.to_csv(OUT / "evidence_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(13, 4.2))
    ax.axis("off")
    blocks = df["block"].unique()
    palette = {b: c for b, c in zip(blocks, ["#DDE6F0", "#E4EFE0", "#CFE8DE",
                                             "#F5E6CC", "#EFE0EC"])}
    tbl = ax.table(
        cellText=df[["line", "modality", "source", "result", "epistemic_status"]].values,
        colLabels=["#", "Modality", "Source", "Result", "Epistemic status"],
        cellLoc="left", colWidths=[0.03, 0.14, 0.17, 0.46, 0.15], loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#333333"); cell.set_text_props(color="w", weight="bold")
        else:
            cell.set_facecolor(palette[df.iloc[r - 1]["block"]])
    # block labels down the left
    y = 0.86
    ax.set_title("IKZF1 evidence — internal SEA-AD analyses are one non-independent block",
                 fontsize=11, weight="bold")
    fig.text(0.01, 0.5, "\n".join(f"- {b}" for b in blocks), fontsize=7, va="center")
    fig.tight_layout()
    fig.savefig(OUT / "evidence_summary.png", dpi=200, bbox_inches="tight")
    print(df.to_string(index=False))
    print(f"\nWrote {OUT}/evidence_summary.csv and .png")


if __name__ == "__main__":
    main()
