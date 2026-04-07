# Table 3 — Target Scorecard

| TF | GRN¹ | AUC² | PT ρ³ | scRNA DGE⁴ | Bulk rep.⁵ | ATAC⁶ | CellOracle⁷ | Drug score⁸ |
|---|:---:|:---:|:---:|---|---|:---:|:---:|:---:|
| PPARG | N | — | — | +1.15 [LAM]*** | -0.021 (ns) | 4.5× | — | 0.996 |
| IKZF1 | Y | 0.1526 | +0.309 | +0.75 [DAM]*** | +0.638** | 3.9× | 0.047 | 0.001 |
| IRF8 | N | — | — | +1.70 [DAM]*** | +0.615* | — | 0.026 | 0.659 |
| BHLHE41 | N | — | — | +2.55 [DAM]*** | +0.268* | — | — | 0.003 |
| SPI1 | N | — | — | +0.12 [IRM]* | +0.489* | — | 0.013 | 0.381 |
| RUNX1 | N | — | — | +0.43 [DAM]*** | +0.399 (ns) | — | — | 0.092 |
| RUNX2 | N | — | — | — | — | — | — | n.d.† |
| MAF | N | — | — | — | +0.220 (ns) | — | — | 0.432 |
| ACSL1 | N | — | — | — | — | — | — | 0.000 |
| PIK3CA | N | — | — | — | — | — | — | 0.012 |

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
