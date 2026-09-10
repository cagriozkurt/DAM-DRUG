#!/usr/bin/env Rscript
# WP2.3 — Per-region CellChat v2. Parameter-identical to the accepted MTG run
# (code/phase2_LR/02_cellchat_nichechat.R): CellChatDB.human, "Secreted
# Signaling" subset, identifyOverExpressedGenes(do.fast=FALSE),
# computeCommunProb(type="triMean"), filterCommunication(min.cells=10).
#
# Usage (inside dam-drug-r.sif):
#   Rscript s3_cellchat.R <REGION>
# Reads  results/phase7/cellchat_multiregion/<REGION>/prep/
# Writes results/phase7/cellchat_multiregion/<REGION>/
#   cellchat_<REGION>.rds, netP_<REGION>.csv, lr_interactions_<REGION>.csv

suppressPackageStartupMessages({
  library(CellChat); library(Matrix); library(hdf5r)
})

args   <- commandArgs(trailingOnly = TRUE)
region <- args[1]
proj   <- Sys.getenv("DAM_DRUG_DIR", unset = getwd())
pdir   <- file.path(proj, "results/phase7/cellchat_multiregion", region)
prep   <- file.path(pdir, "prep")
dir.create(pdir, showWarnings = FALSE, recursive = TRUE)

## ── load counts (CSC h5 -> dgCMatrix) ──────────────────────────────────────
h5 <- H5File$new(file.path(prep, "counts_raw.h5"), mode = "r")
shape <- h5[[".."]]$attr_open("shape")$read()   # (ncells, ngenes)
data_v <- h5[["data"]]$read(); idx_v <- h5[["indices"]]$read(); ptr_v <- h5[["indptr"]]$read()
barcodes <- h5[["barcodes"]]$read(); genes <- h5[["gene_names"]]$read()
h5$close_all()
# CSC over cells x genes: p = column(gene) pointers, i = row(cell) indices
M <- sparseMatrix(i = idx_v + 1, p = ptr_v, x = data_v,
                  dims = c(shape[1], shape[2]),
                  dimnames = list(barcodes, genes))
counts <- t(M)                     # genes x cells for CellChat

meta <- read.csv(file.path(prep, "cell_meta.csv"), row.names = 1)
meta <- meta[colnames(counts), , drop = FALSE]
group <- factor(meta$cell_type_broad)
cat(sprintf("[%s] %d genes x %d cells; groups: %s\n", region,
            nrow(counts), ncol(counts), paste(levels(group), collapse = ", ")))

## ── CellChat (identical params to MTG) ─────────────────────────────────────
cc <- createCellChat(object = counts, meta = data.frame(group = group,
                     row.names = colnames(counts)), group.by = "group")
cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling",
                  key = "annotation")
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
cc <- identifyOverExpressedInteractions(cc)
options(future.globals.maxSize = 8 * 1024^3)
cc <- computeCommunProb(cc, type = "triMean")
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)

saveRDS(cc, file.path(pdir, sprintf("cellchat_%s.rds", region)))

## ── export interaction tables ─────────────────────────────────────────────
lr <- subsetCommunication(cc)                       # per L-R, per sender-receiver
lr$region <- region
write.csv(lr, file.path(pdir, sprintf("lr_interactions_%s.csv", region)),
          row.names = FALSE)

netp <- subsetCommunication(cc, slot.name = "netP") # per pathway
netp$region <- region
write.csv(netp, file.path(pdir, sprintf("netP_%s.csv", region)), row.names = FALSE)

## SLIT2 -> ROBO2 quick look
s <- lr[lr$ligand == "SLIT2" & lr$receptor == "ROBO2", ]
if (nrow(s)) {
  s <- s[order(-s$prob), ]
  cat(sprintf("[%s] SLIT2->ROBO2 top: %s -> %s  prob=%.4g\n",
              region, s$source[1], s$target[1], s$prob[1]))
} else {
  cat(sprintf("[%s] SLIT2->ROBO2 not detected\n", region))
}
cat(sprintf("[%s] done\n", region))
