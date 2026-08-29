# Regenerate CellChat plots from existing cellchat_object.rds
# Run on TRUBA (needs cellchat_r env) or locally if CellChat is installed.
# Skips the ~55-min computeCommunProb step.
#
# Usage:
#   apptainer exec containers/cellchat.sif Rscript --vanilla code/phase2_LR/03_cellchat_plots.R

suppressPackageStartupMessages({
    library(CellChat)
    library(ggplot2)
    library(patchwork)
})

proj_dir <- Sys.getenv("DAM_DRUG_DIR", unset=getwd())
cc_out   <- file.path(proj_dir, "results/phase2/LR/cellchat")

cat("Loading cellchat_object.rds...\n")
cc <- readRDS(file.path(cc_out, "cellchat_object.rds"))
keep_cts <- levels(cc@idents)
cat(sprintf("  %d cell types, %d pathways\n",
            length(keep_cts), length(cc@netP$pathways)))

# ── Bubble plot: L-R pairs → Microglia ───────────────────────────────────
cat("Generating bubble plot...\n")
tryCatch({
    pdf(file.path(cc_out, "bubble_to_microglia.pdf"), width=14, height=10)
    p <- netVisual_bubble(cc,
                          sources.use = keep_cts[keep_cts != "Microglia"],
                          targets.use = "Microglia",
                          remove.isolate = TRUE)
    print(p)
    dev.off()
    cat("  Saved bubble_to_microglia.pdf\n")
}, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat("  bubble failed:", conditionMessage(e), "\n")
})

# ── Heatmap: interaction strength ────────────────────────────────────────
cat("Generating heatmap...\n")
tryCatch({
    pdf(file.path(cc_out, "heatmap_strength.pdf"), width=8, height=7)
    ht <- netVisual_heatmap(cc, color.heatmap="Reds")
    ComplexHeatmap::draw(ht)
    dev.off()
    cat("  Saved heatmap_strength.pdf\n")
}, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat("  heatmap failed:", conditionMessage(e), "\n")
})

# ── Chord diagram ─────────────────────────────────────────────────────────
# Use netVisual_chord_cell (circlize) — draws sector names OUTSIDE the ring so
# they never overlap the arcs. Wide canvas + big outer margin so the outside
# labels are not clipped. Falls back to netVisual_circle if unavailable.
cat("Generating chord diagram...\n")
tryCatch({
    pdf(file.path(cc_out, "chord_count.pdf"), width=10, height=10)
    par(mar = c(6, 6, 6, 6), xpd = TRUE)
    netVisual_chord_cell(cc,
                         net = cc@net$count,
                         lab.cex = 1.1,
                         title.name = "Number of interactions",
                         show.legend = FALSE)
    dev.off()
    cat("  Saved chord_count.pdf (chord_cell, outside labels)\n")
}, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat("  chord_cell failed:", conditionMessage(e), "- falling back to circle\n")
    tryCatch({
        pdf(file.path(cc_out, "chord_count.pdf"), width=10, height=10)
        par(mar = c(2, 2, 3, 2))
        netVisual_circle(cc@net$count,
                         vertex.weight = rowSums(cc@net$count),
                         weight.scale = TRUE, label.edge = FALSE,
                         vertex.label.cex = 1.2,
                         title.name = "Number of interactions")
        dev.off()
        cat("  Saved chord_count.pdf (circle fallback)\n")
    }, error=function(e2) {
        try(dev.off(), silent=TRUE)
        cat("  chord failed:", conditionMessage(e2), "\n")
    })
})

cat("\nDone.\n")
