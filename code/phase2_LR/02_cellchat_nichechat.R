# Phase 2.9–2.10: CellChat + NicheNet — L-R interactions in MTG microglia
# ==========================================================================
# Inputs  (from 01_prep_mtg_for_cellchat.py):
#   results/phase2/LR/prep/counts_raw.h5
#   results/phase2/LR/prep/cell_meta.csv
#   results/phase2/LR/prep/microglia_markers.csv
#
# Outputs (results/phase2/LR/):
#   cellchat/cellchat_object.rds
#   cellchat/bubble_*.pdf          — L-R bubble plots per sender
#   cellchat/chord_*.pdf           — Chord diagrams
#   cellchat/lr_summary.csv        — Top L-R pairs table
#   nichechat/nichenet_object.rds
#   nichechat/ligand_activity.csv  — Top ligands ranked by Pearson
#   nichechat/ligand_target_heatmap.pdf
#   nichechat/ligand_expression_dotplot.pdf
#
# Run: sbatch code/slurm/17_cellchat_final.slurm

suppressPackageStartupMessages({
    library(Matrix)
    library(hdf5r)
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(CellChat)
})
# nichenetr is optional — not installable without Seurat in this env
HAS_NICHENET <- requireNamespace("nichenetr", quietly=TRUE)
if (HAS_NICHENET) suppressPackageStartupMessages(library(nichenetr))
cat(sprintf("nichenetr available: %s\n", HAS_NICHENET))

# ── Paths ─────────────────────────────────────────────────────────────────────
proj_dir  <- Sys.getenv("DAM_DRUG_DIR",
                         unset=getwd())
prep_dir  <- file.path(proj_dir, "results/phase2/LR/prep")
cc_out    <- file.path(proj_dir, "results/phase2/LR/cellchat")
nn_out    <- file.path(proj_dir, "results/phase2/LR/nichechat")
dir.create(cc_out, recursive=TRUE, showWarnings=FALSE)
dir.create(nn_out, recursive=TRUE, showWarnings=FALSE)

cat("=== CellChat + NicheNet ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ── 1. Load preprocessed data ────────────────────────────────────────────────
cat("Loading preprocessed MTG data...\n")
h5f      <- H5File$new(file.path(prep_dir, "counts_raw.h5"), mode="r")
data_vec <- h5f[["data"]][]
indices  <- h5f[["indices"]][]
indptr   <- h5f[["indptr"]][]
dims     <- h5f$attr_open("shape")$read()
barcodes <- h5f[["barcodes"]][]
genes    <- h5f[["gene_names"]][]
h5f$close_all()

# Decode byte strings (Python h5py stores as bytes)
barcodes <- sub("b'(.+)'", "\\1", barcodes)
genes    <- sub("b'(.+)'", "\\1", genes)

# Matrix is CSC (cells × genes): indptr = column pointers (length ngenes+1),
# indices = row indices (cell indices). R sparseMatrix(p=, i=) expects exactly this.
counts_cc_t <- sparseMatrix(
    i = indices + 1L,   # 0-indexed → 1-indexed row indices (cells)
    p = indptr,         # column pointers, length ngenes+1
    x = as.numeric(data_vec),
    dims = dims,        # c(ncells, ngenes)
    dimnames = list(barcodes, genes),
    repr = "C"
)
counts <- t(counts_cc_t)   # genes × cells (CellChat convention)
cat(sprintf("  Counts: %d genes × %d cells\n", nrow(counts), ncol(counts)))

meta <- fread(file.path(prep_dir, "cell_meta.csv"), data.table=FALSE)
rownames(meta) <- meta[[1]]; meta[[1]] <- NULL
# Align to count matrix columns
meta <- meta[colnames(counts), , drop=FALSE]

cat("  Cell type composition (broad):\n")
print(table(meta$cell_type_broad))

# ── 2. CellChat ───────────────────────────────────────────────────────────────
cat("\n--- CellChat ---\n")

# Remove very small groups (< 10 cells)
ct_counts <- table(meta$cell_type_broad)
keep_cts  <- names(ct_counts)[ct_counts >= 10]
keep_cells <- meta$cell_type_broad %in% keep_cts
counts_cc <- counts[, keep_cells]
meta_cc   <- meta[keep_cells, ]
meta_cc$samples <- "sample1"   # required by CellChat; all cells from one merged object
cat(sprintf("  Kept %d cells across %d cell types\n",
            sum(keep_cells), length(keep_cts)))

cc <- createCellChat(object=counts_cc,
                     meta=meta_cc,
                     group.by="cell_type_broad")
CellChatDB <- CellChatDB.human   # built-in human L-R database
# Use the secreted signalling subset (most relevant for microglia)
CellChatDB.use <- subsetDB(CellChatDB, search="Secreted Signaling",
                            key="annotation")
cc@DB <- CellChatDB.use

cc <- subsetData(cc)
# identifyOverExpressedGenes runs sequentially (do.fast=FALSE; parallel adds no benefit)
cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
cc <- identifyOverExpressedInteractions(cc)
# computeCommunProb benefits from parallelism; raise globals limit for large matrix
options(future.globals.maxSize = 8 * 1024^3)   # 8 GiB
future::plan("multisession", workers=4)
cc <- computeCommunProb(cc, type="triMean")
future::plan("sequential")
cc <- filterCommunication(cc, min.cells=10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)

cat("  Saving CellChat object...\n")
saveRDS(cc, file.path(cc_out, "cellchat_object.rds"))

# ── CellChat plots ────────────────────────────────────────────────────────────
cat("  Generating CellChat plots...\n")

# Bubble plot: all L-R pairs where Microglia is receiver
tryCatch({
    pdf(file.path(cc_out, "bubble_to_microglia.pdf"), width=12, height=8)
    p <- netVisual_bubble(cc, sources.use=keep_cts[keep_cts != "Microglia"],
                          targets.use="Microglia",
                          remove.isolate=TRUE)
    print(p)
    dev.off()
}, error=function(e) cat("  bubble plot failed:", conditionMessage(e), "\n"))

# Chord diagram: number of interactions per cell type pair
tryCatch({
    pdf(file.path(cc_out, "chord_count.pdf"), width=8, height=8)
    netVisual_circle(cc@net$count, vertex.weight=rowSums(cc@net$count),
                     weight.scale=TRUE, label.edge=FALSE,
                     title.name="Number of interactions")
    dev.off()
}, error=function(e) cat("  chord plot failed:", conditionMessage(e), "\n"))

# Heatmap of interaction strength
tryCatch({
    pdf(file.path(cc_out, "heatmap_strength.pdf"), width=8, height=6)
    netVisual_heatmap(cc, color.heatmap="Reds")
    dev.off()
}, error=function(e) cat("  heatmap failed:", conditionMessage(e), "\n"))

# Export top L-R pairs table
lr_df <- subsetCommunication(cc)
write.csv(lr_df, file.path(cc_out, "lr_summary.csv"), row.names=FALSE)
cat(sprintf("  Wrote lr_summary.csv (%d L-R pairs)\n", nrow(lr_df)))

# ── 3. NicheNet ───────────────────────────────────────────────────────────────
cat("\n--- NicheNet ---\n")
if (!HAS_NICHENET) {
    cat("  nichenetr not installed — skipping NicheNet analysis.\n")
    cat("  CellChat outputs (lr_summary.csv, bubble, chord) are complete above.\n")
    cat("\n=== Done:", format(Sys.time()), "===\n")
    quit(save="no", status=0)
}

# Load NicheNet prior model (downloads to tempdir if not cached)
nn_model_dir <- file.path(proj_dir, "data/resources/nichenet")
dir.create(nn_model_dir, recursive=TRUE, showWarnings=FALSE)

get_nn_file <- function(fname) {
    local_path <- file.path(nn_model_dir, fname)
    if (!file.exists(local_path)) {
        url <- paste0("https://zenodo.org/record/7074291/files/", fname)
        cat(sprintf("  Downloading %s...\n", fname))
        download.file(url, local_path, quiet=TRUE)
    }
    readRDS(local_path)
}

ligand_target_matrix    <- get_nn_file("ligand_target_matrix_nsga2r_final.rds")
lr_network              <- get_nn_file("lr_network_human_21122021.rds")
weighted_networks       <- get_nn_file("weighted_networks_nsga2r_final.rds")

# Receiver = Microglia; senders = all other types
receiver_ct  <- "Microglia"
sender_cts   <- keep_cts[keep_cts != receiver_ct]

# Expressed genes per cell type (≥ 10% cells)
expressed_genes <- function(cts, min_pct=0.10) {
    cell_idx <- which(meta_cc$cell_type_broad %in% cts)
    mat      <- counts_cc[, cell_idx, drop=FALSE]
    genes_expressed <- rownames(mat)[rowMeans(mat > 0) >= min_pct]
    return(genes_expressed)
}
receiver_expressed <- expressed_genes(receiver_ct)
sender_expressed   <- expressed_genes(sender_cts)
cat(sprintf("  Receiver expressed: %d genes\n", length(receiver_expressed)))
cat(sprintf("  Sender expressed:   %d genes\n", length(sender_expressed)))

# DAM geneset (from preprocessing)
marker_file <- file.path(prep_dir, "microglia_markers.csv")
if (file.exists(marker_file)) {
    dam_markers <- read.csv(marker_file)$gene
} else {
    # Fallback hardcoded DAM signature
    dam_markers <- c("TREM2", "APOE", "SPP1", "LGALS3", "LPL", "GPNMB",
                     "IKZF1", "IRF8", "PPARG", "SPI1", "RUNX1",
                     "CST7", "CTSL", "CTSD", "AXL", "CD9")
}
geneset_oi <- intersect(dam_markers, receiver_expressed)
cat(sprintf("  DAM geneset (intersection with expressed): %d genes\n",
            length(geneset_oi)))

# Background = all expressed receiver genes
background_genes <- receiver_expressed

# Potential ligands (sender-expressed, have receptor in receiver)
all_ligands      <- unique(lr_network$from)
all_receptors    <- unique(lr_network$to)
pot_ligands      <- intersect(sender_expressed, all_ligands)
pot_receptors    <- intersect(receiver_expressed, all_receptors)
lr_pairs         <- lr_network %>%
    filter(from %in% pot_ligands & to %in% pot_receptors)
pot_ligands_filt <- unique(lr_pairs$from)
cat(sprintf("  Potential ligands: %d\n", length(pot_ligands_filt)))

# Ligand activity analysis
if (length(geneset_oi) >= 5 && length(pot_ligands_filt) >= 5) {
    ligand_activities <- predict_ligand_activities(
        geneset=geneset_oi,
        background_expressed_genes=background_genes,
        ligand_target_matrix=ligand_target_matrix,
        potential_ligands=pot_ligands_filt
    )
    ligand_activities <- ligand_activities %>%
        arrange(desc(aupr_corrected)) %>%
        mutate(rank=row_number())

    write.csv(ligand_activities, file.path(nn_out, "ligand_activity.csv"),
              row.names=FALSE)
    cat(sprintf("  Wrote ligand_activity.csv (%d ligands ranked)\n",
                nrow(ligand_activities)))

    # Top 20 ligands
    top_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>%
        pull(test_ligand)

    # Ligand–target heatmap
    tryCatch({
        active_ligand_target <- ligand_target_matrix[geneset_oi, top_ligands] %>%
            .[rowSums(.) > 0, ]
        pdf(file.path(nn_out, "ligand_target_heatmap.pdf"), width=12, height=8)
        p <- make_heatmap_ggplot(
            t(active_ligand_target),
            y_name="Ligand", x_name="Target gene",
            color=colorRampPalette(c("white", "#D55E00"))(100),
            legend_title="Regulatory potential"
        )
        print(p)
        dev.off()
    }, error=function(e) cat("  heatmap failed:", conditionMessage(e), "\n"))

    # Ligand expression dotplot across sender cell types
    tryCatch({
        # Mean expression per sender type
        expr_df <- lapply(sender_cts, function(ct) {
            cell_idx <- which(meta_cc$cell_type_broad == ct)
            mat      <- counts_cc[top_ligands, cell_idx, drop=FALSE]
            data.frame(
                ligand    = top_ligands,
                cell_type = ct,
                mean_expr = rowMeans(mat),
                pct_expr  = rowMeans(mat > 0) * 100
            )
        }) %>% bind_rows()

        p <- ggplot(expr_df, aes(x=cell_type, y=ligand,
                                  size=pct_expr, color=mean_expr)) +
            geom_point() +
            scale_color_gradient(low="white", high="#D55E00",
                                 name="Mean expr") +
            scale_size_continuous(name="% cells", range=c(1, 6)) +
            theme_bw(base_size=9) +
            theme(axis.text.x=element_text(angle=35, hjust=1)) +
            labs(x="Sender cell type", y="Ligand",
                 title="Top NicheNet ligand expression in sender cells")
        ggsave(file.path(nn_out, "ligand_expression_dotplot.pdf"),
               p, width=10, height=8)
    }, error=function(e) cat("  dotplot failed:", conditionMessage(e), "\n"))

    saveRDS(list(ligand_activities=ligand_activities,
                 top_ligands=top_ligands,
                 lr_pairs=lr_pairs),
            file.path(nn_out, "nichenet_results.rds"))
} else {
    cat("  WARNING: too few genes for NicheNet — skipping activity analysis\n")
}

cat("\n=== Done:", format(Sys.time()), "===\n")
