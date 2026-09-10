#!/usr/bin/env Rscript
# Section 2A.3 - LMM sensitivity analyses
# ======================================
# (a) leave-one-region-out: refit the full model dropping each region once;
#     does the DAM / LateAD-DAM state coefficient stay positive, CI excluding 0?
# (b) region as a FIXED effect instead of random (robustness to few region levels).
# (c) NB-GLMM on raw counts with log(total_umi) offset (IKZF1 only, best effort).
#
# Outputs (results/phase7/lmm/):
#   lmm_leave_one_region_out.csv
#   lmm_region_fixed_effect.csv
#   lmm_nbglmm_ikzf1.csv
#
# Usage: Rscript 04_lmm_sensitivity.R

suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

proj <- Sys.getenv("DAM_DRUG_DIR", unset = getwd())
lmm_dir <- file.path(proj, "results/phase7/lmm")

long <- read.csv(file.path(lmm_dir, "pseudobulk_multiregion.csv"), stringsAsFactors = FALSE)
meta <- read.csv(file.path(lmm_dir, "pseudobulk_group_meta.csv"), stringsAsFactors = FALSE)
long <- merge(long, meta[, c("group_id", "age_z", "sex")], by = "group_id", all.x = TRUE)
long$passes_min10 <- long$passes_min10 %in% c("True", TRUE, "TRUE")

GENES <- c("IKZF1", "IRF8", "SPI1", "BHLHE41", "RUNX1", "PPARG", "CEBPB")
FOCUS <- c("stateDAM", "stateLateAD-DAM")

prep <- function(g, min10 = TRUE) {
  d <- long[long$gene == g, ]
  if (min10) d <- d[d$passes_min10, ]
  d <- d[is.finite(d$log1p_cpm) & !is.na(d$age_z) & !is.na(d$sex), ]
  d$state  <- relevel(factor(d$state), ref = "Homeostatic")
  d$sex <- factor(d$sex); d$donor <- factor(d$donor); d$region <- factor(d$region)
  d
}

grab <- function(m, terms) {
  tt <- as.data.frame(coef(summary(m))); tt$term <- rownames(tt)
  ci <- tryCatch(confint(m, method = "Wald", parm = terms), error = function(e) NULL)
  out <- list()
  for (trm in terms) {
    if (!trm %in% tt$term) next
    lo <- hi <- NA
    if (!is.null(ci) && trm %in% rownames(ci)) { lo <- ci[trm, 1]; hi <- ci[trm, 2] }
    out[[trm]] <- data.frame(term = trm, estimate = tt[tt$term == trm, "Estimate"],
                             se = tt[tt$term == trm, "Std. Error"], ci_lo = lo, ci_hi = hi)
  }
  do.call(rbind, out)
}

## (a) leave-one-region-out ------------------------------------------------
loro <- list()
for (g in GENES) {
  d0 <- prep(g)
  regions <- levels(droplevels(d0$region))
  for (rr in c("<none>", regions)) {
    d <- if (rr == "<none>") d0 else droplevels(d0[d0$region != rr, ])
    if (nlevels(droplevels(d$region)) < 2) next
    m <- tryCatch(lmer(log1p_cpm ~ state + age_z + sex + (1 | donor) + (1 | region),
                       data = d, REML = FALSE), error = function(e) NULL)
    if (is.null(m)) next
    cf <- grab(m, FOCUS)
    if (is.null(cf)) next
    cf$gene <- g; cf$region_dropped <- rr; cf$n_obs <- nrow(d)
    cf$is_singular <- isSingular(m)
    loro[[length(loro) + 1]] <- cf
  }
}
loro_df <- do.call(rbind, loro)
loro_df$sig_positive <- loro_df$ci_lo > 0
write.csv(loro_df, file.path(lmm_dir, "lmm_leave_one_region_out.csv"), row.names = FALSE)

cat("=== leave-one-region-out: fraction of refits with CI_lo > 0 ===\n")
agg <- aggregate(sig_positive ~ gene + term, loro_df, function(x) sprintf("%d/%d", sum(x), length(x)))
print(agg, row.names = FALSE)

## (b) region as fixed effect --------------------------------------------
rfx <- list()
for (g in GENES) {
  d <- prep(g)
  m <- tryCatch(lmer(log1p_cpm ~ state + age_z + sex + region + (1 | donor),
                     data = d, REML = FALSE), error = function(e) NULL)
  if (is.null(m)) next
  cf <- grab(m, FOCUS); cf$gene <- g; cf$is_singular <- isSingular(m)
  rfx[[length(rfx) + 1]] <- cf
}
rfx_df <- do.call(rbind, rfx)
write.csv(rfx_df, file.path(lmm_dir, "lmm_region_fixed_effect.csv"), row.names = FALSE)
cat("\n=== region as FIXED effect: DAM / LateAD-DAM state coefficients ===\n")
print(rfx_df[, c("gene", "term", "estimate", "ci_lo", "ci_hi")], row.names = FALSE, digits = 3)

## (c) NB-GLMM on raw counts, IKZF1 -------------------------------------
d <- prep("IKZF1")
grp_umi <- meta[, c("group_id", "total_umi")]
d <- merge(d, grp_umi, by = "group_id", all.x = TRUE)
d <- d[d$total_umi > 0 & is.finite(d$count), ]
d$off <- log(d$total_umi)
nb <- tryCatch(
  glmer.nb(count ~ state + age_z + sex + (1 | donor) + (1 | region) + offset(off),
           data = d, control = glmerControl(optimizer = "bobyqa")),
  error = function(e) { cat("NB-GLMM failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(nb)) {
  tt <- as.data.frame(coef(summary(nb))); tt$term <- rownames(tt)
  keep <- tt[grepl("^state", tt$term), c("term", "Estimate", "Std. Error", "Pr(>|z|)")]
  names(keep) <- c("term", "log_rate_ratio", "se", "p_value")
  write.csv(keep, file.path(lmm_dir, "lmm_nbglmm_ikzf1.csv"), row.names = FALSE)
  cat("\n=== NB-GLMM IKZF1 (log rate ratio vs Homeostatic) ===\n")
  print(keep, row.names = FALSE, digits = 3)
}
cat("\nDone.\n")
