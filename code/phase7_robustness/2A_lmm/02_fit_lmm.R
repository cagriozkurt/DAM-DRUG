#!/usr/bin/env Rscript
# Section 2A.2 - Linear mixed-effects models for regional confounding
# ==================================================================
# For each target gene, on log1p(CPM) pseudobulk:
#   full     : expr ~ state + age_z + sex + (1|donor) + (1|region)
#   reduced  : expr ~ state + age_z + sex + (1|donor)
#   no_state : expr ~ age_z + sex + (1|donor) + (1|region)
# Fit under two data variants: 'min10' (>=10 cells/group) and 'unconstrained'.
#
# Outputs (results/phase7/lmm/):
#   <gene>_lmm.csv          fixed-effect state coefficients + CI + p
#   lmm_model_comparison.csv LRT full-vs-reduced, full-vs-no_state; ICC; R2; singular flags
#
# Usage: Rscript 02_fit_lmm.R

suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(broom.mixed)
})

proj <- Sys.getenv("DAM_DRUG_DIR", unset = getwd())
lmm_dir <- file.path(proj, "results/phase7/lmm")

long <- read.csv(file.path(lmm_dir, "pseudobulk_multiregion.csv"),
                 stringsAsFactors = FALSE)
meta <- read.csv(file.path(lmm_dir, "pseudobulk_group_meta.csv"),
                 stringsAsFactors = FALSE)

meta_cols <- c("group_id", "age_z", "sex", "braak", "adnc", "cognitive_status")
long <- merge(long, meta[, meta_cols], by = "group_id", all.x = TRUE)

genes  <- sort(unique(long$gene))
target_genes <- setdiff(genes, c("ACTB", "GAPDH", "B2M"))  # housekeeping = controls, still fit
variants <- c("min10", "unconstrained")

icc <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  tot <- sum(vc$vcov)
  setNames(vc$vcov / tot, ifelse(is.na(vc$var1), vc$grp, paste0(vc$grp)))
}

coef_rows <- list(); cmp_rows <- list()

for (g in genes) {
  for (v in variants) {
    d <- long[long$gene == g, ]
    if (v == "min10") d <- d[d$passes_min10 == "True" | d$passes_min10 == TRUE, ]
    d <- d[is.finite(d$log1p_cpm) & !is.na(d$age_z) & !is.na(d$sex) & !is.na(d$state), ]
    d$state  <- relevel(factor(d$state), ref = "Homeostatic")
    d$sex    <- factor(d$sex)
    d$donor  <- factor(d$donor)
    d$region <- factor(d$region)
    n_states <- nlevels(droplevels(d$state))
    if (nrow(d) < 30 || n_states < 2 || nlevels(droplevels(d$region)) < 2) {
      cmp_rows[[length(cmp_rows) + 1]] <- data.frame(
        gene = g, variant = v, status = "skipped_insufficient_data",
        n_obs = nrow(d), n_states = n_states)
      next
    }

    ok <- TRUE
    full <- tryCatch(
      lmer(log1p_cpm ~ state + age_z + sex + (1 | donor) + (1 | region),
           data = d, REML = FALSE,
           control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))),
      error = function(e) { ok <<- FALSE; NULL })
    reduced <- tryCatch(
      lmer(log1p_cpm ~ state + age_z + sex + (1 | donor),
           data = d, REML = FALSE), error = function(e) NULL)
    nostate <- tryCatch(
      lmer(log1p_cpm ~ age_z + sex + (1 | donor) + (1 | region),
           data = d, REML = FALSE), error = function(e) NULL)
    if (!ok || is.null(full)) {
      cmp_rows[[length(cmp_rows) + 1]] <- data.frame(
        gene = g, variant = v, status = "full_model_failed",
        n_obs = nrow(d), n_states = n_states)
      next
    }

    ci <- tryCatch(confint(full, method = "Wald", parm = "beta_"),
                   error = function(e) NULL)
    tt <- as.data.frame(coef(summary(full)))
    tt$term <- rownames(tt)
    for (trm in grep("^state", tt$term, value = TRUE)) {
      lo <- hi <- NA
      if (!is.null(ci) && trm %in% rownames(ci)) { lo <- ci[trm, 1]; hi <- ci[trm, 2] }
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        gene = g, variant = v, term = trm,
        estimate = tt[trm, "Estimate"], se = tt[trm, "Std. Error"],
        df = tt[trm, "df"], t = tt[trm, "t value"],
        p_value = tt[trm, "Pr(>|t|)"], ci_lo = lo, ci_hi = hi,
        n_obs = nrow(d))
    }

    lrt_region <- if (!is.null(reduced)) anova(reduced, full) else NULL
    lrt_state  <- if (!is.null(nostate)) anova(nostate, full) else NULL
    ic <- icc(full)
    r2 <- tryCatch(suppressWarnings(as.data.frame(broom.mixed::glance(full))),
                   error = function(e) data.frame())

    cmp_rows[[length(cmp_rows) + 1]] <- data.frame(
      gene = g, variant = v, status = "ok", n_obs = nrow(d), n_states = n_states,
      is_singular = isSingular(full),
      aic_full = AIC(full), aic_reduced = if (!is.null(reduced)) AIC(reduced) else NA,
      lrt_region_chisq = if (!is.null(lrt_region)) lrt_region$Chisq[2] else NA,
      lrt_region_p     = if (!is.null(lrt_region)) lrt_region$`Pr(>Chisq)`[2] else NA,
      lrt_state_chisq  = if (!is.null(lrt_state)) lrt_state$Chisq[2] else NA,
      lrt_state_p      = if (!is.null(lrt_state)) lrt_state$`Pr(>Chisq)`[2] else NA,
      icc_donor  = ifelse("donor"  %in% names(ic), ic[["donor"]], NA),
      icc_region = ifelse("region" %in% names(ic), ic[["region"]], NA),
      var_resid  = attr(VarCorr(full), "sc")^2)
  }
}

coef_df <- do.call(rbind, coef_rows)
cmp_df  <- do.call(rbind, cmp_rows)

write.csv(cmp_df, file.path(lmm_dir, "lmm_model_comparison.csv"), row.names = FALSE)
for (g in unique(coef_df$gene)) {
  write.csv(coef_df[coef_df$gene == g, ],
            file.path(lmm_dir, paste0(g, "_lmm.csv")), row.names = FALSE)
}

cat("\n=== state coefficients (full model, min10 variant) ===\n")
sub <- coef_df[coef_df$variant == "min10", ]
print(sub[, c("gene", "term", "estimate", "se", "ci_lo", "ci_hi", "p_value")],
      row.names = FALSE, digits = 3)
cat("\n=== model comparison (does (1|region) matter? does state survive?) ===\n")
print(cmp_df[cmp_df$status == "ok",
             c("gene", "variant", "is_singular", "lrt_region_p", "lrt_state_p",
               "icc_donor", "icc_region")], row.names = FALSE, digits = 3)
cat("\nWrote", file.path(lmm_dir, "lmm_model_comparison.csv"), "and per-gene _lmm.csv\n")
