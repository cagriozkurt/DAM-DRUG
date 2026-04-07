# DAM-DRUG | PICO 4 — Meta-Analysis: MM-GBSA Accuracy vs. Experimental Affinity
# Pools Spearman ρ (MM-GBSA ΔG_bind vs. experimental Kexp) across studies
# Output: forest plot + subgroup analysis by target class and MD duration
# Requires: metafor, readxl

library(metafor)
library(readxl)

# ── Load extraction table ──────────────────────────────────────────────────
dat <- read_excel("../evidence_tables/PICO4_MD_freeenergy.xlsx",
                  sheet = "extraction")
# Expected columns:
#   study_id | target_class | md_duration_ns | spearman_rho | n_compounds
#   | fep_rho (if reported) | GB_model | force_field_protein

# ── Fisher z-transformation of Spearman ρ ─────────────────────────────────
dat$zi  <- 0.5 * log((1 + dat$spearman_rho) / (1 - dat$spearman_rho))
dat$sei <- 1 / sqrt(dat$n_compounds - 3)

# ── Overall pooled estimate ────────────────────────────────────────────────
res_overall <- rma(yi = zi, sei = sei, data = dat, method = "REML",
                   slab = study_id)

# Back-transform to ρ scale
pooled_rho    <- tanh(res_overall$beta)
pooled_rho_lb <- tanh(res_overall$ci.lb)
pooled_rho_ub <- tanh(res_overall$ci.ub)
cat(sprintf("Pooled Spearman ρ (MM-GBSA): %.3f [%.3f, %.3f]  I²=%.1f%%\n",
            pooled_rho, pooled_rho_lb, pooled_rho_ub, res_overall$I2))

# ── Forest plot ────────────────────────────────────────────────────────────
pdf("forest_MMGBSA_rho.pdf", width = 10, height = 7)
forest(res_overall,
       atransf = tanh,
       at      = atanh(c(0.3, 0.5, 0.7, 0.9)),
       main    = "MM-GBSA Spearman ρ vs. Experiment (random-effects)")
dev.off()

# ── Subgroup: target class ─────────────────────────────────────────────────
# Classes: kinase | TF | GPCR | nuclear_receptor | other
res_subgroup_class <- rma(yi = zi, sei = sei, data = dat,
                           mods = ~ target_class, method = "REML")
cat("\n--- Subgroup by target class ---\n")
print(summary(res_subgroup_class))

# ── Subgroup: MD duration ──────────────────────────────────────────────────
dat$duration_bin <- cut(dat$md_duration_ns,
                        breaks = c(0, 100, 500, Inf),
                        labels = c("<=100ns", "101-500ns", ">500ns"))
res_subgroup_dur <- rma(yi = zi, sei = sei, data = dat,
                         mods = ~ duration_bin, method = "REML")
cat("\n--- Subgroup by MD duration ---\n")
print(summary(res_subgroup_dur))

# ── Plan revision trigger check ────────────────────────────────────────────
# Trigger: pooled ρ < 0.5 → flag for FEP escalation
if (pooled_rho < 0.5) {
  warning("TRIGGER: Pooled MM-GBSA rho < 0.5. Flag for FEP escalation on top 5 compounds.")
} else {
  cat("MM-GBSA accuracy acceptable (rho >= 0.5). No plan revision triggered.\n")
}

# ── Save results ───────────────────────────────────────────────────────────
write.csv(
  data.frame(pooled_rho, pooled_rho_lb, pooled_rho_ub,
             I2 = res_overall$I2, k = res_overall$k),
  "pooled_MMGBSA_rho.csv", row.names = FALSE
)
cat("Done. Results saved to pooled_MMGBSA_rho.csv\n")