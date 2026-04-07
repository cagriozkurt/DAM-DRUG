# DAM-DRUG | PICO 2a — Meta-Analysis: scRNA-seq Integration Benchmarking
# Pools kBET acceptance rates across studies comparing scVI vs. Harmony
# Output: forest plots + pooled estimates (random-effects model)
# Requires: metafor, ggplot2, readxl

library(metafor)
library(ggplot2)
library(readxl)

# ── Load extraction table ──────────────────────────────────────────────────
dat <- read_excel("../evidence_tables/PICO2a_integration_benchmark.xlsx",
                  sheet = "extraction")
# Expected columns:
#   study_id | tool | n_cells_log10 | kBET_mean | kBET_sd | n_batches | tissue | runtime_hrs

# ── Subset by tool ─────────────────────────────────────────────────────────
scvi    <- subset(dat, tool == "scVI")
harmony <- subset(dat, tool == "Harmony")

# ── Random-effects meta-analysis (kBET) ────────────────────────────────────
res_scvi    <- rma(yi = kBET_mean, sei = kBET_sd, data = scvi,
                   method = "REML", slab = study_id)
res_harmony <- rma(yi = kBET_mean, sei = kBET_sd, data = harmony,
                   method = "REML", slab = study_id)

# ── Forest plots ───────────────────────────────────────────────────────────
pdf("forest_kBET_scVI.pdf", width = 10, height = 6)
forest(res_scvi, main = "kBET Acceptance Rate — scVI (random-effects)")
dev.off()

pdf("forest_kBET_Harmony.pdf", width = 10, height = 6)
forest(res_harmony, main = "kBET Acceptance Rate — Harmony (random-effects)")
dev.off()

# ── Subgroup analysis by dataset size ──────────────────────────────────────
# Strata: small (<100K), medium (100K–1M), large (>1M)
dat$size_stratum <- cut(10^dat$n_cells_log10,
                        breaks = c(0, 1e5, 1e6, Inf),
                        labels = c("<100K", "100K-1M", ">1M"))

res_subgroup <- rma(yi = kBET_mean, sei = kBET_sd, data = dat,
                    mods = ~ tool + size_stratum,
                    method = "REML")
summary(res_subgroup)

# ── Runtime comparison (median by tool × size stratum) ─────────────────────
cat("\n--- Runtime Summary (median hours) ---\n")
print(aggregate(runtime_hrs ~ tool + size_stratum, data = dat, FUN = median))

# ── Save pooled estimates ──────────────────────────────────────────────────
results <- data.frame(
  tool        = c("scVI", "Harmony"),
  pooled_kBET = c(res_scvi$beta,   res_harmony$beta),
  CI_lower    = c(res_scvi$ci.lb,  res_harmony$ci.lb),
  CI_upper    = c(res_scvi$ci.ub,  res_harmony$ci.ub),
  I2_pct      = c(res_scvi$I2,     res_harmony$I2),
  k_studies   = c(res_scvi$k,      res_harmony$k)
)
write.csv(results, "pooled_kBET_estimates.csv", row.names = FALSE)
cat("Done. Results saved to pooled_kBET_estimates.csv\n")