#!/usr/bin/env bash
# Section 2C.1 - Open-access external microglia cohort acquisition
# ==============================================================
# No credentials. Downloads to results/phase7/external/raw/.
# Mathys 2019/2023 and Sun 2023 microglia are Synapse-credentialed and are
# NOT attempted here. Grubman 2019 (GSE138852) is the primary open-access
# fallback: entorhinal cortex snRNA-seq with microglia subclusters + AD/control
# labels. Olah 2020 (GSE146639) full matrix is inside GSE146639_RAW.tar.
set -euo pipefail

DIR="${DAM_DRUG_DIR:-$(pwd)}/results/phase7/external/raw"
mkdir -p "$DIR"; cd "$DIR"
LOG="${DAM_DRUG_DIR:-$(pwd)}/results/phase7/external/ACQUISITION_LOG.md"

fetch () {  # url
  local url="$1" fn; fn="$(basename "$url")"
  if [[ -s "$fn" ]]; then echo "  [skip] $fn exists"; return 0; fi
  if curl -fsSL --max-time 180 -O "$url"; then
    echo "- OK   \`$fn\` ($(du -h "$fn" | cut -f1)) <- $url" >> "$LOG"
  else
    echo "- FAIL $url" >> "$LOG"
  fi
}

echo "# Section 2C acquisition log"           >  "$LOG"
echo ""                                        >> "$LOG"
echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"    >> "$LOG"
echo ""                                        >> "$LOG"

# Grubman 2019 - GSE138852 (primary open-access fallback)
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/GSE138852_counts.csv.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/GSE138852_covariates.csv.gz"

# Olah 2020 - GSE146639 (full matrix is in RAW.tar; large)
# fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146639/suppl/GSE146639_RAW.tar"

echo "" >> "$LOG"
echo "Not attempted (credentialed): Mathys 2019 (syn18485175), Mathys 2023, Sun 2023." >> "$LOG"
cat "$LOG"
