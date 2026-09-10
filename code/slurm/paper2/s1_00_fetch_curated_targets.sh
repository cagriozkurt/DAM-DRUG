#!/bin/bash
# WP1.0 — Build the curated IKZF1 target ground-truth set for the s1_04
# hypergeometric test. Run on a login node (no SLURM). Verified 2026-09-10.
#
#   export DAM_DRUG_DIR=/arf/scratch/mozkurt/DAM-DRUG
#   bash code/slurm/paper2/s1_00_fetch_curated_targets.sh
#
# Sources (experimental / literature TF-target ground truth):
#   ChEA3 Literature ChIP-seq   138 IKZF1 targets
#   ChEA3 ENCODE ChIP-seq      1841
#   ChEA3 ReMap ChIP-seq       1453
#   DoRothEA (via OmniPath, levels A/B/C)  40
# Union (2026-09-10 snapshot): 3,238 unique HGNC symbols.
#
# NOTE: dorothea.opentargets.io/statics/dorothea-data.zip is the 2017
# drug-response dataset, NOT the regulon DB — do not use it. OmniPath serves
# DoRothEA as plain TSV with no R dependency.
set -eo pipefail

PROJDIR=${DAM_DRUG_DIR:-$(pwd)}
OUT=$PROJDIR/data/references
TMP=$OUT/targets_raw
mkdir -p "$TMP"
FINAL=$OUT/curated_ikzf1_targets.txt
PROV=$OUT/curated_ikzf1_targets.provenance.txt

: > "$PROV"
for lib in Literature_ChIP-seq ENCODE_ChIP-seq ReMap_ChIP-seq; do
    gmt="$TMP/$lib.gmt"
    [ -s "$gmt" ] || curl -fsSL -o "$gmt" "https://maayanlab.cloud/chea3/assets/tflibs/$lib.gmt"
    awk -F'\t' '$1 ~ /^IKZF1(_|$)/ {for (i=3; i<=NF; ++i) if ($i != "") print $i}' \
        "$gmt" | tr -d '\r' > "$TMP/ikzf1_$lib.txt"
    echo "ChEA3 $lib : $(wc -l < "$TMP/ikzf1_$lib.txt")" | tee -a "$PROV"
done

# DoRothEA A/B/C via OmniPath (plain TSV, no R)
curl -fsSL "https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B,C&genesymbols=1" \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
                $h["source_genesymbol"]=="IKZF1"{print $h["target_genesymbol"]}' \
  | tr -d '\r' > "$TMP/ikzf1_dorothea.txt"
echo "DoRothEA (OmniPath A/B/C) : $(wc -l < "$TMP/ikzf1_dorothea.txt")" | tee -a "$PROV"

cat "$TMP"/ikzf1_*.txt | grep -E '^[A-Za-z0-9][A-Za-z0-9._-]*$' | sort -u > "$FINAL"
echo "UNION : $(wc -l < "$FINAL")" | tee -a "$PROV"
echo "date  : $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$PROV"
echo
echo "Wrote $FINAL"
echo "s1_benchmark.py reads this path for the IKZF1 hypergeometric test."
