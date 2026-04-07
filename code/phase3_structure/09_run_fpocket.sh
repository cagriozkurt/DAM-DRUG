#!/bin/bash
# DAM-DRUG Phase 3 — Run fpocket on all prepared structures
# ==========================================================
# Loops over data/structures/prepared/*.pdb and runs fpocket on each.
# Output directories are moved to results/phase3/fpocket/.
#
# Usage (inside container):
#   conda run -n base bash code/phase3_structure/09_run_fpocket.sh

set -e

PROJDIR=${DAM_DRUG_DIR:-/Volumes/PortableSSD/untitled\ folder/DAM-DRUG}
PREP_DIR="${PROJDIR}/data/structures/prepared"
OUT_DIR="${PROJDIR}/results/phase3/fpocket"

mkdir -p "${OUT_DIR}"

echo "fpocket version: $(fpocket --version 2>&1 | head -1)"
echo "Structures to process:"
ls "${PREP_DIR}"/*.pdb 2>/dev/null | wc -l

total=0; success=0; failed=0

for pdb in "${PREP_DIR}"/*.pdb; do
    name=$(basename "${pdb}" .pdb)
    dest="${OUT_DIR}/${name}_out"

    if [ -d "${dest}" ]; then
        echo "  [skip] ${name} — output already exists"
        ((success++)) || true
        continue
    fi

    echo "  [run]  ${name}  ($(date +%H:%M:%S))"
    # fpocket creates {stem}_out/ in the same directory as the input PDB
    cd "${PREP_DIR}"
    # Run fpocket; show output if no _out dir is produced
    fpocket_log=$(fpocket -f "${pdb}" 2>&1)
    fpocket_exit=$?
    if [ -d "${PREP_DIR}/${name}_out" ]; then
        mv "${PREP_DIR}/${name}_out" "${dest}"
        ((success++)) || true
        echo "         → done"
    else
        echo "         → FAILED (exit=${fpocket_exit}, no _out dir)"
        echo "${fpocket_log}" | tail -10 | sed 's/^/         LOG: /'
        ((failed++)) || true
    fi
    ((total++)) || true
done

echo ""
echo "=== fpocket complete ==="
echo "  Processed : ${total}"
echo "  Success   : ${success}"
echo "  Failed    : ${failed}"
ls "${OUT_DIR}"
