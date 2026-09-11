#!/bin/bash
# WP3.0 — Set up the ligand-prep toolchain on TRUBA (no SLURM; login node).
# Idempotent. Verified end-to-end 2026-09-11.
#
#   ssh truba
#   bash /arf/scratch/mozkurt/DAM-DRUG/code/slurm/paper2/s4_00_setup_truba_toolchain.sh
#
# Neither scenic.sif nor lipogate-env.sif has the full docking-prep chain
# (rdkit+openbabel+vina+meeko+acpype+pdb2pqr). What IS usable:
#   - conda env `pld3` (/arf/home/mozkurt/miniconda3/envs/pld3, someone else's
#     project env, python 3.11): rdkit 2025.09.5, openbabel 3.1.0 (+ obabel
#     binary), vina 1.2.7, numpy, pandas. Used READ-ONLY here — nothing is
#     installed INTO this env.
#   - meeko / acpype / pdb2pqr: not present anywhere -> pip install --user
#     against pld3's python (writes only to ~/.local, not the shared env).
#   - acpype 2023.10.27 bundles its own AmberTools (amber_linux/{bin,lib}/) —
#     antechamber, parmchk2, teLeap, sqm all present. Its PyPI wheel is
#     missing SONAME symlinks for 3 bundled libs, and its bundled libcrypto.so.3
#     is incompatible with the system libssl.so.3 that gets pulled in
#     transitively. Both are host/package-layout bugs, fixed by symlinks
#     below (no compilation, no root, nothing outside ~/.local).
set -eo pipefail

PLD3=/arf/home/mozkurt/miniconda3/envs/pld3/bin/python
ACP="$HOME/.local/lib/python3.11/site-packages/acpype/amber_linux/lib"

echo "[1/3] pip install --user meeko acpype pdb2pqr gemmi (against pld3's python)"
"$PLD3" -m pip install --user --quiet meeko acpype pdb2pqr gemmi

echo "[2/3] verify imports"
"$PLD3" - << 'PY'
import importlib
for m in ("rdkit", "openbabel", "vina", "meeko", "acpype", "pdb2pqr"):
    mod = importlib.import_module(m)
    print(f"  {m:10s} OK {getattr(mod, '__version__', '')}")
PY

echo "[3/3] fix acpype's bundled AmberTools shared-library layout"
cd "$ACP"
# missing SONAME symlinks (the wheel ships the fully-versioned .so files but
# not the SONAME the dynamic linker looks up)
[ -e libhdf5_hl.so.310 ] || ln -s libhdf5_hl.so.310.0.2 libhdf5_hl.so.310
[ -e libhdf5.so.310 ]    || ln -s libhdf5.so.310.2.0    libhdf5.so.310
[ -e libzip.so.5 ]       || ln -s libzip.so.5.5          libzip.so.5
# libbz2: point at the system copy (SONAME differs: system ships .so.1.0.8,
# teLeap wants .so.1.0)
[ -e libbz2.so.1.0 ]     || ln -s /usr/lib64/libbz2.so.1.0.8 libbz2.so.1.0
# libcrypto.so.3: acpype's own copy is incompatible with the system libssl.so.3
# that gets pulled in transitively (OPENSSL_3.0.1 symbol not found). acpype
# code hard-sets LD_LIBRARY_PATH to exactly this directory before calling
# teLeap, so the only way to make it fall through to the system's compatible
# libcrypto is to remove (rename, not delete) the bundled one from here.
[ -e libcrypto.so.3.bak_moved_by_dam_drug ] || \
    { [ -e libcrypto.so.3 ] && mv libcrypto.so.3 libcrypto.so.3.bak_moved_by_dam_drug; }

ls -la libhdf5* libzip.so.5 libbz2.so.1.0 libcrypto.so.3* 2>&1
echo
echo "Setup complete. Smoke-tested 2026-09-11: acpype -c gas -a gaff2 -o gmx"
echo "produced complete GROMACS topology (test_GMX.{gro,itp,top}) end to end."
echo "acpype -c bcc (AM1-BCC/sqm) also runs but needs a clean-valence input"
echo "mol2 (an odd-electron test case failed on sqm) -- s4_01 uses -c gas."
