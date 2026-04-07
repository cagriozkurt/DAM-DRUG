"""
DAM-DRUG Supplementary Figure S5 — MD Simulation RMSD Traces
=============================================================
Panels (2×2 grid):
  A. PPARG (1FM9-LBD) + Diflunisal (CHEMBL898): ligand RMSD 100 ns
  B. IRF8 (AF2-DBD) + Tafamidis (CHEMBL2103837): ligand RMSD 100 ns
  C. IRF8 (AF2-DBD) + Chlorpromazine (CHEMBL42): ligand RMSD 100 ns
  D. IRF8 (AF2-DBD) + Paroxetine (CHEMBL490): ligand RMSD 100 ns

All traces: core-RMSD (ligand after fit to binding-site backbone).
Running average (window=50 ps) overlaid on raw trace.

Runs locally:
  python code/phase6_figures/supp_fig_s5_md_rmsd.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

PROJECT = Path("/Volumes/PortableSSD/untitled folder/DAM-DRUG")
MD_DIR  = PROJECT / "results/phase4/md"
OUT     = PROJECT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})

# ── XVG parser ────────────────────────────────────────────────────────────
def parse_xvg(path):
    """Return (time_ns array, rmsd_nm array) from a GROMACS .xvg file."""
    times, vals = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                vals.append(float(parts[1]))
    return np.array(times), np.array(vals)

def running_avg(arr, window=50):
    """Simple uniform running average."""
    kernel = np.ones(window) / window
    return np.convolve(arr, kernel, mode="same")

# ── Panel definitions ─────────────────────────────────────────────────────
PANELS = [
    {
        "label": "A",
        "target_dir": "PPARG_1FM9_LBD_prep",
        "chembl": "CHEMBL898",
        "drug": "Diflunisal",
        "target": "PPARG (1FM9-LBD)",
        "color": "#56B4E9",
        "xvg": "core_ligand_rmsd.xvg",   # prefer core RMSD; fall back to ligand_rmsd.xvg
    },
    {
        "label": "B",
        "target_dir": "IRF8_AF2_DBD_prep",
        "chembl": "CHEMBL2103837",
        "drug": "Tafamidis",
        "target": "IRF8 (AF2-DBD)",
        "color": "#E69F00",
        "xvg": "core_ligand_rmsd.xvg",
    },
    {
        "label": "C",
        "target_dir": "IRF8_AF2_DBD_prep",
        "chembl": "CHEMBL42",
        "drug": "Chlorpromazine",
        "target": "IRF8 (AF2-DBD)",
        "color": "#009E73",
        "xvg": "core_ligand_rmsd.xvg",
    },
    {
        "label": "D",
        "target_dir": "IRF8_AF2_DBD_prep",
        "chembl": "CHEMBL490",
        "drug": "Paroxetine",
        "target": "IRF8 (AF2-DBD)",
        "color": "#CC79A7",
        "xvg": "core_ligand_rmsd.xvg",
    },
]

# ── Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
axes_flat = axes.flatten()

for ax, panel in zip(axes_flat, PANELS):
    run_dir = MD_DIR / panel["target_dir"] / panel["chembl"]
    xvg_path = run_dir / panel["xvg"]
    if not xvg_path.exists():
        xvg_path = run_dir / "ligand_rmsd.xvg"

    t, rmsd = parse_xvg(xvg_path)

    # Convert time to ns (GROMACS outputs in ps)
    t_ns = t / 1000.0
    # Convert RMSD to Å (GROMACS outputs in nm)
    rmsd_A = rmsd * 10.0

    # Running average
    window = min(50, len(rmsd_A) // 10)
    avg_A  = running_avg(rmsd_A, window=window) if window > 1 else rmsd_A

    # Plot raw (thin, low alpha) + average (thick)
    ax.plot(t_ns, rmsd_A, color=panel["color"], alpha=0.25, linewidth=0.6, rasterized=True)
    ax.plot(t_ns, avg_A,  color=panel["color"], alpha=0.95, linewidth=1.8,
            label=f"{window}-frame avg")

    ax.set_xlabel("Time (ns)", fontsize=8)
    ax.set_ylabel("Ligand RMSD (Å)", fontsize=8)
    ax.set_title(f"{panel['drug']}\n{panel['target']}", fontsize=8.5, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, t_ns[-1])

    # Stability annotation
    last20 = rmsd_A[int(len(rmsd_A) * 0.8):]
    mean20 = last20.mean()
    std20  = last20.std()
    ax.axhline(mean20, color=panel["color"], linestyle="--", linewidth=0.8, alpha=0.7)
    ax.text(0.97, 0.95, f"Last 20 ns: {mean20:.2f}±{std20:.2f} Å",
            transform=ax.transAxes, ha="right", va="top", fontsize=7,
            color=panel["color"])

    # Panel label
    ax.text(-0.12, 1.04, panel["label"], transform=ax.transAxes,
            fontsize=12, fontweight="bold")

fig.suptitle(
    "Supplementary Figure S5 — MD Simulation Ligand RMSD Traces (100 ns)",
    fontsize=10
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"supp_fig_S5_md_rmsd.{ext}",
                bbox_inches="tight", dpi=300)
    print(f"Saved → {OUT}/supp_fig_S5_md_rmsd.{ext}")

plt.close(fig)
