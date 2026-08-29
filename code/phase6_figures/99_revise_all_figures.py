"""
DAM-DRUG — Regenerate every manuscript figure with enlarged fonts
================================================================
Reviewer request:
    "In general, the font size in figures is too small. Please, increase the
     font size employed by 2-3 fold."

This driver imports each figure script in ``code/phase6_figures/`` and runs its
``main()`` *unchanged*, with a set of matplotlib monkeypatches active:

  1. SCALE       — every font size is multiplied at the single choke point
                   ``FontProperties.set_size`` (titles, axis/tick labels,
                   annotations, legends, colorbars, rcParams defaults).
  2. TABLE_SCALE — matplotlib ``Table`` text + column widths use a gentler
                   factor (fixed-width tables cannot reflow), still a clear bump.
  3. SPACE       — gridspec / subplot wspace+hspace are widened so the larger
                   titles and outside-axes legends do not collide.
  4. FIGSCALE    — the canvas is enlarged a little for extra breathing room.
  5. savefig     — output is redirected from results/figures/ to
                   results/revised_figures/ (originals untouched) and raster dpi
                   is capped (vector PDF unaffected).

Usage
-----
    python code/phase6_figures/99_revise_all_figures.py             # all
    python code/phase6_figures/99_revise_all_figures.py 02 03 05    # subset
    SCALE=2.5 SPACE=2.0 python code/phase6_figures/99_revise_all_figures.py

Figures 01, 09, 10, 12, 14 need TRUBA-only inputs (microglia_trajectory.h5ad,
MD .xvg traces). Run this same script on TRUBA with DAM_DRUG_DIR set to produce
those into results/revised_figures/ as well.
"""
import importlib.util
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.figure as mfigure
import matplotlib.font_manager as mfont
import matplotlib.gridspec as mgridspec
import matplotlib.table as mtable

SCALE       = float(os.environ.get("SCALE", "2.0"))        # general font factor

# Per-script SCALE overrides. The graphical abstract is a hand-placed schematic
# with fixed-size text boxes; 2x overflows them, so it gets a gentler bump.
SCALE_OVERRIDES = {"06_graphical_abstract.py": 1.35}
TABLE_SCALE = float(os.environ.get("TABLE_SCALE", "1.6"))  # matplotlib tables
SPACE       = float(os.environ.get("SPACE", "1.7"))         # panel gap factor
# Per-script SPACE overrides — scripts that already tune their own panel gaps
# (and would just get big dead bands from the global inflation).
SPACE_OVERRIDES = {"01_fig1_atlas.py": 1.15}
FIGSCALE    = float(os.environ.get("FIGSCALE", "1.15"))     # canvas factor
DPI_CAP     = int(os.environ.get("DPI_CAP", "200"))

CODE_DIR = Path(__file__).resolve().parent
PROJECT  = Path(os.environ.get("DAM_DRUG_DIR", str(CODE_DIR.parents[1])))
SRC_DIR  = (PROJECT / "results" / "figures").resolve()
DST_DIR  = (PROJECT / "results" / "revised_figures").resolve()
DST_DIR.mkdir(parents=True, exist_ok=True)

# figure scripts read PROJECT from DAM_DRUG_DIR (or cwd); pin it
os.environ.setdefault("DAM_DRUG_DIR", str(PROJECT))

# small render buffer — font metrics are in points, independent of figure dpi
matplotlib.rcParams["figure.dpi"] = 110

# ── 1. font scaling ──────────────────────────────────────────────────────────
# FontProperties.set_size resolves None / named / float sizes to a plain float
# in self._size; scaling there catches every matplotlib text path exactly once.
_orig_set_size = mfont.FontProperties.set_size


def _scaled_set_size(self, size):
    _orig_set_size(self, size)
    self._size *= SCALE


mfont.FontProperties.set_size = _scaled_set_size

# ── 2. tables — gentler factor + wider columns ───────────────────────────────
# Cell text still flows through FontProperties.set_size (*SCALE), so pre-divide
# to land on *TABLE_SCALE net. Fixed column widths get the same treatment.
_orig_tbl_fontsize = mtable.Table.set_fontsize


def _tbl_set_fontsize(self, size):
    _orig_tbl_fontsize(self, size * TABLE_SCALE / SCALE)


mtable.Table.set_fontsize = _tbl_set_fontsize

_orig_cell_set_width = mtable.Cell.set_width


def _cell_set_width(self, w):
    _orig_cell_set_width(self, w * TABLE_SCALE)


mtable.Cell.set_width = _cell_set_width

# ── 3. panel spacing ────────────────────────────────────────────────────────
_orig_gs_init = mgridspec.GridSpec.__init__


def _gs_init(self, nrows, ncols, figure=None, *, wspace=None, hspace=None, **kw):
    if wspace is not None:
        wspace *= SPACE
    if hspace is not None:
        hspace *= SPACE
    _orig_gs_init(self, nrows, ncols, figure=figure,
                  wspace=wspace, hspace=hspace, **kw)


mgridspec.GridSpec.__init__ = _gs_init
matplotlib.rcParams["figure.subplot.wspace"] *= SPACE
matplotlib.rcParams["figure.subplot.hspace"] *= SPACE

# ── 4. canvas enlargement ───────────────────────────────────────────────────
# Scale figsize once, at creation. Patching set_size_inches compounds because
# savefig(bbox_inches="tight") and layout engines call it again.
_orig_fig_init = mfigure.Figure.__init__


def _scaled_fig_init(self, *args, figsize=None, **kwargs):
    if figsize is None:
        figsize = matplotlib.rcParams["figure.figsize"]
    figsize = (figsize[0] * FIGSCALE, figsize[1] * FIGSCALE)
    _orig_fig_init(self, *args, figsize=figsize, **kwargs)


mfigure.Figure.__init__ = _scaled_fig_init

# ── 5. redirect output ──────────────────────────────────────────────────────
_orig_savefig = mfigure.Figure.savefig


def _redirect_savefig(self, fname, *args, **kwargs):
    try:
        p = Path(os.fspath(fname)).resolve()
        if p.parent == SRC_DIR:
            fname = str(DST_DIR / p.name)
            kwargs["dpi"] = min(kwargs.get("dpi", 300) or 300, DPI_CAP)
    except TypeError:
        pass  # file-like object
    return _orig_savefig(self, fname, *args, **kwargs)


mfigure.Figure.savefig = _redirect_savefig

# ── run the figure scripts ──────────────────────────────────────────────────
SCRIPTS = [
    "01_fig1_atlas.py",
    "02_fig2_targets.py",
    "03_fig3_docking.py",
    "04_fig4_validation.py",
    "05_fig_cellchat.py",
    "06_graphical_abstract.py",
    "07_supp_fig_s1_regulon_heatmap.py",
    "08_supp_fig_s2_celloracle.py",
    "09_supp_fig_s3_qc.py",
    "10_supp_fig_s4_md_rmsd.py",
    "11_supp_fig_s5_af2_quality.py",
    "12_supp_fig_s6_bhlhe_coexpr.py",
    "13_supp_fig_s7_slit2_robo2_expr.py",
    "14_supp_fig_s8_umap_facets.py",
]


def run_script(fname: str) -> bool:
    global SCALE, SPACE
    path = CODE_DIR / fname
    spec = importlib.util.spec_from_file_location(f"_figmod_{path.stem}", path)
    mod = importlib.util.module_from_spec(spec)
    saved_scale, saved_space = SCALE, SPACE
    SCALE = SCALE_OVERRIDES.get(fname, SCALE)
    SPACE = SPACE_OVERRIDES.get(fname, SPACE)
    if (SCALE, SPACE) != (saved_scale, saved_space):
        print(f"  (overrides: SCALE={SCALE} SPACE={SPACE})")
    try:
        spec.loader.exec_module(mod)
        mod.main()
        return True
    except Exception:
        print(f"  FAILED: {fname}")
        traceback.print_exc()
        return False
    finally:
        SCALE, SPACE = saved_scale, saved_space
        plt.close("all")


def main():
    sel = sys.argv[1:]
    todo = [s for s in SCRIPTS
            if not sel or s in sel or s.split("_")[0] in sel]

    print(f"SCALE={SCALE} TABLE_SCALE={TABLE_SCALE} SPACE={SPACE} "
          f"FIGSCALE={FIGSCALE}  ->  {DST_DIR}")
    results = {}
    for fname in todo:
        print(f"\n=== {fname} ===")
        results[fname] = run_script(fname)

    print("\n" + "=" * 60)
    for fname, ok in results.items():
        print(f"  {'ok  ' if ok else 'FAIL'}  {fname}")
    print(f"{sum(results.values())}/{len(results)} figures -> {DST_DIR}")


if __name__ == "__main__":
    main()
