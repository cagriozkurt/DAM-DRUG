"""
Graphical Abstract — DAM-DRUG
==============================
Landscape-oriented summary figure for BIB submission.
Redesigned for visual clarity and journal quality.

Output: results/figures/graphical_abstract.pdf
         results/figures/graphical_abstract.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
from pathlib import Path
import numpy as np

PROJECT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Palette ────────────────────────────────────────────────────────────────────
BG       = "#FFFFFF"
PANEL_BG = "#F4F6FA"
C_INPUT  = "#1A5276"   # dark navy  — data
C_TRAJ   = "#2E86C1"   # blue       — trajectory
C_GRN    = "#1E8449"   # green      — GRN
C_SCREEN = "#B7770D"   # amber      — virtual screening
C_CELL   = "#7D3C98"   # purple     — CellChat
C_IKZF   = "#C0392B"   # red        — IKZF1 finding
C_DRUG   = "#117A65"   # teal       — drug finding
C_SLIT   = "#6C3483"   # violet     — CellChat finding
C_ARROW  = "#555555"
C_STEP   = "#ECF0F1"   # step label bg


# ── Helper: rounded box with header stripe ─────────────────────────────────────
def rbox(ax, x, y, w, h, color, title, body_lines,
         title_size=10.5, body_size=8.5, alpha=1.0):
    """Box with coloured header band and white body."""
    # Shadow
    shadow = FancyBboxPatch((x + 0.06, y - 0.06), w, h,
                            boxstyle="round,pad=0.12",
                            linewidth=0, facecolor="#CCCCCC",
                            zorder=2, alpha=0.5)
    ax.add_patch(shadow)
    # White body
    body = FancyBboxPatch((x, y), w, h,
                          boxstyle="round,pad=0.12",
                          linewidth=1.5, edgecolor=color,
                          facecolor="white", zorder=3, alpha=alpha)
    ax.add_patch(body)
    # Coloured header stripe (top 32% of box)
    stripe_h = h * 0.36
    stripe = FancyBboxPatch((x, y + h - stripe_h), w, stripe_h,
                            boxstyle="round,pad=0.12",
                            linewidth=0, facecolor=color,
                            zorder=4, alpha=alpha)
    ax.add_patch(stripe)
    # Clip stripe to box top — overdraw white at bottom of stripe to hide seam
    ax.add_patch(mpatches.Rectangle((x, y + h - stripe_h),
                                    w, stripe_h * 0.18,
                                    color=color, zorder=4))
    # Title text
    ax.text(x + w / 2, y + h - stripe_h / 2,
            title, ha="center", va="center",
            fontsize=title_size, fontweight="bold",
            color="white", zorder=5,
            path_effects=[pe.withStroke(linewidth=1, foreground=color)])
    # Body lines
    n = len(body_lines)
    for i, line in enumerate(body_lines):
        ypos = y + (h - stripe_h) * (1 - (i + 1) / (n + 1))
        ax.text(x + w / 2, ypos, line,
                ha="center", va="center",
                fontsize=body_size, color="#2C3E50", zorder=5,
                linespacing=1.4)


def output_box(ax, x, y, w, h, color, title, body_lines,
               title_size=11.5, body_size=9):
    """Solid-colour output box for findings."""
    shadow = FancyBboxPatch((x + 0.08, y - 0.08), w, h,
                            boxstyle="round,pad=0.15",
                            linewidth=0, facecolor="#AAAAAA",
                            zorder=2, alpha=0.4)
    ax.add_patch(shadow)
    rect = FancyBboxPatch((x, y), w, h,
                          boxstyle="round,pad=0.15",
                          linewidth=0, facecolor=color,
                          zorder=3)
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h * 0.70,
            title, ha="center", va="center",
            fontsize=title_size, fontweight="bold",
            color="white", zorder=5)
    n = len(body_lines)
    for i, line in enumerate(body_lines):
        ypos = y + h * 0.70 - (i + 1) * (h * 0.56 / (n + 0.5))
        ax.text(x + w / 2, ypos, line,
                ha="center", va="center",
                fontsize=body_size, color="white",
                alpha=0.93, zorder=5)


def harrow(ax, x1, x2, y, color=C_ARROW, lw=2.2):
    ax.annotate("", xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=lw, mutation_scale=16),
                zorder=6)


def varrow(ax, x, y1, y2, color=C_ARROW, lw=2):
    ax.annotate("", xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=lw, mutation_scale=14),
                zorder=6)


def main():
    fig = plt.figure(figsize=(18, 8), facecolor=BG)
    ax  = fig.add_axes([0, 0, 1, 1], facecolor=BG)
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 8)
    ax.axis("off")

    # ── Title bar ─────────────────────────────────────────────────────────────
    title_rect = mpatches.Rectangle((0, 7.25), 18, 0.75,
                                     color=C_INPUT, zorder=1)
    ax.add_patch(title_rect)
    ax.text(9, 7.625,
            "Integrative Single-Cell Analysis of AD Microglia: "
            "IKZF1 Identification and Drug Repurposing",
            ha="center", va="center", fontsize=13.5,
            fontweight="bold", color="white", zorder=2)

    # ── Section label: PIPELINE ───────────────────────────────────────────────
    ax.text(0.30, 7.10, "PIPELINE", fontsize=8, color="#888888",
            fontweight="bold", va="center")
    ax.plot([0.28, 17.7], [6.98, 6.98], color="#DDDDDD", lw=1, zorder=1)

    # ── Pipeline boxes (y: 5.15 – 6.90, box height 1.65) ─────────────────────
    PY   = 5.10
    PH   = 1.72
    GAP  = 0.28
    PW   = 2.95
    PX0  = 0.28

    pipeline = [
        (C_INPUT,  "Data Input",
         ["236,002 nuclei · 84 donors", "SEA-AD multi-regional atlas",
          "Braak 0–VI continuum"]),
        (C_TRAJ,   "Trajectory",
         ["PAGA + Diffusion Pseudotime", "HM→IRM→DAM→LAM→LateAD-DAM",
          "5 ordered substates"]),
        (C_GRN,    "GRN Inference",
         ["pySCENIC · 5-seed GRNBoost2", "1.46 M edges · 46 regulons",
          "AUCell + pseudotime scoring"]),
        (C_SCREEN, "Virtual Screening",
         ["1,677 FDA drugs · Vina+Gnina", "6 TF structural targets",
          "MM-GBSA · 100 ns MD"]),
        (C_CELL,   "Cell Communication",
         ["CellChat · 3,718 L-R pairs", "15 cell types · MTG",
          "Pathway-level analysis"]),
    ]

    pipe_centers = []
    for i, (col, title, body) in enumerate(pipeline):
        x = PX0 + i * (PW + GAP)
        rbox(ax, x, PY, PW, PH, col, title, body,
             title_size=10, body_size=8.3)
        cx = x + PW / 2
        pipe_centers.append(cx)
        if i < len(pipeline) - 1:
            harrow(ax, x + PW + 0.04, x + PW + GAP - 0.04,
                   PY + PH / 2, color=col)

    # ── Section label: KEY FINDINGS ───────────────────────────────────────────
    ax.plot([0.28, 17.7], [4.85, 4.85], color="#DDDDDD", lw=1, zorder=1)
    ax.text(0.30, 4.72, "KEY FINDINGS", fontsize=8, color="#888888",
            fontweight="bold", va="center")

    # ── Output boxes (y: 0.85 – 4.45, box height 3.30) ───────────────────────
    OY   = 0.85
    OH   = 3.55
    OW   = 5.10
    OGap = 0.48
    OX   = [0.28, 0.28 + OW + OGap, 0.28 + 2 * (OW + OGap)]

    findings = [
        (C_IKZF, "IKZF1",
         ["Primary late-disease regulator",
          "LateAD-DAM · AUCell = 0.153",
          "Spearman ρ = +0.309 with DPT",
          "Six independent modalities:",
          "DGE · AUCell · DPT · Bulk · ATAC · KO",
          "PROTAC/molecular glue candidate"]),
        (C_DRUG, "Tafamidis  ·  Diflunisal",
         ["Repurposing candidates",
          "Tafamidis \u2192 IRF8  (\u0394G = \u22129.5 kcal/mol)",
          "Diflunisal \u2192 PPARG  (\u0394G = \u22122.8 kcal/mol)",
          "MM-GBSA · 100 ns MD validation",
          "Both: established TTR stabilisers",
          "Convergent pharmacophore hypothesis"]),
        (C_SLIT, "SLIT2 → ROBO2",
         ["Inhibitory interneuron → Microglia",
          "Dominant extrinsic signal in AD MTG",
          "Cumulative probability = 0.318",
          "Pvalb · Lamp5 · Chandelier · Sst",
          "GAS6-MERTK → IRF8 efferocytosis",
          "Novel axis — unreported in human AD"]),
    ]

    out_centers = []
    for i, (col, title, body) in enumerate(findings):
        ox = OX[i]
        output_box(ax, ox, OY, OW, OH, col, title, body,
                   title_size=11.5, body_size=8.8)
        out_centers.append(ox + OW / 2)

    # ── Connecting arrows: pipeline source → finding ──────────────────────────
    # GRN (idx 2) → IKZF1
    # Virtual Screening (idx 3) → Tafamidis/Diflunisal
    # CellChat (idx 4) → SLIT2-ROBO2
    connections = [(2, 0, C_GRN), (3, 1, C_SCREEN), (4, 2, C_CELL)]
    for pi, oi, col in connections:
        sx = pipe_centers[pi]
        dx = out_centers[oi]
        # vertical from pipeline bottom
        varrow(ax, sx, PY, PY - 0.18, color=col, lw=1.8)
        # angled line to output top
        ax.annotate("", xy=(dx, OY + OH + 0.04),
                    xytext=(sx, PY - 0.18),
                    arrowprops=dict(arrowstyle="-|>", color=col,
                                    lw=1.8, mutation_scale=13,
                                    connectionstyle="arc3,rad=0.0"),
                    zorder=6)

    # ── Footer ────────────────────────────────────────────────────────────────
    ax.text(9, 0.30,
            "Data: SEA-AD atlas (Allen Institute for Brain Science) · GSE95587 · "
            "ChEMBL (max_phase = 4) · AlphaFold2 DB · RCSB PDB",
            ha="center", va="center", fontsize=7.5, color="#999999")

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"graphical_abstract.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=BG)
        print(f"Saved {out}")

    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
