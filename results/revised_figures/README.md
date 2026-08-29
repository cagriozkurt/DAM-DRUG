# Revised figures — enlarged fonts (reviewer response)

> Reviewer: *"In general, the font size in figures is too small. Please, increase
> the font size employed by 2-3 fold."*

Every figure here is regenerated from the **unchanged** panel logic in
`code/phase6_figures/`, with all font sizes scaled up. Originals in
`results/figures/` are untouched.

## How it was done

`code/phase6_figures/99_revise_all_figures.py` imports each figure script and
runs its `main()` with matplotlib monkeypatched:

| Factor | Default | Effect |
|--------|---------|--------|
| `SCALE`       | **2.0×** | all fonts (titles, axis + tick labels, annotations, legends, colorbars) |
| `TABLE_SCALE` | 1.6×    | matplotlib `Table` text + column widths (fixed-width tables cannot reflow) |
| `SPACE`       | 1.7×    | panel gaps (`wspace`/`hspace`) so larger titles/legends do not collide |
| `FIGSCALE`    | 1.15×   | canvas size, for breathing room |

`graphical_abstract` uses `SCALE = 1.35×` (per-script override) — it is a
hand-placed schematic whose fixed-size boxes overflow at 2×.

Regenerate:

```bash
# local figures (needs synced CSVs; env with matplotlib + seaborn + pandas)
python code/phase6_figures/99_revise_all_figures.py

# tune if a venue wants more/less
SCALE=2.5 SPACE=2.0 python code/phase6_figures/99_revise_all_figures.py

# subset
python code/phase6_figures/99_revise_all_figures.py 02 03
```

## `tiff/` — submission TIFFs

TIFF of every figure, rendered from the vector PDF, 8-bit RGB, lossless
**AdobeDeflate (zip) + horizontal predictor**. `ResolutionUnit` = pixels/inch.

**Main figures (fig1–4, fig_cellchat, graphical_abstract):** the JAD upload
platform rejects TIFFs over **40 megapixels**, so each is rendered at the
highest dpi that stays just under 40 MP (≈ 38.8 MP):

| figure | dpi | pixels | MP |
|--------|-----|--------|----|
| fig1_atlas          | 355 | 7882×4939 | 38.9 |
| fig2_target_prioritization | 412 | 7192×5401 | 38.8 |
| fig3_docking        | 366 | 7464×5200 | 38.8 |
| fig4_validation     | 375 | 7287×5325 | 38.8 |
| fig_cellchat        | 337 | 6494×5989 | 38.9 |
| graphical_abstract  | 445 | 9301×4183 | 38.9 |

(all ≥ 300 dpi — fine for combination art). `supp_fig_S*` stay at **600 dpi**;
they ship inside the supplementary-information PDF, not through the figure
uploader.

Regenerate a main figure at its dpi cap:
```bash
d=$(pdfinfo fig3_docking.pdf | awk '/Page size/{print int(sqrt(39e6/(($3/72)*($5/72))))}')
magick -density $d 'fig3_docking.pdf[0]' -flatten -depth 8 \
  -define tiff:predictor=2 -compress zip tiff/fig3_docking.tiff
```
Supp figures: `magick -density 600 ...` as before.

## Status — all 14 done

fig2, fig3, fig4, fig_cellchat, graphical_abstract, supp_fig_S1–S8.

- **Local (9):** fig2, fig3, fig4, graphical_abstract, S1, S2, S5, S7, and
  fig_cellchat (rebuilt locally after pulling the fresh chord PDF + `lr_summary.csv`).
- **TRUBA (5):** S4 (MD `.xvg`), S6 (raw SEA-AD `.h5ad`), and fig1 / S3 / S8
  (needed `microglia_trajectory.h5ad`).

### Intermediates rebuilt on TRUBA (2026-08-28)

Both were missing from disk and were regenerated under `/arf/scratch/mozkurt/DAM-DRUG`:

| file | job | walltime |
|------|-----|----------|
| `results/phase1/trajectory/microglia_trajectory.h5ad` (1.7 GB) | `05_trajectory_paga.slurm` (6296960) | 24 min |
| `results/phase2/LR/cellchat/cellchat_object.rds` (485 MB) | `17_cellchat_final.slurm` (6296970) | 1 h 01 m |

`02_cellchat_nichechat.R` was patched to read `counts_raw.h5` via raw-binary
sidecars (`cc_*.f32/.i32/.txt`, written by a small python step) — the
`dam-drug-r.sif` image has no `hdf5r`. `03_cellchat_plots.R` chord switched to
`netVisual_chord_cell` (labels outside the ring).

> These intermediates live on TRUBA **scratch** — copy them somewhere durable
> (`/arf/home/...`) if you want to keep them.

## fig_cellchat Panel A

Now a proper circlize chord diagram (`netVisual_chord_cell`) — sector names sit
outside the ring, no overlap with the arcs.

## Source edits

Small layout fixes so the 2× text does not overlap (they also improve the originals):

- `05_fig_cellchat.py` — Panel A given its own wide area, Panel B left-shifted for
  its long y-labels and made taller (`height_ratios`, `figsize` 16×15); dropped
  redundant `[PATHWAY]` suffix on Panel B labels
- `07_supp_fig_s1_regulon_heatmap.py` — footnote folded into the x-axis label
- `09_supp_fig_s3_qc.py` — panel letters moved clear of titles, Braak x-ticks rotated, `n=` labels dropped below the axis
- `12_supp_fig_s6_bhlhe_coexpr.py` — `figsize` 14×20 so the 30 gene y-labels per bar panel don't collide
- `13_supp_fig_s7_slit2_robo2_expr.py` — moved the group bracket + suptitle clear of the panels
- `code/phase2_LR/02_cellchat_nichechat.R` — reads `counts_raw.h5` via raw-binary sidecars (no `hdf5r`)
- `code/phase2_LR/03_cellchat_plots.R` — chord switched to `netVisual_chord_cell` (labels outside the ring)
