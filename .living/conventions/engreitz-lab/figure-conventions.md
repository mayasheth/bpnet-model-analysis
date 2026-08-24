# Figure Conventions — Engreitz Lab

**Every figure follows Nature-family styling AND is saved in a vector format
(PDF) alongside the PNG. Never PNG-only.** This is a standing lab rule (see
`CC-Perturb-seq/CLAUDE.md`), enforced here so every mycelium analysis obeys it —
not just ones that reach the report step.

## Always save both PDF and PNG

- **PDF (vector)** is the archival/paper format — editable in Illustrator, fonts
  preserved as text, geometry preserved. **PNG** is for quick viewing / embedding
  in the HTML report.
- A figure saved only as PNG is a bug. If a script writes `fig1.png`, it must
  also write `fig1.pdf`.
- Optional `.eps` when a downstream tool needs it (some projects keep
  `figures/*.{png,pdf,eps}`).

## Use the lab helper — don't hand-roll rcParams

The `/analysis-report` skill ships `nature_style.py`. Copy it into the project
(e.g. `analysis/nature_style.py` or `src/`) and use it:

```python
from nature_style import apply_rcparams, save_fig

apply_rcparams()                         # Nature rcParams — call once in main()
fig, ax = plt.subplots(figsize=(3.5, 2.5))
# ... plot ...
save_fig(fig, "outputs/fig1_coverage.png")   # writes BOTH .png and .pdf
```

`save_fig(fig, path)` writes the `.png` and a sibling `.pdf` by default
(`also_eps=True` for EPS). Do not call `fig.savefig(...png)` directly for a
figure that belongs in a report or the paper.

R users: use the ggplot2 theme in the skill's `r_ggplot_theme.md`, and
`ggsave(..., device = "pdf")` **and** `png` for every figure.

## Nature styling (what `apply_rcparams` sets — and what to match if hand-rolling)

- **Font**: Arial / Helvetica only. 7 pt axis/tick labels, ≤8 pt panel titles;
  no bold except where essential.
- **On-plot titles: report yes, manuscript no.** For **report / quick-view**
  figures keep a short on-figure title — the `/analysis-report` HTML hides the
  legend by default, so the figure must be self-explanatory. For a **final
  manuscript** figure, move the title into the legend (`ax.set_title` off). Same
  `save_fig` call; the distinction is the figure's destination, not two files.
- **Width**: single-column 89 mm (3.5"), double-column 183 mm (7.2"); ≤247 mm tall.
- **Spines**: remove top and right. **Lines**: 0.5–1 pt (0.75 pt default) for axes/ticks.
- **Resolution**: `savefig(dpi=300)`; display dpi 150.
- **Color**: colorblind-friendly palettes; never rely on red–green alone.
- **Panel labels**: bold lowercase a, b, c… at 8 pt, placed at ~(-0.12, 1.05) axes coords.

## Composite / multi-panel figures must stay vector

Do **not** build a composite by reading panel PNGs into matplotlib and saving as
PDF — that bakes a raster into a PDF (no editable text/geometry). Assemble
composites from the panels' own **vector PDFs** (e.g. PyMuPDF
`page.show_pdf_page()`), add panel letters as real text with an embedded font,
and render the preview PNG *from* the vector page so the two never diverge.
Verify vectorness by opening the PDF and checking for selectable text +
vector drawing ops. (Full worked lesson: `CC-Perturb-seq/tasks/lessons.md`.)

## Why this lives in the convention pack

Figures are produced during `/analyze`, before any report step. Putting the rule
here (read into context on every analysis) — not only in `/engreitzlab-report` —
is what makes "always vector + Nature style" hold for every figure, every time.
The `reproducibility-env` reviewer in `/engreitzlab-review` checks it after the fact.
