# Analysis Conventions — Engreitz Lab

Entry point for the `engreitz-lab` convention pack. Read this first, then the
sibling files (`environment-conventions.md`, `data-conventions.md`) as needed.

These conventions layer the lab's own tools onto mycelium's analysis workflow.
Mycelium still owns the repo structure (`analysis/<name>/`, manifests, `.living/`,
the post-action hook protocol); this pack redirects three things to lab tooling.

## 1. Reports → `/analysis-report` (not mycelium's LaTeX report-generator)

When you would run mycelium's `/report`, use the lab's **`/engreitzlab-report`**
overlay instead — it routes to the **`/analysis-report`** skill:

- Quarto (`.qmd`) source rendering to a single self-contained HTML file.
- Nature-style figure legends (headline figure + supporting figures).
- Regeneration command, date, status, goals, findings, open questions.
- A render+lint step that catches broken image links and stale figures.

Do not hand-roll a LaTeX PDF via mycelium's `report-generator` pack — the lab
standard is the Quarto HTML report.

## 2. Genome browser views → `/igv`

For any locus screenshot, session, or track view (ENCODE + lab datasets), use
the **`/igv`** skill rather than describing tracks in prose. Reference the IGV
session/artifact from the analysis doc.

## 3. Every analysis is registered in the EngreitzLabVault (automatic)

At analysis launch, a **lean pointer note** is written to the `EngreitzLabVault`
so the analysis is discoverable lab-wide. This happens automatically via the
lab's vault bridge (`mycelium-vault-bridge/log_analysis_to_vault.py`), invoked
from the `/analyze` post-action step. The note:

- is a `project`-type note (sensitivity `normal` by default);
- records the **question**, the **repo + Sherlock path**, and **datasets/owner
  as `[[wikilinks]]`** — it points back to the repo, it never duplicates data;
- carries `status: active`.

The write is routed through the fail-closed **vault-ingest gate**. If the note
trips the gate (credentials, Stanford High-Risk, handle-with-care/HR language),
it is **skipped and logged** to `.living/`, never written. You do not manage
this by hand — but if `.living/` reports skipped registrations, register those
analyses manually via `/vault-ingest` after resolving the flagged content.

## 4. Figures: Nature-styled, always vector + PNG

Every figure follows Nature styling **and is saved as both a vector PDF and a
PNG** — never PNG-only. Use the `nature_style.py` helper from `/analysis-report`
(`apply_rcparams()` + `save_fig()`, which writes both formats). Full rules,
including composite-figure vectorness, are in `figure-conventions.md` (this pack).
This applies at analysis time (`/analyze`), not just at report time.

## Everything else follows mycelium + robust-analysis

- Every analysis: its own `analysis/<name>/` folder, `UPPER_SNAKE_CASE.md` doc,
  `scripts/`, `outputs/`, `reports/`, and a `run.sh`/`run.py` that reproduces
  final outputs.
- **Figures:** Nature-styled, saved as **both `.pdf` (vector) and `.png`** via
  the `nature_style` helper — see `figure-conventions.md`.
- Follow the **robust-analysis** core pack: fail loudly, assert shapes/types,
  log row counts, sensitivity-sweep every decision, test nulls by
  permutation/bootstrap. Do not subset data without explicit confirmation.
- Run the environment under **pixi/mise+renv** (`environment-conventions.md`);
  keep raw/sensitive bytes on Sherlock (`data-conventions.md`).
