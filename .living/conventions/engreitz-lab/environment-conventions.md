# Environment Conventions — Engreitz Lab

**These conventions REPLACE mycelium's `skills/core/references/environment-setup.md`.**
That upstream file recommends `uv` / `conda`. **In this lab we do not use `uv`,
system Homebrew, system Python, or global `pip`/`conda`.** Everything lives in
an isolated per-project environment under `~/Claude/projects/<project>/`.

The decision is by project type, not by taste:

| Project type | Manager | Why |
|---|---|---|
| Python, or **mixed Python + R**, or Snakemake/Nextflow pipelines | **pixi** | One lockfile spans conda-forge + bioconda; handles compiled bio tools and R together |
| **Pure R** (CRAN/Bioconductor only) | **mise** (R) + **renv** | renv lockfile; Bioconductor packages that aren't on conda-forge |
| Unix CLI tools (jq, ripgrep, …) | **mise** | Version-pinned CLIs without touching the system |

## Python / mixed Python+R / pipelines → pixi

```bash
cd ~/Claude/projects/my-project
pixi init                                   # creates pixi.toml (+ pixi.lock)
pixi add python=3.12 numpy pandas           # conda-forge
pixi add --channel bioconda samtools bwa    # bioinformatics tools
pixi add r-base r-tidyverse                  # R in the same env for mixed pipelines
pixi run python script.py
pixi run Rscript script.R
pixi run snakemake --cores 4
```

`pixi.toml` gotchas (these bite silently):
- Use **`[workspace]`**, not the deprecated `[project]`, table.
- **Quote package names containing dots**: `"r-data.table"`, `"bioconductor-org.hs.eg.db"`.
- Version syntax: `"==3.12.0"` (exact) or `"3.12.*"` (minor range); a bare `"3.12.0"` warns.
- On osx-arm64, **`r-base = "4.*"`** (R 3.x has no Apple-Silicon build).
- Bioconductor packages must all come from one release — pin `r-base` to anchor it
  (e.g. `r-base = "4.4.*"` → Bioconductor 3.20).
- `pixi install` writes to `~/Library/Caches/rattler/`; under a sandbox it needs the
  sandbox disabled.

## Pure R (CRAN/Bioconductor) → mise + renv

```bash
cd ~/Claude/projects/my-r-project
mise exec -- Rscript -e "install.packages('renv'); renv::init()"
mise exec -- Rscript -e "renv::install('bioc::DESeq2'); renv::snapshot()"
mise exec -- Rscript script.R
```

If a package fails citing a missing **system library**, STOP and report the
library name — do NOT run `brew install`. The user installs system libs.

## Unix CLI tools → mise

```bash
mise use jq@latest ripgrep@latest
mise exec -- jq '.foo' data.json
```

## Sherlock HPC + Oak/Group storage

Raw data and heavy compute live on Sherlock, **not** in the repo. The repo holds
pointers (see `data-conventions.md`).

- **Storage tiers** — `$GROUP_HOME` is lab-shared and backed up: **non-sensitive
  lab data only**. `$SCRATCH` is per-user, purged, and the place for **sensitive
  or high-risk data**. Never put Stanford High-Risk data on `$GROUP_HOME` or in
  any Git repo / GitHub.
- **Transfers** — always via the data-transfer node:
  `scp`/`rsync` to `engreitz@dtn.sherlock.stanford.edu`, **never the login node**.
- **Sandbox on Sherlock** — do NOT enable Claude Code's `/sandbox` on Sherlock;
  its bubblewrap backend fails on compute nodes. See the lab repo's
  `sandbox-on-sherlock.md`.

## The one rule that overrides everything

**No global installs.** Never `brew install`, never global `pip install`, never
`conda install`/`mamba install` outside pixi, never `install.packages()` outside
an renv project, never `pixi global install`. If a tool is missing, add it to the
project's pixi/mise env. Missing **system** libraries are reported to the user,
not installed.

## ENVIRONMENTS_INSTALLATIONS.md

Mycelium still expects `ENVIRONMENTS_INSTALLATIONS.md` at the repo root — but in
this lab it documents the **pixi/mise/renv** setup, not uv/conda. Record: the
manager, the exact `pixi add` / `renv::install` commands, and any gotcha hit
(e.g. a dot-quoted package, a Bioconductor pin, a reported system lib). Update it
the moment you add a dependency — don't batch.
