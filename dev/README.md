# Development workspace

This directory holds occJSDM development guides, planning documents, scientific validation evidence and teaching-build tools. Maintained source documents are tracked in Git. Large local simulation outputs and other scratch files are ignored. The whole directory is excluded from the R package build and the documentation website; tracked files remain readable in the public repository.

## Directory guide

| Location | Purpose | Tracking |
| --- | --- | --- |
| [parallelisation/](parallelisation/README.md) | Historical parallelisation proposals and reference guides | Tracked |
| [superpowers/specs/](superpowers/specs/) | New general design documents produced through Superpowers | Tracked |
| [superpowers/plans/](superpowers/plans/) | New general Superpowers implementation plans | Tracked |
| [simstudy/](simstudy/) | Study scripts, plans, reports, audits and compact evidence | Tracked, with explicit output exceptions |
| [simstudy/vignette-lesson/](simstudy/vignette-lesson/README.md) | Build and verify the numerical teaching examples | Tracked source and compact artifacts |
| [simstudy/lesson-site/](simstudy/lesson-site/DESIGN.md) | Development and tests for the teaching website | Tracked source; generated preview ignored |
| `simstudy/results/` | Large saved simulation runs and fits | Ignored local output |

Existing study-specific plans and specs remain beside their scripts and evidence in `simstudy/`. The new Superpowers location does not require moving that history. The older `test_sample_l.R` and `test_samplel.R` directly under `dev/` remain tracked historical scripts; other loose scratch files are ignored.

## Superpowers and the workshop

The repository's preferred locations override Superpowers' default `docs/superpowers/` paths. Keep the `superpowers` name so the workflow is visible. [AGENTS.md](../AGENTS.md#development-documents-and-superpowers-locations) records the rule for future agents. The hidden `.superpowers/` working state is separate from these human-readable documents and is not being reorganised.

The [MCP class pilot spec](superpowers/specs/2026-10-04-occjsdm-mcp-class-pilot-design.md) and [implementation plan](superpowers/plans/2026-10-05-occjsdm-paper2agent-mcp.md) describe the Posit Cloud workshop, Paper2Agent conversion, model fallback, Indonesian explanations and adaptive exercises. They are planning documents; the server is not implemented. The optional supported server and student setup will live at root `mcp/`. Disposable MCP feasibility experiments may use `dev/mcp-pilot/`; retaining them in Git requires an explicit tracking decision.

## Generated files and archive review

Root `docs/` is generated pkgdown website output, not a source-document directory. Teaching sources live under `vignettes/`; development utilities that prepare their verified results stay here. Preserve existing generated results until their dependencies are understood.

The current reorganisation leaves `simstudy/results/` untouched. [TODO.md](../TODO.md#review-and-maintenance) records a deferred review to identify which local outputs are needed for rebuilding lessons, reports and validation results, and which can be archived or regenerated. Verify an archive and preserve reproducibility records before removing local copies.
