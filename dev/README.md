# Development workspace

This directory holds occJSDM development guides, planning documents, scientific validation evidence and teaching-build tools. Maintained source documents are tracked in Git. Large local simulation outputs and other scratch files are ignored. The whole directory is excluded from the R package build and the documentation website; tracked files remain readable in the public repository.

## Directory guide

| Location | Purpose | Tracking |
| --- | --- | --- |
| [parallelisation/](parallelisation/README.md) | Historical parallelisation proposals and reference guides | Tracked |
| [superpowers/specs/](superpowers/specs/) | New general design documents produced through Superpowers | Tracked |
| [superpowers/plans/](superpowers/plans/) | New general Superpowers implementation plans | Tracked |
| [release/](release/2026-10-07-beta/REPORT.md) | Dated release refits, provenance, check reports and announcement drafts | Explicitly tracked |
| [maintenance/](maintenance/post-beta-maintenance-20261007/REPORT.md) | Package cleanup audits and local simulation-archive retention review | Explicitly tracked compact records; bulky evidence ignored |
| [simstudy/](simstudy/) | Study scripts, plans, reports, audits and compact evidence | Tracked, with explicit output exceptions |
| [simstudy/vignette-lesson/](simstudy/vignette-lesson/README.md) | Build and verify the numerical teaching examples | Tracked source and compact artifacts |
| [simstudy/lesson-site/](simstudy/lesson-site/DESIGN.md) | Development and tests for the teaching website | Tracked source; generated preview ignored |
| `simstudy/results/` | Large saved simulation runs and fits | Ignored local output |

Existing study-specific plans and specs remain beside their scripts and evidence in `simstudy/`. The new Superpowers location does not require moving that history. The older `test_sample_l.R` and `test_samplel.R` directly under `dev/` remain tracked historical scripts; other loose scratch files are ignored.

## Superpowers and the workshop

The repository's preferred locations override Superpowers' default `docs/superpowers/` paths. Keep the `superpowers` name so the workflow is visible. [AGENTS.md](../AGENTS.md#development-documents-and-superpowers-locations) records the rule for future agents. The hidden `.superpowers/` working state is separate from these human-readable documents and is not being reorganised.

The [MCP class pilot spec](superpowers/specs/2026-10-04-occjsdm-mcp-class-pilot-design.md) and [implementation plan](superpowers/plans/2026-10-05-occjsdm-paper2agent-mcp.md) describe the Posit Cloud workshop, Paper2Agent conversion, model fallback, Indonesian explanations and adaptive exercises. The optional four-tool server and pilot setup now live at root [mcp/](../mcp/README.md). Local verification and user-run Cloud execution passed on the recorded pinned runtime. The [rehearsal review](../mcp/paper2agent/reports/comparison-review.md) supports MCP as an execution interface, while controlled comparison and workshop adoption remain pending. Users learn from the vignettes and lessons; scientific and Indonesian review must precede classroom use. Disposable MCP feasibility experiments may use `dev/mcp-pilot/`; retaining them in Git requires an explicit tracking decision.

## Generated files and archive review

Root `docs/` is generated pkgdown website output, not a source-document directory. Teaching sources live under `vignettes/`; development utilities that prepare their verified results stay here. Preserve existing generated results until their dependencies are understood.

The [7 October archive review](maintenance/archive-review-20261007/REPORT.md) inventories 87.53 GiB in `simstudy/results/`, maps the dependencies needed for lessons, reports and numerical audits, and recommends retaining complete scientific studies. All raw outputs remain untouched. Targeted integrity checks are limited to 1.87% of raw bytes; no replacement archive or restore has been verified. [TODO.md](../TODO.md#review-and-maintenance) keeps the unresolved Lesson 7 raw-archive location open. Verify a complete archive and preserve reproducibility records before removing local copies.
