# Post-beta maintenance implementation plan

> **For agentic workers:** Use superpowers:executing-plans for the dependent code tasks. The independent archive review uses superpowers:dispatching-parallel-agents. Check each deliverable against the approved scope before integration.

**Goal:** Retire confirmed dead package helpers, declare genuine data-masked columns, and document the local simulation archive's retention needs.

**Architecture:** Audit the current call graph rather than relying on historical dead-function lists. Preserve public exports, model defaults, native samplers and reference routines used by tests. Keep archive inspection read-only and publish compact metadata and a retention report.

**Tech stack:** R, Rcpp, codetools, testthat, roxygen2, Python metadata inspection and Git.

**Spec:** The user's approval on 7 October 2026 of the three Review and maintenance tasks in TODO.md. The accepted archive deliverable is an inventory and retention plan; moving or deleting raw archives is a separate step.

## Constraints

- Use `.worktrees/post-beta-maintenance` on `codex/post-beta-maintenance`, based on `0c5b1ff`; preserve the primary checkout's existing TODO edits.
- Preserve `computePredictiveProbs`, `partition_r2`, `returnSpatialEffectMean`, `plotSpatialEffect`, `thinOutput` and `.onLoad`. The four-function decision is recorded in `e90a5b2`.
- Preserve all public R exports, live C++ samplers, priors, seed handling and thread behavior. Keep native reference routines that existing tests call.
- Move only helpers whose package callers are all retired or absent. Check non-deprecated tests, development scripts, MCP and vignettes for direct or dynamic calls.
- Remove native export annotations only after checking internal C++ use, then regenerate both Rcpp export files with `Rcpp::compileAttributes()`.
- Declare only real data-masked columns. Preserve reports of genuine undefined symbols in intentionally retained code.
- Raw simulation outputs, frozen runtimes, historical worktrees and fit/provenance records remain untouched.
- Save compact records under `dev/maintenance/`; put bulky checks and metadata under ignored `dev/simstudy/results/`. Use unwrapped prose and no em-dashes.

## Review focus

- Hidden roots: package hooks, name-based test invocation, native reference calls and fixed MCP runtimes must survive the cleanup.
- Shared native callees: a removable R wrapper does not make its C++ body dead.
- Historical scientific evidence: initial and longer fits, selected draws and source/library versions have distinct roles and must not be treated as interchangeable duplicates.
- Data-masking declarations: a column name must not conceal an unrelated executable-code defect with the same spelling.
- Build environment: verify the R version and installed package namespace, avoiding the current framework launcher selecting another R version.

## Task 1: Retire confirmed dead helpers

**Files:** R/mcmcfun.R, R/jsdmfun.R, R/diagnostics.R, src/functions.cpp, src/jsdm.cpp, regenerated R/RcppExports.R and src/RcppExports.cpp; historical bodies under deprecated/R and deprecated/; a compact audit under dev/maintenance/post-beta-maintenance-20261007/.

**Interface:** Input is the current package plus its repository-local caller evidence. Output is the same public API and live model behavior with confirmed unreachable helpers removed.

- [x] Record baseline namespace exports, reachable R functions, direct/dynamic test and research callers, native export callers and package-check diagnostics.
- [x] Classify each candidate as retire, de-export reference body, or preserve; record the exact reason and historical exceptions.
- [x] Move confirmed unreachable R helper groups and wholly dead native bodies to deprecated files, preserving the original bodies and documenting provenance. Preserve all existing tests; name-based sampler tests and research-script callers were found, so their reference helpers stay live.
- [x] De-export serial reference bodies that have no R callers but remain useful native references, leaving those bodies unchanged.
- [x] Regenerate native wrappers with `Rcpp::compileAttributes()`, verify the public export set is unchanged and run the existing meaningful sampler, spatial-algebra, diagnostics and API tests.

## Task 2: Declare genuine data-masked columns

**Files:** R/globals.R if needed, R/occJSDM-package.R and DESCRIPTION if a utils import is needed, generated NAMESPACE, and the compact maintenance report.

**Interface:** Input is the post-cleanup codetools report. Output is a documented list of genuine data-masked columns with remaining actual defects still reported.

- [x] Run codetools/package checking before declarations and classify each remaining diagnostic by exact function and expression.
- [x] Add one `utils::globalVariables()` declaration with only verified column names; record each declared name's real data-masking use.
- [x] Re-run usage checks and inspect residual diagnostics, explicitly preserving errors in the deliberately retained helpers rather than adding their undefined symbols to the declaration.
- [x] Regenerate documentation if import tags change; run the complete source suite, a fresh installed suite and full package check, comparing warnings with the existing baseline.

## Task 3: Review the simulation archive

**Files:** dev/maintenance/archive-review-20261007/ for compact inventory, scripts and report; ignored per-file metadata under dev/simstudy/results/post-beta-maintenance-20261007/.

**Interface:** Input is read-only primary `dev/simstudy/results/` plus existing scripts, manifests, reports and named dependency locations. Output is an inventory and practical retention plan, with evidence and unresolved provenance/location gaps.

- [x] Inventory every file's path, size, type and study grouping without changing raw outputs.
- [x] Map dependencies from lesson rendering, compact-bundle exporters, scientific audits and refitting to inputs, fits, random states, source and library records.
- [x] Verify selected identities against existing manifests and state exactly what was checked; distinguish recorded historical verification from a fresh complete recheck.
- [x] Document what to retain locally, archive with its runtime/provenance, or regenerate. Identify demonstrated duplicates separately from distinct initial/longer fit histories.
- [x] Preserve all raw outputs and require a verified replacement before any later removal.

## Completion

- [x] Independently review the code audit and archive recommendations; resolve important findings.
- [x] Reconcile TODO with actual completion, leaving intentional retention decisions and model extensions open.
- [x] Preserve validation evidence and record exact checks, warning/NOTE comparisons and any limitation. Commit a concrete reviewable branch without overwriting the primary TODO edits.

## Execution rulings

- Preserve `sample_BBsL()` and `loglik_spatialEffect()` as regression/research roots. No test is obsolete solely because a helper is unreachable from public exports.
- Use an isolated R 4.5 launcher and linker configuration for both baseline and final checks; the current framework is R 4.6. Do not change the production spatial kernel to address the mixed-runtime baseline crash.
- The independent review found silent traversal-error handling in the archive inventory script. Capture these errors, fail the inventory on incomplete traversal/stat reads and verify unreadable-directory and empty-archive behavior. A fresh real inventory has identical counts/bytes and no errors.

- Both complete source/private-installed suites pass 1,075 expectations; baseline and corrected final package checks have zero errors and the same three warnings/three NOTEs. The initial extra utils NOTE was fixed with a regenerated explicit import. Both package vignettes build/rebuild successfully.

- Final independent review confirms no unresolved material findings and reconciles the corrected checks, source hashes and completion wording. Keep the committed maintenance branch and worktree for review; integration was not part of this approval. The primary checkout's earlier TODO edits are unchanged and included in this branch.
