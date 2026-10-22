# PR #14 conflict resolution implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task by task. Steps use checkbox (`- [ ]`) syntax for tracking. Recommended execution is one implementer in the existing PR worktree, followed by one independent review of the integrated diff.

**Goal:** Make the approved PR #14 mergeable against current main, preserving its observed-data site-WAIC calculation and main's newer model interfaces, lesson structure and provenance protections.

**Architecture:** Merge main into the existing PR branch with a merge commit, preserving the reviewed commit and avoiding a force-push. Regenerate the native interface from the combined C++ sources, retain main's lesson sources as the starting point, and apply the WAIC teaching changes to their current locations. Keep archived numerical results and the pending scientific validation separate from the integration checks.

**Tech stack:** Git, R, Rcpp/RcppArmadillo, testthat, roxygen2, rmarkdown/knitr, Node's built-in test runner and GitHub CLI.

**Spec:** [Approved PR #14](https://github.com/AlexDiana/occJSDM/pull/14), especially its mathematical explanation and API scope; the existing [site-WAIC proposal](https://github.com/AlexDiana/occJSDM/blob/3eb9cd33c2a7bee1980eb39ba8260e7fe86a3302/dev/simstudy/site-waic-plan.md) and [historical validation report](https://github.com/AlexDiana/occJSDM/blob/3eb9cd33c2a7bee1980eb39ba8260e7fe86a3302/dev/simstudy/site-waic-validation.md); current `AGENTS.md` and `dev/README.md`.

## Inspected state, 7 October 2026

- Remote main and local main: `66ffa783991f8b19f9eb85e01ba499ae1db5508c`.
- Remote PR branch `codex/site-waic`: `3eb9cd33c2a7bee1980eb39ba8260e7fe86a3302`.
- Common ancestor: `3a9726760f692ffe0ef21ca125d0c7acd02c7532`.
- Alex approved the PR head on 6 October 2026. GitHub reports the PR open, approved and conflicting; no status checks are listed.
- A `git merge-tree` simulation found exactly five conflicted paths. Its generated objects were confined to a temporary directory; no branch, working tree or index was merged during planning.
- The existing PR worktree is `/Users/douglasyu/Documents/Codex/2026-09-10/fam/worktrees/occJSDM-site-waic`, on the approved head and clean when inspected.
- The primary checkout has unrelated modifications to `dev/superpowers/plans/2026-10-05-occjsdm-paper2agent-mcp.md` and `dev/superpowers/specs/2026-10-04-occjsdm-mcp-class-pilot-design.md`. Preserve them.

## Global constraints

- This request is for a plan. Implement only when the user asks to execute it. Execution ends with a verified, updated PR ready for merge; it does not include merging into main or publishing the lessons.
- Keep `R/site-waic.R` and `src/site_waic.cpp` identical to the approved PR unless a demonstrated integration failure requires a separately reviewed change. Do not alter the sampler, priors, likelihood target, quadrature policy or supported model types as part of conflict cleanup.
- Preserve the distinction between default observed-data site WAIC, explicit `extractWAIC(..., type = "legacy")`, and the legacy scalar in `results_output$WAIC`.
- Preserve main's lesson numbering, teaching explanations, plot changes, publication exclusions and saved three-sample teaching bundles. Alex's code approval does not replace Doug's lesson review.
- Do not modify recorded archive hashes or expand the plot-only source-compatibility allowlist to make historical archive verification pass. Do not refit models or overwrite teaching bundles in this work.
- Write planning-document paragraphs on single lines. Use no em-dashes. Render lesson Markdown with `--wrap=none`.
- Do not stash, reset or incorporate unrelated work from the primary checkout. Recheck both remote tips and worktree cleanliness before execution; revisit this inventory if either tip has changed.

## Conflict decisions

| Conflicted file | Resolution |
| --- | --- |
| `src/RcppExports.cpp` | Regenerate from combined C++ sources. Retain main's 17-argument registrations for both `sample_BBsL_cpp` and `sample_BBsL_parallel`, and PR #14's 6-argument `site_loglik_cpp`. Regenerate `R/RcppExports.R` at the same time. |
| `TODO.md` | Start from main's current structure. Update its existing WAIC entries for Alex's approval and pending integration; do not restore the PR's duplicate maintenance item or discard newer archive/release tasks. |
| `vignettes/LESSON-PLAN.md` | Start from main. Keep Lessons 0-7 and all completed teaching work; update only the WAIC item and its current Lesson 4 location. |
| `vignettes/occJSDM-lesson-3.Rmd` | Start from the complete main version, not the textual auto-merge. Retain its outputs/diagnostics lesson and make only the small WAIC cross-reference correction. Move the PR's substantive prediction/WAIC additions to Lesson 4. |
| `vignettes/occJSDM-lesson-3.md` | Regenerate from the resolved Rmd. Do not resolve this large generated-text conflict by hand. Also regenerate Lesson 4's Markdown after its source changes. |

## Review focus

1. Native registration can silently retain the obsolete 16-argument sampler interface. Task 1 checks the compiled registration counts and Task 3 runs the existing intercept-prior tests.
2. Cleanly merged R files can still lose main's prior metadata, covariate fixes or threshold behavior. Tasks 1 and 3 inspect the integrated diff and test threshold conversion and the current sampler interface.
3. Keeping both sides of the lesson conflict would reinsert obsolete teaching material. Task 2 uses main's complete sources, then checks the small intended edits and rendered lessons.
4. September's two-sample WAIC evidence must not be presented as a validation of October's three-sample fits. Task 2 labels historical evidence and retains unevaluated examples and unchanged saved results.
5. Modified exporter code invalidates recorded source hashes even when its numerical output is unchanged. Task 2 preserves the provenance gate and documents the separate historical reproduction route; Task 3 runs the source-check regression tests.

## Task 1: Integrate main and regenerate the native interface

**Files:** Resolve `src/RcppExports.cpp`; regenerate it and `R/RcppExports.R`; inspect the auto-merges in `R/runOccJSDM.R`, `R/output.R`, `NAMESPACE`, `_pkgdown.yml` and `man/runOccJSDM.Rd`. Preserve the new scorer, its help files and its tests from the PR.

**Interfaces:** `sample_BBsL_cpp(..., sigma_b0)` and `sample_BBsL_parallel(..., sigma_b0)` retain main's 17 arguments. `site_loglik_cpp(eta, log_g0, log_g1, loadings, nodes, log_weights)` has 6. The public scorer signatures and site/legacy behavior remain those approved in PR #14.

- [ ] Check the existing PR worktree's ownership, status and head, and reuse it if still clean and available. Follow the worktree skill and native attachment workflow where applicable. If it cannot be reused, create an isolated checkout from the approved PR head; never perform this merge in the dirty primary checkout.
- [ ] Fetch `origin`, verify the inspected heads against this plan, and run `git merge --no-commit --no-ff origin/main` from the PR branch. Expect the five listed conflicts, and inspect any additional conflict before proceeding. All subsequent relative paths and commands refer to this isolated worktree.
- [ ] Use main's `src/RcppExports.cpp` as the temporary generated-file baseline, retain all combined annotated C++ sources, then run `Rscript -e 'Rcpp::compileAttributes(".")'`. Inspect both generated files and stage the resolved generated conflict. Do not hand-maintain the registration table.
- [ ] Inspect the auto-merged R code against main. Confirm that `infos$threshold` is added without replacing raw `infos$OTU`; retain main's simultaneous threshold masks, `noise_prior`, `spatial_sd_prior` and `intercept_prior`. Retain current covariate-response and prediction code. Confirm the only new public exports are `computeSiteWAIC` and `compareSiteWAIC`, alongside the revised `extractWAIC` interface.
- [ ] Build the combined native library and verify its registrations with the command below. Expected result: successful compilation and all three assertions pass.

```sh
Rscript -e 'pkgload::load_all(".", recompile = TRUE, quiet = TRUE); r <- getDLLRegisteredRoutines("occJSDM")[[".Call"]]; stopifnot(r[["_occJSDM_sample_BBsL_cpp"]]$numParameters == 17L, r[["_occJSDM_sample_BBsL_parallel"]]$numParameters == 17L, r[["_occJSDM_site_loglik_cpp"]]$numParameters == 6L)'
```

Do not create an intermediate merge commit with unresolved documentation. Task 4 creates the single integration merge commit after all checks.

## Task 2: Reconcile teaching, status and archive provenance

**Files:** Resolve `TODO.md`, `vignettes/LESSON-PLAN.md`, `vignettes/occJSDM-lesson-3.Rmd` and its Markdown; edit `vignettes/occJSDM-lesson-4.Rmd` and regenerate its Markdown; reconcile `dev/simstudy/vignette-lesson/prediction-export.R`, `prediction-verify.R` and `README.md`; annotate `dev/simstudy/site-waic-plan.md` and `site-waic-validation.md`. Preserve main's publication configuration and all files in `vignettes/teaching-data/`.

**Interfaces:** Existing `prediction_examples$manifests[[arm]]$waic` values continue to mean the archived legacy scalar. Corrected-score teaching uses `computeSiteWAIC(fit, threshold = 1)` and `compareSiteWAIC(first, second)` in `eval=FALSE` examples. No new numerical teaching bundle is produced.

- [ ] Replace both conflicted Lesson 3 files with main's versions before editing. In the Rmd's `Can WAIC compare these fits?` subsection, preserve the warning that perfect-observation and PCR fits have different responses. Explain that the legacy values are unsuitable for new-site selection, while the corrected site criterion is described in Lesson 4 and still needs reliability checks. Retain the current link anchor and surrounding teaching.
- [ ] In Lesson 4's `Check the additional fit and understand the WAIC limitation` section, change the displayed archived-value calls to `extractWAIC(..., type = "legacy")` and label the table `Legacy stored score`. Preserve its existing numerical values. Add the approved `computeSiteWAIC`/`compareSiteWAIC` examples with the actual old-fit threshold `1`, in unevaluated chunks. Explain that the returned difference is first minus second, lower WAIC is preferable for identical observations, and the paired SE describes variation across sites rather than MCMC or integration error. Keep factor-count selection unfinished.
- [ ] Label the September two-sample check explicitly as historical if mentioned. Do not copy its 98/100 warning count or 2.20-point difference into claims about the current three-sample lesson fits. Current-fit rescoring and held-out observed-survey validation remain separate work; this integration does not claim those have been completed.
- [ ] Resolve `LESSON-PLAN.md` from main, retaining the current sequence and completed trait/simulator work. Update the WAIC task to say the implementation was approved on 6 October and selection validation remains open. Resolve `TODO.md` from main, updating its existing approval/status and future-validation entries without reinstating the PR's duplicate maintenance paragraph. Do not mark the PR merged before it actually is.
- [ ] Preserve the cleanly merged explicit `fit$results_output$WAIC` reads in the prediction exporter and verifier. Retain main's `lesson_source_check()` calls, `export_source_hashes`, and current Lesson 4 lookup. Change the verifier's `Current augmented WAIC` message to `Legacy stored score`, and update the build README's corresponding description.
- [ ] Document the provenance consequence in the build README: the new namespace/C++/R code is outside the existing plot-only allowlist, and `prediction-export.R` itself is included in the compact bundle's `source_md5`. Existing archive verifiers must be run with their matching historical source and library; they are not a merged-source compatibility check. Keep their gates intact. Defer re-exporting/refitting under merged source to the existing release-rebuild task, without changing recorded hashes or teaching bundles here.
- [ ] Add dated approval/integration notes to the existing site-WAIC plan and report. Preserve September's hashes, test totals and numerical results as dated evidence; identify the former Lesson 3 as current Lesson 4 where needed. Append the new integration checks after running them, with their actual outcomes. Keep `validate_site_waic.R`'s historical archive layout unchanged because this plan does not undertake a new scoring study.
- [ ] Preserve `_config.yml`, `.Rbuildignore`, `vignettes/.install_extras` and the lesson `published` flags from main. Retain the scorer entries added to `_pkgdown.yml`. No website deployment is part of this work.

## Task 3: Verify integrated behavior and rendered documentation

**Files:** Extend `tests/testthat/test-site-waic.R` only with the focused integration check below if no equivalent check exists by execution time. Use the existing package tests, lesson link tests and source-provenance tests. Regenerate roxygen output only from the resolved sources, reviewing any unrelated generated changes before retaining them.

**New integration test:** `new thresholded fits score the same observations as pre-binarized fits`. For both `occupancy` and `two_stage`, generate a non-spatial fixture, ensure its raw response matrix contains boundary values `0, 1, 2, 3`, and make a second copy with `OTU = 1 * (OTU >= 2)`. Fit the raw data at threshold `2` and the binary copy at threshold `1` with the same fixed seed, fixture MCMC settings, covariates, one factor and `listPriors = list(sigma_b0 = 2)`. Call `runOccJSDM()` directly because `fit_fixture()` has no threshold argument. Assert the two fits have identical posterior results, raw reads are retained in the first fit, thresholds are saved as `2` and `1`, the intercept-prior metadata retains SD `2`, matched-draw site log likelihoods agree, and explicit legacy extraction equals each fit's stored scalar. Compare `log_lik` directly; `compareSiteWAIC()` intentionally rejects different threshold provenance. The test protects main's threshold fix and changed sampler interface together with the PR's post-processing contract; it does not assert a stochastic fitted estimate equals truth.

- [ ] Add that targeted test using the existing fixtures, without altering the shared fixture interface. Run the focused scorer tests first and investigate failures before broadening verification. Existing tests already cover exact state enumeration, shared-factor integration, rotations/scale, missing data, matched chains, old-fit threshold requirements and site identity.
- [ ] Run the full source suite once, with opt-in long recovery studies left off. It includes current intercept-prior, factor-reparameterisation, RNG, covariate-response and API regression tests. Expected result: no failures or unexpected warnings; record actual counts rather than repeating September's total.

```sh
Rscript -e 'testthat::test_local(".", filter = "site-waic", reporter = "summary", stop_on_failure = TRUE)'
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)'
```

- [ ] Generate/reconcile public help with the repository's roxygen2 version (`8.0.0` in `DESCRIPTION`). Confirm that main's new fitting arguments and all site-WAIC help survive; inspect generated changes rather than accepting a broad documentation rewrite.
- [ ] Render Lessons 3 and 4 to Markdown and HTML in separate fresh R sessions, loading the integrated source first. Use the pattern below for each lesson. Existing chunks must read the tracked compact teaching bundles and must not fit models or execute the new WAIC examples. Inspect both rendered lessons for misplaced sections, broken links, raw math, changed numeric tables and unintended figure changes. Generated Markdown is canonical output of the resolved Rmd, not an independently edited source.

```sh
Rscript -e 'pkgload::load_all(".", quiet = TRUE); suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble); library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format = rmarkdown::github_document(html_preview = FALSE, pandoc_args = "--wrap=none"))'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble); library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format = "rmarkdown::html_vignette")'
```

- [ ] Run the three independent documentation/provenance checks below. All must pass. These test the source-check mechanism, not compatibility of the historical full-fit archives with changed package code.

```sh
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
Rscript dev/simstudy/vignette-lesson/test_source_check.R
node --test dev/simstudy/lesson-site/test_lessons.js
```

- [ ] Run an installed-package check, with manual and vignette builds disabled because the changed lessons have been rendered explicitly. Use `rcmdcheck::rcmdcheck(".", args = c("--no-manual", "--no-build-vignettes"), error_on = "error")` if available, otherwise equivalent `R CMD build` and `R CMD check` commands. Require zero errors and investigate new warnings/notes. If a failure appears pre-existing, reproduce it on the recorded main in a separate clean checkout before classifying it that way.
- [ ] Inspect the final diff against main. `R/site-waic.R` and `src/site_waic.cpp` must match the approved head, the publication configuration and teaching bundles must match main, and no conflict markers or unintended generated/build files may remain. Run both `git diff --check` and `git diff --cached --check`, and confirm `git diff --name-only --diff-filter=U` is empty.

## Task 4: Review, update the existing PR and report readiness

**Files:** Final resolved changes from Tasks 1-3, including the dated integration record. No changes to unrelated MCP planning documents or primary-checkout work.

- [ ] Have one independent reviewer inspect the final diff against main and the resolution against the approved PR. Ask specifically about native argument counts, main-side preservation, Lesson 3/4 placement, threshold provenance and historical versus current teaching evidence. Resolve material findings and rerun only the affected checks.
- [ ] Confirm the remote heads have not changed since the merge began. Commit the completed merge on the existing PR branch, with a message such as `Merge main into site WAIC proposal`. Do not squash or rebase away the approved commit as part of conflict resolution.
- [ ] Push normally to `origin/codex/site-waic`, then inspect GitHub's mergeability, checks and approval state for the new head. Preserve Alex's recorded approval of the old head without describing it as a fresh review of the integration commit; honour any approval requirement GitHub now reports.
- [ ] Update the PR description's stale review-status wording and add the actual integration validation results. Preserve its mathematical explanation and use `\mathrm` rather than `\operatorname`, which failed in GitHub's renderer in this PR. Verify the saved description and any changed rendered content.
- [ ] Report the updated head, whether GitHub considers it mergeable, test/check outcomes, and the still-open scientific validation. Leave merging into main and lesson publication outside this plan.

## Completion criteria

All five conflicts are resolved with main's newer work retained; the approved likelihood implementation is unchanged; regenerated native registrations have argument counts 17, 17 and 6; source tests and package checks pass subject only to explicitly verified existing diagnostics; Lessons 3 and 4 render with unchanged archived results and correct API labels; historical evidence and provenance checks remain honest; the existing PR is updated and its current merge/review status is reported. A ready-to-merge PR is not a completed factor-count validation or permission to publish the teaching lessons.
