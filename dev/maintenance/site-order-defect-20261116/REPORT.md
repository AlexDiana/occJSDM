# Site-order defect in fits with one row per site, 16 November 2026

Written by Claude for Doug. Branch `claude/vigorous-bhaskara-f80e14` (worktree `.claude/worktrees/vigorous-bhaskara-f80e14`), based on `main` at `48b2b3a`. Nothing is pushed.

## Summary

- **The defect.** In a fit with one row per site (the models the package calls `binary` and `continuous`), `runOccJSDM()` sorted the occupancy covariates and the site coordinates by `Site` but left the response in input order. If `data$info$Site` was not already in sorted order, every site was fitted with another site's covariates. There was no error and no warning.
- **Who is affected.** Any one-row-per-site fit whose `Site` column was present and not in sorted order. Names in file order such as `FHT1`, `FHT2`, ..., `FHT10` are not in sorted order (`FHT10` sorts before `FHT2`). A fit is not affected if `data$info` had no `Site` column, or if `Site` was `1:n`, increasing numbers, or zero-padded names in order.
- **Who is not affected.** Occupancy and two-stage fits (repeated rows per site). They sort the response together with the covariates. Verified by test, not assumed.
- **Releases.** All four occJSDM betas (`v0.1.0-beta`, `v0.1.1-beta`, `v0.2.0-beta`, `v0.3.0-beta`) and `main` have the defect, reproduced by running each one. There is no `v0.1.0` tag in occJSDM: `v0.1.0` of 13 November 2026 is the workshop repository's release, which pins occJSDM `v0.3.0-beta`.
- **Shipped and taught results.** None found affected. Every one-row-per-site fit in the lessons, the saved teaching results, the simulation studies and the tests either has no `Site` column or has `Site` as increasing integers. The Quickstart, the shipped example fits and every workshop fit are two-stage. The MCP agent's validation step already refuses one-row-per-site data with an unsorted `Site`. The one known affected analysis is the AlphaEco shakedown, where the defect was found.
- **The fix.** One-row-per-site data now stay in input order throughout. Fits whose `Site` was absent or already sorted have identical draws (22 of 22 compared fits match `main`; see *Verification*).
- **Site labels, added later the same day at Doug's request.** A one-row-per-site fit now stores its `Site` values in `infos$siteNames`, not `1:n`.

## What went wrong

`process_covariates()` in `R/runOccJSDM.R` built the occupancy covariate matrix and the coordinate matrix with `group_by(Site)` followed by `summarise(first)`. dplyr returns one row per group in sorted group order. For occupancy and two-stage data that is correct, because `runOccJSDM()` has already sorted `data$info` and `data$OTU` together with `arrange(Site, Sample, Primer)`, and `group_by()` and `arrange()` use the same order. For data with one row per site there is no `arrange()`: the response `OTU`, `infos$data_info`, `infos$OTU` and the latent factor rows all stay in input order, and only the covariates and coordinates were sorted.

The code has been this way since `process_covariates()` was written (`93f1c62`, 9 July 2026), before the first beta.

Measured on `main` with 150 sites, 6 species and a true slope of 1.5 on one covariate (`reproduce.R` in this folder):

- `Site = 1:n` or `P001`, `P002`, ...: fitted slope 1.17, `fit$X_psi` equals the standardised covariate.
- `Site = FHT1`, `FHT2`, ... (character or factor): fitted slope -0.04, correlation of `fit$X_psi` with the covariate -0.05.
- Shuffled names: slope 0.09. Integers in decreasing order: slope 0.15.
- A factor whose levels are in input order: slope 1.17, because dplyr sorts factors by level.

After the fix all seven cases give a slope of 1.17 and `fit$X_psi` equal to the standardised covariate.

## What was affected, path by path

Each line says whether the path was affected and which test in `tests/testthat/test-site-order.R` covers it. "Affected" means wrong before the fix.

- **Occupancy covariates, binary model: affected.** `X_psi` rows were in sorted site order, the response in input order. Tests: "binary covariates stay on their input rows", "a planted slope is recovered whatever the site labels".
- **Occupancy covariates, continuous model: affected**, same cause. Test: "continuous covariates stay on their input rows".
- **Count responses: not applicable.** `runOccJSDM()` stops with "Counts model not supported yet" on `main`. The count model is a design proposal (PR 42) with no code. When it is written, its one-row-per-site path must use the same input-order rule.
- **Spatial covariates and coordinates (`spatCovariates`), binary and continuous: affected.** `Xs` was built by the same call, so the coordinates were sorted too. Covariates and coordinates agreed with each other but both disagreed with the response, so the spatial field was fitted to the wrong sites. Test: "spatial coordinates stay on their input rows".
- **Occupancy covariates and coordinates, occupancy and two-stage models: not affected.** Test: "occupancy and two-stage matrices follow the sorted site order", with unsorted character, factor and integer names and shuffled rows, and "an occupancy slope is recovered under unsorted site names".
- **Collection covariates (`collCovariates`): not affected.** They exist only in occupancy and two-stage fits, where the sample table is built from the sorted data. Covered by the test above and by the existing `test-collection-alignment.R`.
- **Traits matched to species: not affected.** They are matched by row name with `match()`. Test: "traits are matched to species by name, not by row".
- **Site names stored in `fit$infos`: not wrong, but inconsistent, and now changed.** Occupancy and two-stage fits store the sorted site names in `infos$siteNames`, in the order of every per-site matrix. One-row-per-site fits stored `1:n`, the row number, whatever `Site` held, with the identifiers only in `infos$data_info$Site`. Doug decided on 16 November 2026 to store the `Site` values there instead (`1:n` only when there is no `Site` column). Test: "a one-row-per-site fit names its sites by Site".
- **`predictNewSites()`: affected only through the fit.** It takes the new sites' covariates from the caller and does no site ordering of its own, so its predictions were wrong because the coefficients were wrong. Test: the slope test checks predictions at the fitted sites against the truth (mean absolute error 0.03 after the fix, 0.22 before).
- **`returnOrdinationScores()`, `plotOrdinationScores()`, `plotBiplot()`: affected only through the fit.** They label the rows of the factor scores with `infos$siteNames` and do no ordering of their own. The scores of a misaligned binary fit were estimated from a wrong model. Test: "site labels alone do not change a one-row-per-site fit" checks the scores are identical under sorted and unsorted labels.
- **`returnLatentPresences()` and `plotLatentPresences()`: not affected.** Two-stage only; they read the sorted `infos$data_info` and the stored row indices. Test: "returnLatentPresences() rows carry their own site's probabilities".
- **`computePredictiveOccupancyProbs()`, `computeConditionalOccupancyProbs()`: not affected in occupancy and two-stage fits** (rows named by `infos$siteNames`; tested). In a binary fit they return nothing useful, because the latent-state summaries are `NULL`.
- **`computeSiteWAIC()`, `extractWAIC()`, `compareSiteWAIC()`: affected in binary fits.** The function pairs `infos$OTU` (input order) with `fit$X_psi` (sorted), so it scored each site against another site's covariates, and labelled the result with `infos$data_info$Site` in input order. Occupancy and two-stage fits were not affected. Test: "site WAIC scores each site against its own covariates".
- **Other output functions** (`returnOccupancyCovariates()`, `returnOccupancyGradient()`, `plotCovariateEffect()`, variation partitioning, residual correlations, trait coefficients, convergence diagnostics): no site ordering of their own. In a misaligned binary or continuous fit all of them reported results of the wrong model.

## The fix, and the order chosen

Two changes in `R/runOccJSDM.R`:

- In the "create covariates matrix" block, the occupancy covariates and the coordinates are grouped by `Site` only for occupancy and two-stage data. For one-row-per-site data the rows are taken as they stand.
- In `process_covariates()`, the data frame is made a tibble on both paths, so the stored `infos$X0_psi` keeps the class it had. Without this, an otherwise identical fit differed from `main` in that one attribute.

**The order.** Every per-site result follows the order of the sites in `infos$data_info`, the copy of `data$info` the model was fitted to:

- one row per site: the rows of `data$info` as given, which is order of first appearance, with the sites named in `infos$siteNames` by their `Site` values;
- repeated rows per site: sorted by `Site`, then `Sample`, then `Primer`, with the sites named in `infos$siteNames`.

This is documented in a new section of `?runOccJSDM`, *Row order of the results*.

**Why not order of first appearance everywhere.** The brief preferred it unless the code made sorted order the only consistent choice. For one-row-per-site data it is what the fix does, and it was already what the response, `infos$data_info`, `infos$OTU` and `computeSiteWAIC()` assumed. For occupancy and two-stage data I kept sorted order, for three reasons:

- Those fits are correct now. Reordering their sites would change the site order, and so the seeded numbers, of every existing occupancy and two-stage fit with unsorted names, including the workshop's Gaoligongshan fits, with no gain in correctness.
- Sorted order is built into six places that must agree (`arrange()`, the sample table, `M`, `P`, `K`, the covariates) and is asserted by the existing `test-collection-alignment.R`, which feeds reversed rows and expects sorted results. Changing all of them is not the smallest fix.
- The rows of those models must be regrouped anyway, since replicates of a site need not be adjacent, so "input order" there can only ever mean "sites in order of first appearance", not "rows as given".

If you would rather have one rule, first appearance, for every model, it is a contained follow-up: order by `match(Site, unique(Site))` in those six places and update the alignment tests. It would change existing unsorted occupancy and two-stage fits. **Decided by Doug on 16 November 2026: keep sorted order.**

**One caveat about "sorted".** dplyr sorts character names in the C locale: capitals before lower case, and `S10` before `S2`. Base R's `sort()` and `order()` use the session's locale and can disagree when names mix case. Match by name, through `infos$siteNames`, not by a sort of your own.

## Verification

- **The regression tests fail before the fix and pass after.** Run against an installed copy of `main`, six of the ten tests in `test-site-order.R` fail, 57 of 159 expectations: the binary, continuous, spatial, label-invariance, slope-recovery and site-WAIC tests. The four tests of occupancy, two-stage, latent presences and traits pass before and after. Run against an installed copy of the fix, all 159 expectations pass. The site-label test added afterwards brings the file to 11 tests and 166 expectations; it failed before the label change (6 of its 7 expectations) and passes after.
- **The full suite.** `devtools::test()` on the fixed tree with `NOT_CRAN=true`: 1,971 expectations in 276 tests across 33 files, 0 failures, 0 errors, 0 warnings, 6 skips. The skips are the opt-in coverage study and five identity-across-`ncores` tests that need the reference BLAS (this laptop runs Accelerate); none is new. `test-collection-alignment.R`, `test-site-waic.R`, `test-quickstart-data.R` and `test-regression-bugs.R` pass unchanged. Run on mains power with Low Power Mode off, about 9 minutes. Rerun after the site-label change: 1,978 expectations in 277 tests, 0 failures, 0 errors, 0 warnings, the same 6 skips.
- **Sorted-name fits are unchanged.** `sorted-fits-identity.R` fits the same seeded data under `main` and under the fix: binary, continuous, occupancy and two-stage, each with and without a spatial field; the binary and continuous ones also with an explicit sorted `Site` column as integer, zero-padded character and factor; and short runs of the shipped `sampledata` and `quickstartdata`. With the alignment fix alone (`c7d6cb3`), all 22 fits are `identical()` between the two, comparing the whole returned object except `infos$execution` (process ids and timings). With the site-label change added, 14 of the 22 are still identical, and the 8 one-row-per-site fits with a character or factor `Site` differ only in `infos$siteNames`, which now holds the names: every draw and every other element is identical.
- **Each release reproduces the defect.** Each of the four tags was installed from a `git archive` into its own library and `reproduce.R` was run against it, under R 4.6.1 and dplyr 1.2.1. In all four, sorted labels give a slope of 1.16 to 1.18 and unsorted labels a slope between -0.05 and 0.15, with `fit$X_psi` not matching the covariate. The same script on the fixed code gives 1.17 in every case.
- **Help.** `devtools::document()` regenerated `man/runOccJSDM.Rd`. `tools::checkRd()` passes on it. `NAMESPACE` and `DESCRIPTION` are unchanged.

There are no tests in the repository that refit the shipped `sampleresults` or `quickstartresults` and compare numbers, so "the existing reproduction tests against saved fits" is covered by the identity check above on short runs of the same data, plus the existing structural tests of those objects in the suite.

## Releases

- `v0.1.0-beta` (`ce39b1a`, 7 October 2026): has the defect.
- `v0.1.1-beta` (`cdb29da`, 8 October 2026): has the defect.
- `v0.2.0-beta` (`ef03fa1`, 19 October 2026): has the defect.
- `v0.3.0-beta` (`5cba0db`, 5 November 2026): has the defect.
- `main` (`48b2b3a`): has the defect.
- `v0.1.0`: not a tag of occJSDM, locally or on GitHub (`git ls-remote --tags origin`, `gh release list`). It is the release of the workshop repository `occJSDM-workshop` on 13 November 2026, which pins occJSDM `5cba0db` (`v0.3.0-beta`), so workshop participants have a version with the defect installed.

## Shipped, documented and taught results

Checked by reading the code that builds each `data$info`, and where a saved file holds enough, by reading the file. Verdicts are for one-row-per-site fits; fits with repeated rows per site are listed as not affected by design.

**Package data and the Quickstart**

- `sampledata` and `sampleresults`, `quickstartdata` and `quickstartresults` (`data/`): not affected. Both data sets are two-stage, 1,800 rows over 100 sites, with `Site` numeric in increasing order (read from the `.rda` files). Short seeded refits of both are identical under `main` and under the fix.
- `vignettes/occJSDM-00-quickstart.Rmd`: not affected. Its one fit is the two-stage fit of `quickstartdata`. Line 83 tells readers that data with one row per site are fitted as a JSDM, which is the affected mode, but the page fits none.

**Lessons**

- Lesson 0: no fit. It simulates the two-stage survey.
- Lesson 1: not affected. `fit_one_stage` is an occupancy fit with three samples per site and `Site = 1:n_sites`.
- Lesson 2, the survey fits (`fit`, `alternative_fit`, `unbalanced_fit`): not affected, two-stage.
- Lesson 2, the perfect-observation fit (`perfect_fit`, whose summaries later lessons also read): **not affected**, and this is the one that needed checking, because it is a one-row-per-site fit with an explicit `Site` column. The lesson chunk builds it with `distinct() |> arrange(Site)` and the saved fit was built by `lesson_binary_data()` in `dev/simstudy/vignette-lesson/helpers.R` (lines 89 to 92), which orders the site rows with `order(info$Site)` and takes `OTU` as `z_true[as.character(info$Site), ]`. `Site` comes from the simulator and is numeric: in `vignettes/teaching-data/nonspatial-lesson.rds` it is numeric, 1 to 100, in increasing order. Sorted numeric ids are the safe case.
- Lesson 3: no new fit; it reads Lesson 2's fits.
- Lesson 4: not affected. Its fits are two-stage; the cross-validation helper (`dev/simstudy/waic-cross-validation/cv-lib.R`) subsets rows with a logical mask, which keeps their order.
- Lesson 5 (package comparison): not affected. The occJSDM input is `list(info = input$x, OTU = input$y)` (`dev/simstudy/jsdm-package-comparison/fit-pilot.R`, line 33), and `input$x` holds only the two covariates: in `vignettes/teaching-data/jsdm-comparison.rds`, `training$x` has columns `environment_1` and `environment_2` and no `Site`. With no `Site` column the package uses `1:n`.
- Lesson 6 (traits and calibration): not affected. Same construction: `extension/fit-job.R` line 31 passes `info = input$x`, which `make_input()` in `extension/study.R` (lines 57 to 68) builds from the scaled covariates alone, and the call names every column of it as an occupancy covariate.
- Lesson 7 (spatial sweep), the "true states" arm: not affected. `make_survey()` in `dev/simstudy/spatial-design-sweep/generator.R` (line 116) builds `Site = seq_len(n)`. The eDNA arm is two-stage.
- `vignettes/simulateOccJSDMData.Rmd` and `vignettes/articles/validation.Rmd`: no one-row-per-site fit with a `Site` column. The validation article reports the simulation study below.

**Simulation studies and tests**

- The simulation study's `binary` and `continuous` cells (`tests/testthat/helper-simstudy.R`, `dev/simstudy/run_study.R`), and every other script that fits simulator output directly: not affected. `simulateOccJSDMData()` returns `info` with no `Site` column for `model = "binary"` and `"continuous"` (`R/simulateData.R`, line 379).
- Scripts that build their own one-row-per-site `info`: `validate_factor_correlations.R` (line 58), `validate_spatial_binary.R` (line 122) and `validate_spatial_continuous.R` (line 137) all use `Site = seq_len(n)`; `nonspatial-bias-recheck/run_known_z_control.R` (line 18) orders numeric simulator ids with `order()`. Not affected.
- I did not find a tracked script that fits one row per site with character, factor or shuffled site ids. Studies whose sources live only on archive tags or in ignored `results/` folders were not opened; those built on the simulator are covered by the first point.

**Workshop (`~/src/occJSDM-workshop`, release `v0.1.0`, pins occJSDM `v0.3.0-beta`)**

- Not affected. Every `runOccJSDM()` call in Sessions 2 and 3 and in `analysis/` fits the Gaoligongshan table of 1,818 rows over 101 sites, 18 rows per site. All five saved fits (`data-derived/dhole_fit.rds`, `data-derived/simulated_fit.rds`, their two Posit Cloud copies in `fits/`, and the 108-mammal fit) record `infos$model` as `two_stage`. The exercises (`fit_takin`, `fit_two`, `fit_no_rain`) keep several rows per site. Sessions 1 and 4 do not call `runOccJSDM()`.
- The site names there are not in sorted order in the input (`GLG23033` comes first), which is harmless for a two-stage fit and is exactly the case that would have gone wrong with one row per site.
- What I could not check: anything a participant fitted on their own data during or after the workshop. Participants have `v0.3.0-beta` installed.

**MCP agent (`mcp/paper2agent`, and Doug's copy in `~/src/occJSDM-mcp`)**

- Not affected, and already guarded. The agent's validation step (`mcp/paper2agent/src/r_scripts/quickstart.R`, lines 74 to 77, in the first pilot commit `08381d6` of 5 October 2026) rejects binary data whose `Site` column is not in sorted order, with the message "Pinned binary preparation sorts covariates by unique Site but retains OTU input order". So the behaviour was noticed then and recorded as a property of the pinned package, not reported as a defect. A participant fitting one-row-per-site data through the agent was stopped; one fitting it in plain R was not.
- The repository tracks no MCP example data of its own; its tests and reference runs use the package's two-stage example data.
- After a release with the fix, the guard is no longer needed but does no harm. Remove it when the agent is next re-pinned.

**Known affected**

- The AlphaEco shakedown fits (320 English ponds, site names such as `FHT1`), where the review found the defect. Reported by the session that found it; not re-checked here. Those fits need rerunning on the fixed package.

**Could not tell**

- Fits made by anyone outside this repository: Alex, Qin Pu, workshop participants, or other users of the four betas. Their fits are affected only if they had one row per site and a `Site` column not in sorted order.

## Not done

- **No version bump in these commits.** The version went to 0.3.1 when `v0.3.1-beta` was released on 16 November 2026.
- **No warning for old fits.** The package cannot tell from a saved fit object alone that it was made by an affected version, because fits do not record the package version. A user can check a fit by comparing `fit$X_psi` with their standardised covariates row by row, which the NEWS entry says.
- **`?runOccJSDM` names the returned coordinates `X_s`; the returned element is `Xs`.** A documentation slip noticed in passing, not changed.
- Nothing is pushed and nobody has been told. Doug decided on 16 November 2026 that the fix needs no pull request.

## Decisions for Doug

- **A patch release.** Done: `v0.3.1-beta` was published on 16 November 2026 with the version bumped to 0.3.1 (record in `dev/release/2026-11-16-v0.3.1-beta/`). As written before the release: The fix is small and leaves every correct fit unchanged, and all four betas are wrong for a common kind of input, so I recommend releasing it as `v0.3.1-beta` once Alex has looked at it. That needs the version bumped in `DESCRIPTION` and `CITATION.cff`, and the NEWS heading's "in development" line replaced.
- **Whether to tell workshop participants.** Decided by Doug on 16 November 2026: no need to tell them. As written before the decision: None of the workshop's own fits is affected, since they are all two-stage. The risk is a participant who goes home and fits one row per site (ordinary presence/absence data, the JSDM mode the Quickstart mentions) with site names in file order. A short note with the patch release would cover it: update, and refit any one-row-per-site model whose site names were not sorted. The workshop repository's pin would need moving to the new release for the note to be actionable.
- **Whether to add a line to the README's known limitations now**, (done at the release: the README's known limitations now open with it) before a release, for people on the current betas. The workaround until they update: drop the `Site` column from one-row-per-site data, or sort `info` and `OTU` together by `Site`.
- **Site order for occupancy and two-stage fits.** Decided on 16 November 2026: keep sorted order.
- **`infos$siteNames` for one-row-per-site fits.** Decided on 16 November 2026: store the `Site` values. Done.
- **Alex.** The item is in `TODO.md` under *Needs Alex*. Decided on 16 November 2026: no pull request. Telling him is yours to ask for.

## Files

- `R/runOccJSDM.R`: the fix, the site labels and the new help section.
- `tests/testthat/test-site-order.R`: the regression tests.
- `man/runOccJSDM.Rd`, `NEWS.md`, `TODO.md`: help, changelog and the review-queue item for Alex.
- `reproduce.R` (this folder): the reproduction; run as `Rscript --vanilla reproduce.R <library>` against a library holding the version to test.
- `sorted-fits-identity.R` (this folder): fits the sorted-name cases under one library and saves them; run it once per library and compare the two files with `identical()`.
