# Expanded Quickstart Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn `vignettes/occJSDM-00-quickstart.Rmd` into one long, theory-free tour of every main occJSDM result, built on a new simulated Gaoligongshan-style river survey shipped as `quickstartdata` and `quickstartresults`.

**Architecture:** A dev script simulates the survey with planted effects and saves `data/quickstartdata.rda`; a second script fits it (no spatial field, current default priors), slims the fit and saves `data/quickstartresults.rda`. The vignette loads the precomputed fit and calls only existing package functions, plus short ggplot/dplyr code for the figures and tables the package does not draw. `sampledata`, `sampleresults` and the lessons are untouched.

**Tech Stack:** R 4.6.1 (Accelerate BLAS), occJSDM, testthat 3, roxygen2, rmarkdown/knitr, ggplot2, dplyr.

**Spec:** `dev/superpowers/specs/2026-11-01-quickstart-expansion-design.md` (read it first; this plan argues from it).

## Global Constraints

- No new package functionality; do not call `returnLatentPresences()` (broken) or `plotVariationPartitioning()` (ternary).
- No spatial field in the fit; `b0_prior` left at the current default.
- Survey: 100 sites, 3 water samples per site, 2 primers, 3 PCR replicates per primer (6 PCRs per sample, 1,800 rows of `info`).
- Occupancy covariates: elevation, Human Modification Index (continuous), protected versus not protected (categorical). Collection covariate: water volume filtered. Traits: body mass and Red Listed (0/1). `n_factors = 2`.
- About 20 species named `Class_Order_n` (for example `Mammalia_Carnivora_1`); no real species names.
- Convergence: a handful of flagged rows under `rhat > 1.01 | ess < 400` (this plan pins "a handful" as 2 to 10 rows); the flag the page discusses should be `theta0` of a planted contaminated species if at all possible.
- `data/quickstartresults.rda` well under 25 MB.
- Vignette prose: `AGENTS.md`, *Writing the vignettes and lessons*, plus Doug's uncommitted Lesson 1 rules in `/Users/douglasyu/src/occJSDM/AGENTS.md` (main checkout, read only). No theory, no comparison with truth, no limitations beyond one pointer to the README's *Known limitations*. Every number from the data via inline R. "Hidden environmental covariates", with "latent factors" once to explain `n_factors`; no "factors"/"loadings" otherwise.
- Markdown: one line per paragraph, `-\>` and `\~` escaped in prose, no em-dashes (also in commit messages).
- Computing: at most 5 fitting workers on the laptop; launch fits longer than 2 minutes with `nohup` and check `pmset -g` power mode before timing; never `R CMD build` the main checkout (build from `git archive`); run commands from the worktree `.worktrees/quickstart-expansion` only.
- Do not edit anything in the main checkout `/Users/douglasyu/src/occJSDM`; copy Doug's two Quickstart edits from its working copy by hand (Task 4).

## Review Focus

- **Factor reference level:** `protected` must have levels `c("not_protected", "protected")`, so the reported coefficient is "protected compared with not protected"; a wrong order silently flips the story. Pinned in Task 1's test.
- **New-site prediction with a factor:** `predictNewSites()` must be given `protected` with the same levels as the fit; a character column or a new level must not silently produce a wrong design matrix. Pinned in Task 5 by building the new-site data from the fitted levels and checking the protected and not-protected predictions differ in the planted direction.
- **Slimming must not break a function the page calls:** every package function the vignette uses must run on the slimmed fit. Pinned in Task 2's test, which calls each one.
- **Prose that outlives the numbers:** the fit will be redone with `b0_prior = "community"` later, so no sentence may assert a direction or ranking the data could reverse without inline R behind it. Checked in Task 6's render review.
- **0/1 and continuous traits on different scales:** say in the text whether each trait coefficient is per standard deviation or per unit, as `runOccJSDM()` actually treats them (check its standardisation of `traits`). Task 5.

---

### Task 1: Simulate and ship `quickstartdata`

**Files:**
- Create: `dev/quickstart-data/simulate-quickstartdata.R`
- Create: `data/quickstartdata.rda`, `dev/quickstart-data/truth.rds` (planted values, not shipped; `dev/` is build-ignored)
- Modify: `R/data.R` (add the `quickstartdata` entry)
- Test: `tests/testthat/test-quickstart-data.R`

**Interfaces:**
- Produces: `quickstartdata`, a list with
  - `info`: data.frame, 1,800 rows, columns `Site` (int 1:100), `Sample` (int, unique across survey), `Primer` (int 1:2), `elevation` (numeric, metres, roughly 200 to 5,000), `hmi` (numeric, 0 to 1), `protected` (factor, levels `c("not_protected", "protected")`), `water_volume` (numeric, litres, sample level);
  - `OTU`: numeric read-count matrix, 1,800 rows, about 20 columns named `Class_Order_n`;
  - `traits`: numeric matrix, rows = `colnames(OTU)`, columns `body_mass` (kg, log-normal spread) and `red_listed` (0/1).
- Produces: `dev/quickstart-data/truth.rds`, a list of the planted values (coefficients, `theta0`, `p`, `q`, trait-by-environment matrix) for the builder's own checks only.

- [ ] **Step 1: Write the failing test** `test_that("quickstartdata has the agreed survey design", ...)` asserting:

```r
d <- quickstartdata
expect_named(d, c("info", "OTU", "traits"))
expect_equal(nrow(d$info), 1800L)
expect_equal(length(unique(d$info$Site)), 100L)
expect_true(all(table(unique(d$info[c("Site", "Sample")])$Site) == 3L))
expect_equal(sort(unique(d$info$Primer)), 1:2)
expect_true(all(table(interaction(d$info$Sample, d$info$Primer, drop = TRUE)) == 3L))
expect_identical(levels(d$info$protected), c("not_protected", "protected"))
expect_true(all(c("elevation", "hmi", "water_volume") %in% names(d$info)))
expect_equal(nrow(d$OTU), nrow(d$info))
expect_true(ncol(d$OTU) >= 18L && ncol(d$OTU) <= 22L)
expect_true(all(grepl("^[A-Z][a-z]+_[A-Z][a-z]+_[0-9]+$", colnames(d$OTU))))
expect_identical(rownames(d$traits), colnames(d$OTU))
expect_identical(colnames(d$traits), c("body_mass", "red_listed"))
expect_true(all(d$traits[, "red_listed"] %in% c(0, 1)))
expect_true(any(d$traits[, "red_listed"] == 1) && any(d$traits[, "red_listed"] == 0))
```

- [ ] **Step 2: Run it to see it fail.** Run: `Rscript -e 'devtools::test(filter = "quickstart-data")'`. Expected: FAIL, object `quickstartdata` not found.

- [ ] **Step 3: Write `dev/quickstart-data/simulate-quickstartdata.R`.** Usage `Rscript dev/quickstart-data/simulate-quickstartdata.R`, run from the worktree root, one fixed seed. `simulateOccJSDMData()` cannot plant named covariates, a two-level factor or chosen effects, so either post-process its output or simulate directly with the same model structure (occupancy -\> field collection with `theta0` -\> PCR with `p`, `q`; 2 hidden environmental covariates). Plant, strongly enough to show at 100 sites: Red Listed species respond more positively to `protected` (trait-by-protected); larger-bodied species decline more steeply with `hmi` (trait-by-hmi); an elevation response that varies among species; one or two species with clearly non-zero `theta0`; primers that differ by species, so neither primer is best for all. Save with `save(quickstartdata, file = "data/quickstartdata.rda", compress = "xz")` and the planted values to `truth.rds`.

- [ ] **Step 4: Document `quickstartdata` in `R/data.R`** in the style of the `sampledata` entry (columns, units, the factor's levels, the trait coding, that it is simulated and modelled on Ji et al. 2025), then `Rscript -e 'devtools::document()'`.

- [ ] **Step 5: Run the test.** Run: `Rscript -e 'devtools::test(filter = "quickstart-data")'`. Expected: PASS.

- [ ] **Step 6: Commit.** `git add dev/quickstart-data/simulate-quickstartdata.R dev/quickstart-data/truth.rds data/quickstartdata.rda R/data.R man/quickstartdata.Rd tests/testthat/test-quickstart-data.R` and commit "Add the simulated river survey for the Quickstart".

### Task 2: Fit, slim and ship `quickstartresults`

**Files:**
- Create: `dev/quickstart-data/fit-quickstartresults.R`, `dev/quickstart-data/RECORD.md` (chain lengths tried, timings with power mode, flags, final size)
- Create: `data/quickstartresults.rda`
- Modify: `R/data.R` (add the `quickstartresults` entry)
- Test: `tests/testthat/test-quickstart-data.R` (extend)

**Interfaces:**
- Consumes: `quickstartdata` (Task 1).
- Produces: `quickstartresults`, the output of

```r
runOccJSDM(
  data = quickstartdata,
  listParams = list(n_factors = 2),
  threshold = 1,
  occCovariates = c("elevation", "hmi", "protected"),
  collCovariates = c("water_volume"),
  MCMCparams = list(nchain = 2, nburn = NBURN, niter = NITER, nthin = NTHIN, ncores = 2),
  summarisedLatentPresences = TRUE
)
```

with seed `20261101`, then slimmed. The vignette's live-fit chunk (Task 4) is exactly this call with the final `NBURN`, `NITER`, `NTHIN`.

- [ ] **Step 1: Write the failing test** `test_that("quickstartresults is the agreed non-spatial fit and every Quickstart function runs on it", ...)`:

```r
f <- quickstartresults
expect_identical(f$infos$model, "two_stage")
expect_equal(f$infos$n, 100L)
expect_equal(f$infos$S, ncol(quickstartdata$OTU))
expect_equal(f$infos$n_factors, 2)
rda <- test_path("../../data/quickstartresults.rda")  # absent under R CMD check
if (file.exists(rda)) expect_lt(file.size(rda), 25e6)
diag <- returnConvergenceDiagnostics(f)
n_flag <- sum(diag$rhat > 1.01 | diag$ess < 400)
expect_gte(n_flag, 2L); expect_lte(n_flag, 10L)
expect_silent({
  computePredictiveOccupancyProbs(f); computeConditionalOccupancyProbs(f)
  returnOccupancyRates(f); returnOccupancyCovariates(f); returnTraitsCoeff(f)
  returnCollectionCovariates(f); returnResidualCorrelationMatrix(f)
  returnOrdinationScores(f); returnVariationPartitioning(f)
})
for (p in list(plotOccupancyRates(f), plotDetectionRates(f), plotCollectionRates(f),
               plotStage1FPRates(f), plotStage2FPRates(f), plotFPTPStage2Rates(f),
               plotOccupancyCovariates(f, covName = "hmi"),
               plotTraitsCoefficients(f, covName = "hmi"),
               plotResidualCorrelationMatrix(f), plotBiplot(f)))
  expect_s3_class(p, "ggplot")
```

Adjust each call's required arguments to the function's real signature (read `?fn`); the point is that every package function Tasks 4 to 6 use is called here once. If a call is not silent for a benign reason (a message), use `expect_no_error` for it and say why in a comment.

- [ ] **Step 2: Run it to see it fail.** Expected: FAIL, `quickstartresults` not found.

- [ ] **Step 3: Pilot.** Write `fit-quickstartresults.R` (usage `Rscript dev/quickstart-data/fit-quickstartresults.R NBURN NITER NTHIN OUT.rds`; follow the shape of `dev/release/2026-10-13-parallel-chains-refit/refit-sampleresults.R`: record git revision, source hashes, warnings seen, timings, `sessionInfo()`, BLAS; no `OCCJSDM_REFERENCE_BLAS` gate, the laptop uses Accelerate). Run a short pilot (for example 1,000 + 1,000) detached with `nohup`, after `pmset -g`, and record its time and flag count in `RECORD.md`.

- [ ] **Step 4: Tune and run the real fit.** Choose `NBURN`/`NITER` so the table flags 2 to 10 rows, ideally including `theta0` of a planted contaminated species; if `theta0` is not among them, record which parameters are and pick the one with the best story for the page. Choose `NTHIN` so the saved draws keep the file under 25 MB. Record each attempt in `RECORD.md`.

- [ ] **Step 5: Slim.** Measure component sizes (`object.size` over `results_output` and `results_output$jsdm_output`). Thinning is the main lever. Drop an array only if no function in Step 1's test reads it (check by grepping `R/output.R` for the element name), and confirm the test still passes. Save with `compress = "xz"` to `data/quickstartresults.rda`; write `returnConvergenceDiagnostics()` to `dev/quickstart-data/diagnostics.csv`.

- [ ] **Step 6: Document `quickstartresults` in `R/data.R`** (the exact call, seed, R version, BLAS, what was slimmed), then `devtools::document()`.

- [ ] **Step 7: Run the test.** Run: `Rscript -e 'devtools::test(filter = "quickstart-data")'`. Expected: PASS.

- [ ] **Step 8: Commit** the script, `RECORD.md`, `diagnostics.csv`, `data/quickstartresults.rda`, `R/data.R`, `man/quickstartresults.Rd`, the test: "Fit and ship the Quickstart's example results".

### Task 3: Run the full test suite and a package check

**Files:** none new.

- [ ] **Step 1:** `RCPP_PARALLEL_NUM_THREADS=1 Rscript -e 'devtools::test()'`. Expected: all pass, except tests already known to be skipped. Run interrupt tests in the foreground.
- [ ] **Step 2:** `git archive HEAD | tar -x -C "$TMP/qs"` and `R CMD build --no-build-vignettes` then `R CMD check --no-manual` on the export. Expected: no new ERROR or WARNING compared with the baseline in `AGENTS.md`, *Expected R CMD check output*. Record any new NOTE (for example installed size) in `RECORD.md`.

### Task 4: Vignette, part 1: data, fit, convergence

**Files:**
- Modify: `vignettes/occJSDM-00-quickstart.Rmd` (rewrite sections *The example data* to *Check that the chains agree*)

**Interfaces:**
- Consumes: `quickstartdata`, `quickstartresults` (Tasks 1 and 2). Inside the vignette, name them `info`, `OTU`, `traits` and `fit` as the text does.

- [ ] **Step 1: Opening.** Keep the current first paragraph; replace the second with what this page does (format the data, fit, check convergence, produce every main result as a figure and a table) and say it skips theory; one sentence pointing to the lessons and README *Known limitations*.
- [ ] **Step 2: The example data.** Introduce the survey (simulated, modelled on Ji et al. 2025, cited and linked once). Show `head()` of `info`, `OTU`, `traits`. Explain, as bullets: one row per PCR; the site-level columns repeat within a site; `protected` is a factor and how to make one from your own column (`factor(x, levels = ...)`, reference level first); `water_volume` varies by sample; `red_listed` is 0/1; read counts go in as they are and a PCR with at least `threshold` reads counts as a detection. Keep the `View()` exercise idea.
- [ ] **Step 3: Fit the model.** The live-fit chunk (`eval=FALSE`) is Task 2's exact call; carry over Doug's edits (`ncores = 2` on its own continuation line; "hidden environmental covariates (statistically, latent factors)" for `n_factors`). Bullets explain each argument. Drop the `spatCovariates` line and its note; one sentence says the spatial field is left out. Then `fit <- quickstartresults`. Keep the timing/BLAS paragraph, shortened.
- [ ] **Step 4: Check convergence.** `returnConvergenceDiagnostics()`, `summary()` of `rhat` and `ess`, the flagged rows via the rule of thumb (counts by inline R), a `plotTraceplot()` of the featured flagged parameter beside a healthy one of the same kind, and a short "what to do" bullet list: does it matter for your question; run longer and recheck; if it will not settle, the data may not pin it down; report it. One sentence on `runOccJSDM()`'s end-of-fit warnings, naming the blocks it warned about (from `RECORD.md`; inline R if they can be recomputed, otherwise "for this fit" with the names). State what a longer run did, from `RECORD.md`.
- [ ] **Step 5: Render and check.** `Rscript -e 'rmarkdown::render("vignettes/occJSDM-00-quickstart.Rmd")'` (render from the worktree with the worktree's package installed or `devtools::load_all()`). Expected: renders with no errors or R warnings; every number in the new prose comes from inline R (`grep -nE '[0-9]+\.[0-9]+|[0-9]{2,}' ` on the prose lines, each hit justified).
- [ ] **Step 6: Commit** "Rewrite the Quickstart's data, fit and convergence sections".

### Task 5: Vignette, part 2: results

**Files:**
- Modify: `vignettes/occJSDM-00-quickstart.Rmd` (replace *Look at the results*)

**Interfaces:**
- Consumes: `fit` (Task 4).
- Produces, as named R objects later chunks reuse: `occ_site` (long data frame: `Site`, `species`, `predictive`, `conditional`), `richness` (data frame: `Site`, `group`, `expected_richness`, where `group` is `"All"`, each class and each order), `detection_profile`, `cov_effects`, `trait_env`, `assoc_pairs` (built here or in Task 6, but named exactly so).

Each subsection opens with why it is there and closes with one plain sentence of what it tells the reader, numbers by inline R. Order:

- [ ] **Step 1: Occupancy.** `plotOccupancyRates()` for baseline occupancy with intervals; `computePredictiveOccupancyProbs()` and `computeConditionalOccupancyProbs()` into `occ_site`, with one sentence each on the question it answers (habitat prediction for maps and unsurveyed places; presence given the detections at this site) and that both are posterior means. Show a small map-like plot (site by species tile or one species against elevation) only if it says something the table does not.
- [ ] **Step 2: Environmental effects.** `plotOccupancyCovariates()` per covariate, including the protected contrast (read as "protected compared with not protected"); `plotOccupancyGradient()` or `plotCovariateEffect()` response curves for elevation and `hmi`.
- [ ] **Step 3: Traits.** `plotTraitsCoefficients(fit, covName = ...)` for `hmi` and `protected`; then the added figure: each species' estimated slope (posterior mean from `returnOccupancyCovariates()`) against its trait value, with the line implied by the trait coefficient. State the trait scale (per SD or per unit) as `runOccJSDM()` treats it.
- [ ] **Step 4: Collection, detection and false positives.** One Ji et al. Supp. Fig. 3-style figure: per species, collection probability (at mean water volume) and PCR detection probability for each primer with 95% intervals, false-positive rates as short marks on the left margin, built from the draws the `plot...()` functions read; then `plotCollectionCovariates()` for the water-volume effect. Point to `plotStage1FPRates()`/`plotStage2FPRates()` for the separate views.
- [ ] **Step 5: Residual associations, ordination, variation partitioning.** `plotResidualCorrelationMatrix()`; `plotBiplot()` (and `plotOrdinationScores()` if it adds something); `returnVariationPartitioning()` drawn as one stacked bar per species with ggplot (measured covariates versus hidden environmental covariates; drop all-zero spatial columns).
- [ ] **Step 6: Expected richness.** From `computeConditionalOccupancyProbs()`, `rowSums` over all species, each class and each order (groups parsed from the species names) into `richness`; one figure. Say these are means and that intervals need per-draw presences, which the package does not yet return.
- [ ] **Step 7: Prediction at new sites.** Build a small new-site data frame (for example an elevation gradient, `hmi` low and high, `protected` both levels using `factor(..., levels = levels(info$protected))`) and call `predictNewSites(fit, X_psi = ..., useBiotic = FALSE)` (check the argument names and the expected design of `X_psi` in `?predictNewSites`). Plot one or two species. The protected prediction for a Red Listed species must exceed the not-protected one, as planted; if not, stop and investigate before writing prose.
- [ ] **Step 8: Would more sampling help?** `plotCumulativeSpeciesDetections()` for field samples and for PCR replicates (each primer), one sentence on the answer.
- [ ] **Step 9: Render and check** as in Task 4, Step 5. Also check each "what this tells you" sentence against the rendered figure.
- [ ] **Step 10: Commit** "Add the Quickstart's results tour".

### Task 6: Vignette, part 3: tables, where next, and the record

**Files:**
- Modify: `vignettes/occJSDM-00-quickstart.Rmd` (tables, *Where to go next*)
- Modify: `TODO.md` (one item), `NEWS.md` (one entry)

- [ ] **Step 1: Tables.** One subsection, each table a tidy data frame shown with `head()` and a one-line `write.csv()` example (`eval=FALSE`): `occ_site` plus raw evidence (positive samples out of M, positive PCRs out of all PCRs, computed from `info`/`OTU` at `threshold = 1`); `richness`; `detection_profile` (species: baseline occupancy, collection probability, detection per primer, field and lab false-positive rates, each mean and 95% interval); `cov_effects` (species by coefficient, each non-reference level separately, mean, interval, `excludes_zero`); `trait_env` (trait by covariate, same columns); `assoc_pairs` (pairs whose interval excludes zero); the convergence table. Name columns in the reader's terms.
- [ ] **Step 2: Where to go next.** `?runOccJSDM`; Lesson 3 for reading these outputs against a known truth; the simulation vignette; README *Known limitations* (the page's single pointer to limitations); the drafts note.
- [ ] **Step 3: Full render review.** Render; read the HTML top to bottom against the spec's *Goals*, *Non-goals* and *Anti-goals*, the Review Focus list and the writing rules; fix anything found.
- [ ] **Step 4: TODO and NEWS.** In `TODO.md`, next to the existing item about refitting the lesson fits with the community mean (or under *Studies, package work and documentation* if absent), add one bullet: refit `quickstartresults` with `b0_prior = "community"` once PRs 36 to 38 merge, rerun `dev/quickstart-data/fit-quickstartresults.R`, re-tune for a handful of flags, re-render. In `NEWS.md`, one entry for the expanded Quickstart and the two new datasets.
- [ ] **Step 5: Commit** "Add the Quickstart's tables and record the follow-up refit".

### Task 7: Hand-off render for Doug

- [ ] **Step 1:** Rerun Task 3's test suite and check.
- [ ] **Step 2:** Render the final HTML to the worktree root as `quickstart-review.html` (not committed) and send it to Doug with a short list of: the featured convergence flag, the final chain settings and file size, and anything in `RECORD.md` he should know.
