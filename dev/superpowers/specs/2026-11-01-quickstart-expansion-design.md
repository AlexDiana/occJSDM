# Expanded Quickstart: design

Agreed with Doug in conversation on 1 November 2026. Branch `claude/quickstart-expansion`, worktree `.worktrees/quickstart-expansion`.

## Intent

Lesson 2 covers a lot of theoretical ground and is not much use to someone who wants to explore what occJSDM can do. The Quickstart (`vignettes/occJSDM-00-quickstart.Rmd`) becomes that tour. "Quick" means the reader runs the functions and reads off the results, without going into the theory or the limitations of occJSDM. It is one long page, not split into several vignettes.

The reader has their own two-stage eDNA data and wants a working path from raw data to every main result. By the end they should be able to:

- format the input data for a two-stage model;
- fit the two-stage model with occupancy, collection and trait covariates, including a categorical occupancy covariate, but no spatial field;
- diagnose convergence with R-hat, ESS and traceplots, and know what to do about a flagged parameter;
- generate all the main results occJSDM can produce, as figures and as tables.

Doug expects readers mostly to paste the tables and figures into papers, but he wants occJSDM to be useful to managers too.

## Goals

- Every result section runs package functions and states, in one plain sentence, what the result tells the reader.
- Every main result is available as a table (a tidy data frame) as well as a figure, with a line showing how to write it to CSV.
- The example data feel like a real survey: a simulated river survey modelled on Ji et al. (2025, *Ecology Letters*; in Zotero), the Gaoligongshan aquatic eDNA study.

## Non-goals

- No theory, no comparison with the simulated truth (that is Lesson 3's job), and no discussion of occJSDM's limitations beyond one pointer to the README's *Known limitations*.
- No new package functionality. The Quickstart shows only what the package already returns. Per-site occupancy intervals, richness intervals, richness by higher taxon as a package helper and the broken `returnLatentPresences()` are recorded in `TODO.md` under *Site-level results with intervals* (added 1 November 2026, on main, not yet committed) and are not built here.
- Lesson 3 and the other lessons are not changed. `sampledata` and `sampleresults` are not changed.

## Anti-goals

- A page that turns into a second Lesson 2 or 3: theory, caveats or truth comparisons creeping in.
- Simulated results that could be mistaken for Ji et al.'s real findings (hence generic species names, below).
- Hard-coded numbers in the text (see the writing rules below).

## Decided by Doug

- **Example data, a new dataset used only by the Quickstart.** Not a replacement for `sampledata`, which Lessons 0, 2 and 3, the simulation vignette, the README and the API tests use (Lesson 0 has passed Doug's line-by-line review). Working names `quickstartdata` and `quickstartresults` (Claude's call; Doug did not object).
- **Survey design:** modelled on Ji et al. 2025: about 100 river sites, 3 water samples per site, **2 primers** with **3 PCR replicates per primer** (6 PCRs per sample, as in Ji et al., who used one primer; Doug chose two so the per-primer comparison stays).
- **Occupancy covariates:** elevation, Human Modification Index (both continuous) and **protected versus not protected** (categorical; Ji et al.'s `streamOrigin`, whether a stream's headwaters rise inside a reserve).
- **Collection covariate:** water volume filtered.
- **Traits:** body mass and IUCN Red Listed or not (0/1).
- **No spatial field** in the fit.
- **Prior:** the current default `b0_prior`. Refit with `b0_prior = "community"` once PRs 36 to 38 merge, together with the lessons (add this to `TODO.md`).
- **Convergence:** the example fit deliberately keeps a handful of flagged parameters under the Quickstart's rule of thumb (`rhat > 1.01` or `ess < 400`), so the reader sees what a problem looks like and what to do about it. Doug accepted a handful rather than exactly one or two, to avoid hours of tuning runs. The flag the page discusses should tell a story, ideally the field contamination rate (`theta0`) of a planted contaminated species. Warnings `runOccJSDM()` prints for other blocks get one sentence on the page.
- **Results to cover** (Doug's list, extended in discussion):
  - occupancy at each site, both the predictive probability (what the habitat predicts: `computePredictiveOccupancyProbs()`) and the conditional probability (presence given the detections at that site: `computeConditionalOccupancyProbs()`), each labelled with the question it answers; posterior means only, with each species' baseline occupancy shown with its 95% interval (`plotOccupancyRates()`);
  - environmental covariate effects, including the protected-area contrast, and response curves;
  - trait effects and trait-by-environment effects (JSDM results);
  - collection, PCR detection and false-positive probabilities at both stages, in a figure like Ji et al. 2025 Supplementary Figure 3;
  - residual species associations;
  - ordination and biplot;
  - prediction at new sites;
  - whether more field samples or PCR replicates would help (`plotCumulativeSpeciesDetections()`);
  - tables of all the estimates.
- **Variation partitioning:** include it, shown as one bar per species (measured covariates versus hidden environmental covariates), not the ternary plot, which wastes two of its three sides without a spatial field.
- **Expected species richness per site:** include the mean only, overall and by class and order, computed in the vignette as row sums of the conditional occupancy probabilities over the relevant species. Say plainly that intervals need the package work in `TODO.md`.
- **Do not build any of the TODO items** for this work.

## Claude's calls (Doug did not decide these; check them)

- **Species:** about 20 vertebrates, named by class and order with a running number (for example `Mammalia_Carnivora_1`), so richness by class and order can be computed from the names. No real species names.
- **Planted effects, so each figure tells a story:** Red Listed species favour protected sites more strongly than other species (a planted trait-by-protected effect, echoing Ji et al.'s finding that native species have higher occupancy inside protected areas; all simulated species are native, so Red List status carries the contrast); a planted trait-by-environment effect such as larger-bodied species declining more steeply with human modification; one or two species with noticeable field contamination (`theta0`), echoing the takin in Ji et al. 2025, whose weak southern detections the model discounted as false positives. Effects must be strong enough to show at about 100 sites.
- **Hidden environmental covariates:** `n_factors = 2`, as now.
- **The Supp. Fig. 3-style figure:** per species, collection probability and PCR detection probability (per primer) with 95% intervals, false-positive rates marked on the left margin. Build it from the existing `return...()` outputs with a few lines of ggplot if the existing `plot...()` functions cannot be combined into one panel.
- **The trait-by-environment figures:** `plotTraitsCoefficients(fit, covName = ...)` once per covariate, plus one added figure of each species' estimated slope against its trait value with the fitted trait line, which makes a coefficient tangible.
- **Tables,** each a tidy data frame:
  - site by species occupancy, long format: site, species, predictive and conditional probabilities, and the raw evidence (positive samples out of M, positive PCRs out of all PCRs);
  - expected richness per site, overall and by class and order;
  - per-species detection profile: baseline occupancy, collection probability, detection per primer, field and lab false-positive rates, each with mean and 95% interval;
  - covariate effects: species by coefficient (each non-reference factor level separately), mean, interval, and whether the interval excludes zero;
  - trait-by-environment effects: trait by covariate, same columns;
  - residual associations: only the species pairs whose interval excludes zero;
  - convergence: `returnConvergenceDiagnostics()`.
- **Page outline:**
  1. What occJSDM does (the current opening, adjusted).
  2. The example data: `info`, `OTU` and `traits`, including how the factor sits in `info` and the 0/1 trait in `traits`.
  3. Fit the model: live or precomputed, as now.
  4. Check convergence: R-hat and ESS summaries, the flagged parameters, a traceplot of a flagged parameter beside a healthy one, and what to do (does the parameter matter for your question; run longer and recheck; if it will not settle, the data may not identify it; report it). The text says what a longer run did; no second fit is stored.
  5. Results, in the order of the list above.
  6. Tables.
  7. Where to go next.
- **Carry over Doug's uncommitted edits** to the Quickstart on main: `ncores = 2` in `MCMCparams`, and "hidden environmental covariates (statistically, latent factors)" for `n_factors`.

## Constraints

- **Writing rules:** follow `AGENTS.md`, *Writing the vignettes and lessons*, including the rules from Doug's Lesson 0 review (define every symbol where it first appears; say why each chunk, table or figure is there; one term for one thing; never hard-code a number the data produce, use inline R; check every claim against the rendered numbers). Doug's working copy of `AGENTS.md` on main also holds uncommitted Lesson 1 review rules; read it there. Lesson prose says "hidden environmental covariates", with "latent factors" once per file to explain `n_factors`, and no "factors"/"loadings" otherwise. Do not use `returnLatentPresences()`.
- **Markdown:** the repo's `CLAUDE.md` rules (one line per paragraph, escape `-\>` and `\~` in prose, no em-dashes).
- **Size of `quickstartresults` (Claude's call):** it ships in `data/`, because the Quickstart is the only vignette built into the package. `sampleresults.rda` is already about 50 MB, so slim the new fit: thin the kept draws and drop arrays the page never uses, aiming well under 25 MB. Settle the slimming before the long fit is run.
- **Reproducibility:** a script under `dev/` simulates `quickstartdata` and fits `quickstartresults`, with seeds, so both can be rebuilt. Document both datasets in `R/data.R`. Follow the pattern of `dev/release/2026-10-13-parallel-chains-refit/refit-sampleresults.R`.
- **Simulator:** `simulateOccJSDMData()` draws its own covariates (normal, standard deviation 10; categorical ones with three levels) and coefficients, so it cannot plant named covariates, a two-level factor or chosen effects directly. Either post-process its output or write the simulation in the script; the builder decides.
- **Computing:** R 4.6.1 with Accelerate BLAS. At most 5 fitting workers across all runs on the laptop. Launch long fits detached (nohup), and check power mode before timing. Never `R CMD build` the main checkout; build from a `git archive` export.
- **Review:** Doug reviews the rendered HTML at the end, not the PR diff.

## Left to the builder

- Exact simulation values, chain length (tuned to leave a handful of flags), and object names inside the vignette.
- How to build the richness-by-taxon and Supp. Fig. 3-style figures.
- Which `plot...()` functions need small helper code in the vignette.
