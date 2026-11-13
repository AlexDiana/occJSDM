# Count models for occJSDM

Date: 12 November 2026. Status: a proposal written for Alex's review, at Doug's request. Nothing here is agreed and nothing has been implemented or run. The recommendations are Claude's, made from reading the code on `main` at `1b7ff28`; the open questions in the section *Questions for Alex and Doug* have to be answered before a task-by-task implementation plan is written. Question 1 was answered by Doug on 12 November 2026 and is recorded under *Decided*. This is the `TODO.md` item *Count-data models* under *Later*.

## Decided

- **The data are counts of individuals (Doug, 12 November 2026).** Sequence read counts will not be fitted with this model. For sequence reads the route is a spike-in, as in Alex's eDNAPlus model, which estimates changes in abundance across sites and not absolute abundance. That is the separate `TODO.md` item *Spike-in-based abundance changes*.
- **The first use (Doug, 12 November 2026).** Predict beetles across a landscape from measured covariates and geospatial embeddings, and fit spatial fields to see whether source-sink dynamics can be detected. The dataset is the beetle survey of the embedding study: counts of individuals, several hundred species, a few hundred plots, one trap per plot and so no replicate samples. Claude's reading of what this implies is under *What the first use implies*; it is a recommendation, not a decision.
- **The JSDM sits on abundance: configuration A (Doug, 12 November 2026).** Chosen as the more suitable for the first use, after the comparison under *Where the JSDM structure sits*. The plain count JSDM is built first; it is configuration A with the occupancy layer left out. Configuration B is deferred, not rejected: it stays the natural model for an occupancy analysis made robust to missed detections. This answers question 2 from Doug's side and is open to Alex's review.
- **Taxonomy comes from the grouped-prior work (Doug, 12 November 2026).** The count model adds no taxonomic or phylogenetic structure of its own. Borrowing between related species will come from part 2 of the grouped-prior work on the embedding branch, which changes the prior and not the likelihood. See *Interaction with the embedding covariates*.

## Summary

- `runOccJSDM()` rejects count data today. This proposal adds a negative-binomial count likelihood, with Poisson as a special case, fitted by the same Polya-Gamma augmentation the binary model already uses, so that the existing coefficient, trait, factor and spatial updates run unchanged.
- It serves two uses. First, a zero-inflated count model in which a latent occupancy state separates "absent" from "present but not counted", so that false negatives do not bias the result. Second, a count response that tells the spatial field how much of a species is at each site, where presence and absence told it almost nothing at 100 sites.
- The two uses want the JSDM structure (covariates, traits, latent factors, spatial field) in different places: the first on occupancy, the second on abundance. Doug has chosen abundance, configuration A below, for his first use; Alex is asked whether he sees a statistical objection. The plan builds the plain count JSDM first, which is configuration A without its occupancy layer, and adds that layer only if the data need it.

## Why

- **Robustness to observation error.** A count of zero has two readings: the species is absent, or it is present and the sample missed it. A model with a latent occupancy state and a count distribution for occupied sites gives each zero a probability of being a miss, computed from how large the counts are where the species is seen. This is the count analogue of what the occupancy model already does with repeated binary detections.
- **Information for the spatial field.** The Lesson 7 sweep and the spatial recoverability pilot found that at 100 sites the spatial field recovers the right places but far too faintly, because each site and species reaches the field only as a presence or an absence. `TODO.md` already notes, under *Spatial structure and variation partitioning*, that a quantitative response would give the field more to learn from. A count is the simplest such response. This helps only surveys that count individuals; for eDNA surveys the quantitative response would have to come from the spike-in model, which this proposal does not cover.

## What the code does today

Verified on 12 November 2026 against `1b7ff28`.

- `inferDataModel()` in `R/runOccJSDM.R` labels whole-number data that are not all 0 or 1 as counts, and then stops with "Counts model not supported yet" (line 249). That label is reached only when no site has repeated rows. With repeated rows per site the data are read as occupancy or two-stage data and every count is cut to a detection at `threshold`. So a count model with replicate samples cannot be inferred from the data's shape and must be asked for explicitly.
- Stubs for a count model exist and disagree on its name. `runOccJSDM()` sets `jsdmModel` to "counts" (line 673); `update_jSDMcoef()` and `SD()` in `R/jsdmfun.R` test for "count" (lines 1315 and 1729); `computeModelLoglikJSDM_cpp()` in `src/functions.cpp` tests for "counts" and uses a Poisson density (line 1187). The branch at line 1315 is empty.
- The internal `simulateData()` has a working count branch (`R/jsdmfun.R` line 775) that draws negative-binomial counts with mean `exp(eta)` and one size per species, `rnb`. The exported `simulateOccJSDMData()` does not offer it.
- The retired `sample_rnb()` in `deprecated/R/jsdmfun-maintenance-20261007.R` is a random-walk update of the size on the log scale. It reads the current size from an undefined variable, has no prior and uses a proposal standard deviation of 5.
- Everything in `update_jSDMcoef()` after the weights are drawn is a Gaussian update given a weight matrix `Omega` and a matrix `k`: the collapsed range step, the coefficient and loading draw, the trait regressions and the factor scores. The binary model supplies `k = z - 0.5` and Polya-Gamma weights. Any likelihood that can be written in that form reuses all of it.
- The Polya-Gamma sampler `rpg()` in `src/jsdm.cpp` takes a whole-number shape and sums that many unit draws. A count model needs shapes that are large and not whole numbers.
- The occupancy state is drawn in `sample_z_cpp_parallel()` from the product of the occupancy probability and the probability of the site's sample outcomes under presence and under absence. The collapsed intercept shift in `R/b0_shift.R` uses the same two per-site log-probabilities, `l1` and `l0`. Both generalise to counts by swapping the Bernoulli terms for count terms.
- `computeSiteWAIC()` accepts only the binary, occupancy and two-stage models (`R/site-waic.R` line 121). The covariate-response plots choose between an identity and a logit link (`R/output.R` lines 573 and 624).

## The model

Sites are i, samples within a site are m, species are s. A site may have one sample.

```
z[i,s]            ~ Bernoulli(psi[i,s])            occupancy
y[i,m,s] | z = 1  ~ NegBin(mean = E[i,m] * mu[i,s], size = r[s])
y[i,m,s] | z = 0  = 0

logit(psi[i,s]) = occupancy predictor
log(mu[i,s])    = abundance predictor
```

- `E` is a known effort for each sample (trap-days, hours of observation, area searched), entering as an offset.
- `r` is the size: the variance of a count is its mean plus the squared mean divided by `r`. Large `r` is Poisson.
- Without the occupancy layer the first and third lines drop out and this is an ordinary negative-binomial JSDM.

### Where the JSDM structure sits

The package has one block of JSDM structure: covariates, traits, latent factors and the spatial field, updated together. The model above has two predictors, so there are three configurations.

- **A, on abundance.** The JSDM block models log abundance. Occupancy, if present, is a per-species logistic regression on occupancy covariates, fitted with the sampler that now fits the collection coefficients. This is the configuration for the spatial use: the field sees the counts.
- **B, on occupancy.** The JSDM block stays where it is. The counts replace the binary detections as the observation layer: a per-species count regression on the collection covariates takes the place of the collection probability. This is the smallest change to the package and the configuration for the robustness use. Every existing output keeps its meaning. The spatial field still sees only occupancy, although the occupancy states are better determined than with binary detections.
- **C, on both.** Two JSDM blocks. `update_jSDMcoef()` takes and returns lists, so the sampler could run it twice, but the stored output, diagnostics, plots and WAIC all assume one block. This roughly doubles the surface of the package.

Doug chose A on 12 November 2026 (see *Decided*). The plan therefore builds A without the occupancy layer first and adds A's occupancy layer only if needed. B and C are deferred.

### What the first use implies

The first use (see *Decided*) points to the JSDM on abundance, and possibly to no occupancy layer at all.

- **The source-sink signal is in the numbers.** A sink shows as a few individuals in poor habitat near good habitat. On the abundance scale that is a positive spatial residual where the environmental prediction is low. As presence and absence it is one more presence, and the spatial pilot found presence-and-absence fields too faint to read at a few hundred sites. So the spatial field has to see the counts: the plain count JSDM or configuration A, not B.
- **Prediction comes from the same model.** A count JSDM predicts expected catch per unit effort at any location, and from it the probability of catching at least one individual, which is an occurrence map. Both come with posterior draws.
- **One trap per plot makes the occupancy layer hard to identify** (see *Cautions*). The plain count JSDM already treats a zero at a site of low expected abundance as a likely miss and not as an absence, which is the false-negative behaviour wanted where a species is scarce. The occupancy layer earns its place only if the fitted model predicts too few zeros, which a posterior predictive check on the zero counts will show.
- **Embeddings must then work on log abundance.** This turns the scale check in *Interaction with the embedding covariates* from optional into required work, and reverses the recommendation in question 15.
- **Many rare species.** With several hundred species, many seen a handful of times, a separate size and an unpooled intercept per species will be poorly determined. This bears on questions 5 and 8.

Recommendation for this use: stages 1 and 2, then the embedding scale check and the spatial simulation, and the occupancy layer in configuration A only if the zero check asks for it. Configuration B is not needed for this use.

## How it is fitted

- **Negative binomial by Polya-Gamma augmentation** (Pillow and Scott 2012; Polson, Scott and Windle 2013). Write `phi = eta + o`, where `eta` is the log-mean predictor and the offset `o` is the log effort minus the log size. The likelihood of a count `y` is proportional to `exp(phi)^y / (1 + exp(phi))^(y + r)`. With a weight drawn from a Polya-Gamma distribution of shape `y + r` at `phi`, it becomes Gaussian in `eta` with `Omega` equal to the weight and `k` equal to `(y - r)/2` minus the weight times the offset. That is exactly the form the binary branch already consumes, so the range step, the coefficient draw and the factor-score draw need no new algebra.
- **Replicate samples** add up: a site's `k` and `Omega` are the sums over its samples, because all samples at a site share the site's predictor and differ only by the known offset. A cell with `z = 0` contributes nothing, so its `k` and `Omega` are zero. The JSDM block therefore always sees a site-by-species matrix, as now.
- **The size** is updated per species by a random walk on the log scale with an adaptive step during burn-in and a proper prior, replacing `sample_rnb()`. The predictor stays the log mean, which matches the simulator's existing parameterisation.
- **Poisson** is offered as the negative binomial with the size fixed at a large value. This is an approximation: at a size of 1,000 a mean of 1,000 has twice the Poisson variance. See question 4.
- **The Polya-Gamma sampler** needs a version for any positive shape. Proposed: the existing exact sampler for the whole-number part up to a cut-off, a truncated sum of gammas for the fractional part, and above the cut-off a normal draw with the exact Polya-Gamma mean and variance, as `BayesLogit` does above a shape of 170. Written in C++ on the package's own random stream, to keep fits reproducible across thread counts.
- **The occupancy state.** A site with any positive count is occupied. A site whose counts are all zero is occupied with probability proportional to the occupancy probability times the product over its samples of the probability of a zero count. The collapsed intercept shift uses the same two quantities.

### Approaches considered

- **Polya-Gamma negative binomial (recommended).** Reuses the whole Gaussian block. Costs: a general-shape sampler, a Metropolis step for the size, and Poisson only as a limit.
- **Poisson-lognormal with a latent log intensity.** Each cell gets a latent log intensity, the JSDM runs in its existing continuous mode on it, and the residual standard deviation plays the part of overdispersion. Exact Poisson, and nearly no change to the JSDM block. Costs: one Metropolis update per cell per iteration, slow mixing when overdispersion is small, and a site-level likelihood that needs numerical integration. Its structure is the same as the latent log biomass in the eDITH draft, so if Qin builds that machinery the two could share it. See question 3.
- **Conjugate size update** (Zhou and others 2012). Puts the regression on the logit of the negative-binomial probability, so the size has a gamma full conditional. Rejected for now because the intercept then stops being the log mean and trades off against the size.
- **Metropolis or Hamiltonian updates of the coefficients.** Rejected: it discards the conjugate structure the package is built on, including the collapsed range step.

## Cautions

- **Zero inflation is weakly identified without replicates.** With one sample per site, absence and overdispersion both produce extra zeros, and only the assumed shape of the count distribution separates them. A negative binomial with a small size can absorb most zeros, and a Poisson reads any unmodelled overdispersion as absence. With replicate samples the two separate properly, because an absent site is zero in every sample while misses are independent between samples. This is the same reason occupancy models need repeat visits.
- **Latent factors on log abundance also create overdispersion**, so the factors, the size and the zero inflation compete for the same variation. Recovery has to be checked with all three present.
- **Embeddings and the spatial field compete.** Embeddings are numerous and vary smoothly in space, so they can absorb the pattern the spatial field is meant to find, and the field can take variation that belongs to the embeddings. A source-sink reading needs the environmental part of the model to be the estimate of habitat quality and the field to be what habitat does not explain. How well the two separate has to be measured by simulation before any field is interpreted. This is the count-model case of the `TODO.md` item *Keep environmental effects stable when spatial or latent-factor terms are added*.
- **A spatial residual is evidence of a pattern, not of demography.** A positive field in poor habitat is what a sink would produce, and also what a missing, spatially smooth covariate would produce. Traps that intercept flying insects also catch individuals in transit. The model can show where numbers exceed what habitat predicts; calling that a sink needs an argument from outside the model.
- **Sequence read counts must not be fitted with this model.** They are not counts of organisms: they are compositional within a sequencing run and carry species-specific amplification bias, so their size says nothing absolute about abundance. This is why the model is restricted to counts of individuals (see *Decided*). The hazard that remains is practical: a table of read counts is also a table of whole numbers, so nothing in the data stops a user from fitting one. The help page should say so, and the model should never be chosen automatically from the data's type (question 10).

## Proposed interface

Open to change; see question 10.

- A new argument to `runOccJSDM()`, `countModel`, `NULL` by default. When `NULL`, behaviour is as now, except that the "not supported" error tells the user how to ask for a count model.
- When set, a list: `family` ("negbin" or "poisson"), `occupancy` (`FALSE` for a plain count JSDM, `TRUE` to add the occupancy layer) and `offset` (the name of an effort column in `data$info`, on the natural scale). With configuration A chosen, the JSDM is always on abundance and the list needs no element to say so. A `jsdm` element, "abundance" or "occupancy", would be added only if configuration B is built later.
- With `countModel` set, repeated rows per site are replicate counts and `threshold` is not applied.
- New stored output: the size draws, with their convergence diagnostics; the abundance-scale predictor; for the occupancy layer, the occupancy state and probability as the occupancy model stores them now.
- `simulateOccJSDMData()` gains a count model, and its `rnb` input, documented today as having no effect, becomes live.

### Example calls

A presence and absence or eDNA fit is written exactly as it is today:

```r
fit <- runOccJSDM(data, occCovariates = c("elevation", "canopy"))
```

A plain count JSDM, with no occupancy layer:

```r
fit <- runOccJSDM(
  data,
  occCovariates = c("elevation", "canopy"),
  countModel = list(family = "negbin", offset = "trap_days")
)
```

Configuration A with its occupancy layer adds `occupancy = TRUE` to the call above; how its occupancy covariates are named is question 14.

Deferred, for comparison only: a zero-inflated count model in configuration B, with the JSDM on occupancy and the counts as the observation layer, would look like this:

```r
fit <- runOccJSDM(
  data,
  occCovariates  = c("elevation", "canopy"),
  collCovariates = c("trap_type"),
  countModel = list(family = "negbin", occupancy = TRUE,
                    jsdm = "occupancy", offset = "trap_days")
)
```

### What the covariate arguments mean

- In configuration B the existing names fit: `occCovariates` act on occupancy and `collCovariates` act on the counts, in the place of the collection probability.
- In the plain count JSDM and in configuration A the JSDM block models abundance, so `occCovariates` act on log abundance despite the name. Configuration A with the occupancy layer also needs somewhere to name the covariates of the per-species occupancy regression, which no existing argument provides. See question 14.

### Keeping it out of the way of presence and absence users

The feature follows the rules Doug set for the embedding prior on 12 November 2026.

- It is off unless `countModel` is set. With `countModel = NULL` no fit changes, and stage 2 tests that existing fits are identical draw for draw.
- It is never chosen from the data's type, so data of 0 and 1 cannot reach the count code.
- It has its own help page. `?runOccJSDM` gains the argument and one line pointing there.
- The extra stored output exists only in count fits.
- The lessons and the Quickstart do not change.
- The one change everyone sees is the wording of the existing "not supported" error.
- The count help page must warn about one existing behaviour: without `countModel`, counts with repeated rows per site are cut to detections at `threshold`, silently.

## Interaction with the embedding covariates

Checked on 12 November 2026 against the branch `claude/grouped-covariate-prior` at `d2c19bb`, which is not yet merged. That feature uses large numbers of geospatial embedding columns as occupancy covariates under a grouped prior, switched on by naming `embeddingVars`.

- **The statistics do not interfere.** The grouped prior changes the prior on the coefficients. The count model changes the likelihood, which is how `k` and `Omega` are formed. Both feed the same Gaussian update and neither changes the other's inputs.
- **The code will conflict if built in parallel.** Both features edit `update_jSDMcoef()` and the same files: `R/jsdmfun.R`, `R/chains.R`, `R/runOccJSDM.R`, `R/simulateData.R`, `R/output.R` and `src/jsdm.cpp`. The embedding branch changes about 1,600 lines there and is close to done. Proposed order: merge the embedding branch first, then start stage 2 from `main`. Stage 1 adds only new functions and can proceed at any time.
- **Configuration B leaves embeddings as they are**: occupancy covariates on the logit scale.
- **The plain count JSDM and configuration A put embeddings on log abundance.** The grouped prior's default axis scale was calibrated for a contribution of 1 to 2 logits and has not been checked on the log scale. The first use needs exactly this combination, so the check is required work: repeat the embedding branch's prior check with a log link and choose a default axis scale for count fits. Until it is done, `countModel` with the JSDM on abundance stops with a clear error when `embeddingVars` is set. See question 15.
- **Taxonomy depends on that branch.** Count fits get taxonomic borrowing only from part 2 of the grouped-prior work (see *Decided*). Until part 2 exists, rare species in a count fit borrow strength through the trait regression on their coefficients, which the plain count JSDM and configuration A keep, and through whatever is decided for the intercepts and the size in questions 5 and 8. When part 2 is built, it has to be tested on a count fit as well as on the binary and two-stage fits, because its scales, like the embedding prior's, will have been set on the logit scale.
- **A guard that already exists.** The embedding branch ships a reference fit that must be reproduced when its prior is off. The count work must keep that test passing.

## Questions for Alex and Doug

1. **(Doug) Which data come first? Answered 12 November 2026:** counts of individuals (for example sorted trap catches, point counts, camera-trap detections). Sequence read counts are excluded; see *Decided*. Consequences: the offset is sampling effort, the test data are simulated counts of individuals, and question 6 carries less weight.
2. **Where should the JSDM structure sit: A, B or C?** The recommendation is A without occupancy first, then one of A or B with occupancy. If the robustness use matters most, B is the smaller and safer step. If the spatial use matters most, A with occupancy is the target. Doug's first use (see *What the first use implies*) is the spatial one, so the recommendation for now is the plain count JSDM, then A with occupancy only if needed, with B deferred. **Doug's answer, 12 November 2026: A.** For Alex: is there a statistical objection to A, in particular to an occupancy layer that is a plain per-species regression beside an abundance JSDM?
3. **Polya-Gamma negative binomial or Poisson-lognormal?** In particular, does the latent log-biomass work planned for the eDITH fork make the Poisson-lognormal route cheaper overall, or is it better to keep counts independent of that work?
4. **Is Poisson as a large fixed size acceptable**, or should Poisson be left out of the first version, or be exact through the Poisson-lognormal route?
5. **The size.** One per species, as the simulator has it, or a shared or partially pooled size? Which prior: a gamma on the size, or a half-normal on its inverse square root? The first use has several hundred species, many of them rare, which favours partial pooling.
6. **Counts at unoccupied sites.** Exactly zero, as written above, so that any positive count fixes the site as occupied? Or a small false-positive count rate per species, the analogue of `theta0`? With counts of individuals a positive count at an unoccupied site can only be a misidentification, not contamination, so the recommendation is exactly zero in the first version.
7. **Zero inflation with one sample per site.** Allow it with a warning, or require replicates? The first use has one trap per plot, so this decides whether the occupancy layer is available to it at all.
8. **The abundance intercepts.** Fixed prior, or the community prior that occupancy intercepts now use by default? The many rare species of the first use favour the community prior.
9. **The Polya-Gamma approximation.** Is a moment-matched normal draw above a shape cut-off acceptable in a sampler otherwise exact, and is 170 the right cut-off?
10. **(Doug and Alex) Interface.** One `countModel` list, as proposed, or separate arguments? Should integer data without site replication fit a count model automatically, as the current inference half-intends? The recommendation is no: a table of sequence reads is also whole numbers, and automatic inference would fit it as counts of individuals.
11. **Site-level WAIC.** Leave counts out of `computeSiteWAIC()` at first, with a clear error, and keep only the complete-data WAIC?
12. **Variation partitioning.** The stub takes standard deviations of `exp(eta)`, which the sites with the largest means will dominate. Partition on the log scale instead?
13. **(Doug) Where and when.** Pull requests to `main` from this branch, one per stage, or in Qin's fork? Before or after the beta work in `TODO.md`?
14. **(Doug and Alex) Covariate argument names.** When the JSDM block models abundance, should `occCovariates` simply be documented as the covariates of whatever the JSDM block models, or should `countModel` carry its own covariate names? And where are the occupancy-regression covariates named in configuration A with the occupancy layer?
15. **Embeddings with the JSDM on abundance.** The first use needs it, so the recommendation is now to recheck the grouped prior's scale on the log scale and support the combination, as stage 3. Is the half-Cauchy axis prior still the right form with a log link, where a large coefficient multiplies the expected count instead of saturating a probability?

## Staged plan

Each stage is one pull request to Alex and leaves `main` working. Stages 1 and 2 do not depend on questions 2, 6, 7 or 13. Stage 2 starts from `main` after the embedding branch has merged; see *Interaction with the embedding covariates*. A task-by-task plan for each stage is written in `dev/superpowers/plans/` once its questions are answered.

- **Stage 1: likelihood pieces, tested alone.** The general-shape Polya-Gamma sampler, checked against the exact mean and variance across shapes and against the existing sampler at whole-number shapes. A negative-binomial log-likelihood with offsets. The size update with its prior, checked on simulated counts with a known predictor. One name for the model, "counts", in R and C++. No change to any fit.
- **Stage 2: the plain count JSDM.** Wire the count branch into `update_jSDMcoef()`; accept `countModel`; offsets and replicate samples; size draws in the output, the posterior-draw helpers and the diagnostics; a log link in the covariate-response plots; the complete-data WAIC with the negative-binomial density; `simulateOccJSDMData()` with a count model; a posterior predictive check that compares observed and predicted numbers of zeros per species, which stage 5 depends on. Done when a simulate-and-recover test passes for coefficients, size and residual correlations, existing fits of the other four models are unchanged draw for draw, and a comparison with an independent negative-binomial JSDM (Hmsc or gllvm) on one simulated dataset agrees within Monte Carlo error.
- **Stage 3: embeddings on log abundance.** Repeat the embedding branch's prior check with a log link, choose the default axis scale for count fits, and lift the error that stage 2 raises for the combination. Done when a simulate-and-recover test with an embedding block passes for a count fit and the prior's implied contribution to the predictor is reported on the log scale.
- **Stage 4: does the spatial field gain from counts, and can it show a sink?** A simulation, not a feature. Part one: on communities from the spatial recoverability pilot's design, fit the same sites as counts and as presence and absence and compare recovered amplitude, range and held-out field. Part two, the `TODO.md` item *Source-sink inference scenario*: simulate species whose environmental and spatial effects oppose each other, at the first use's number of plots, with and without an embedding block, and ask whether the fitted field is positive where the environmental prediction is low. Part two measures the competition between embeddings and the field described under *Cautions*. This stage decides whether the first use's spatial aim is achievable. Run detached, at most five workers, power mode recorded.
- **Stage 5: the occupancy layer, if needed.** Built only if a posterior predictive check on the zero counts of a plain count fit shows too few predicted zeros, or if another use needs it. In configuration A: the count version of the occupancy draw; the per-species occupancy regression, fitted with the sampler that now fits the collection coefficients; zero weight in the abundance block for cells drawn as unoccupied, after checking that no update divides by `Omega`; the simulator's zero-inflated case. The collapsed intercept shift is not involved, because it belongs to an occupancy JSDM. Done when recovery and interval coverage are reported at 1, 3 and 5 samples per site, so that the identifiability caution above is measured and not assumed, and when the result at one sample per site settles question 7.
- **Stage 6: documentation.** `?runOccJSDM`, `NEWS.md`, the README line "Count responses are not supported", and a lesson section once Doug has a dataset to teach from.

## Non-goals

- No model of sequence read counts, with or without spike-ins; that remains the eDNAPlus item, and it estimates changes in abundance, not absolute abundance.
- No two-block configuration C in this plan.
- No sample-level covariates on abundance in configuration A beyond the offset.
- No taxonomic or phylogenetic structure; that belongs to part 2 of the grouped-prior work.
- No change to the binary, continuous, occupancy or two-stage models, their defaults or their stored output.

## References

- Pillow and Scott (2012), Fully Bayesian inference for neural models with negative-binomial spiking, *Advances in Neural Information Processing Systems* 25.
- Polson, Scott and Windle (2013), Bayesian inference for logistic models using Polya-Gamma latent variables, *Journal of the American Statistical Association* 108.
- Windle, Polson and Scott (2014), Sampling Polya-Gamma random variates: alternate and approximate techniques, arXiv 1405.0506.
- Zhou, Li, Dunson and Carin (2012), Lognormal and gamma mixed negative binomial regression, *Proceedings of the 29th International Conference on Machine Learning*.
- Diana and others (2024), eDNAPlus, *Journal of the American Statistical Association*, as cited in `TODO.md`.

These citations are from memory and have not been checked against Zotero.
