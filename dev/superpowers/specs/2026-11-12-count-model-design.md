# Count models for occJSDM

Date: 12 November 2026, revised 15, 16, 23 and 25 November 2026. Status: a design under review by Alex, written at Doug's request. The model's structure is agreed between Doug and Alex (see *Decided*); the remaining questions in *Questions for Alex and Doug* have to be answered before a task-by-task implementation plan is written. Nothing has been implemented or run. The recommendations are Claude's, made from reading the code on `main` at `1b7ff28`. This is the `TODO.md` item *Count-data models* under *Later*.

## Decided

- **The data are counts of individuals (Doug, 12 November 2026).** Sequence read counts will not be fitted with this model. For sequence reads the route is a spike-in, as in Alex's eDNAPlus model, which estimates changes in abundance across sites and not absolute abundance. That is the separate `TODO.md` item *Spike-in-based abundance changes*.
- **The first use (Doug, 12 November 2026).** Predict beetles across a landscape from measured covariates and geospatial embeddings, and fit spatial fields to see whether source-sink dynamics can be detected. The dataset is the beetle survey of the embedding study: counts of individuals, several hundred species, a few hundred plots, one trap per plot and so no replicate samples. Each trap was emptied once, as far as Doug recalls (16 November 2026), so there are no repeat collections in time either.
- **One JSDM with a shared predictor: configuration D (proposed by Alex and adopted by Doug, 15 November 2026).** A single JSDM predictor gives the occupancy probability on the logit scale and, through two further parameters, the expected count where the species is present. This replaces Doug's choice of 12 November, configuration A, in which the JSDM modelled abundance and occupancy was a plain regression. The configurations are compared under *Configurations considered*.
- **The plain count JSDM is built first.** It has no occupancy layer. It is the first step towards D, because D needs every piece of it, and it stays in the package as the comparison model.
- **Taxonomy comes from the grouped-prior work (Doug, 12 November 2026).** The count model adds no taxonomic or phylogenetic structure of its own. Borrowing between related species will come from part 2 of the grouped-prior work on the embedding branch, which changes the prior and not the likelihood. See *Interaction with the embedding covariates*.

## Summary

- `runOccJSDM()` rejects count data today. This design adds a negative-binomial count likelihood, with Poisson as a special case, fitted by the same Polya-Gamma augmentation the binary model already uses, so that the existing coefficient, trait, factor and spatial updates run unchanged.
- It serves two uses. First, a latent occupancy state separates "absent" from "present but not counted", so that false negatives do not bias the result. Second, the counts tell the spatial field how much of a species is at each site, where presence and absence told it almost nothing at 100 sites.
- Both uses are met by one model. The JSDM predictor is the occupancy predictor, as it is now, and the expected count where present rises with it. Occupancy output keeps its meaning, and the counts inform every part of the JSDM.
- The work is staged: the likelihood pieces, the plain count JSDM, the shared-predictor model, then a simulation that asks whether the spatial aim of the first use is achievable.

## Why

- **Robustness to observation error.** A count of zero has two readings: the species is absent, or it is present and the sample missed it. A model with a latent occupancy state and a count distribution for occupied sites gives each zero a probability of being a miss, computed from how large the counts are where the species is seen. This is the count analogue of what the occupancy model already does with repeated binary detections.
- **Information for the spatial field.** The Lesson 7 sweep and the spatial recoverability pilot found that at 100 sites the spatial field recovers the right places but far too faintly, because each site and species reaches the field only as a presence or an absence. `TODO.md` already notes, under *Spatial structure and variation partitioning*, that a quantitative response would give the field more to learn from. A count is the simplest such response. This helps only surveys that count individuals; for eDNA surveys the quantitative response would have to come from the spike-in model, which this design does not cover.

## What the code does today

Verified on 12 November 2026 against `1b7ff28`.

- `inferDataModel()` in `R/runOccJSDM.R` labels whole-number data that are not all 0 or 1 as counts, and then stops with "Counts model not supported yet" (line 249). That label is reached only when no site has repeated rows. With repeated rows per site the data are read as occupancy or two-stage data and every count is cut to a detection at `threshold`. So a count model with replicate samples cannot be inferred from the data's shape and must be asked for explicitly.
- Stubs for a count model exist and disagree on its name. `runOccJSDM()` sets `jsdmModel` to "counts" (line 673); `update_jSDMcoef()` and `SD()` in `R/jsdmfun.R` test for "count" (lines 1315 and 1729); `computeModelLoglikJSDM_cpp()` in `src/functions.cpp` tests for "counts" and uses a Poisson density (line 1187). The branch at line 1315 is empty.
- The internal `simulateData()` has a working count branch (`R/jsdmfun.R` line 775) that draws negative-binomial counts with mean `exp(eta)` and one size per species, `rnb`. The exported `simulateOccJSDMData()` does not offer it.
- The retired `sample_rnb()` in `deprecated/R/jsdmfun-maintenance-20261007.R` is a random-walk update of the size on the log scale. It reads the current size from an undefined variable, has no prior and uses a proposal standard deviation of 5.
- Everything in `update_jSDMcoef()` after the weights are drawn is a Gaussian update given a weight matrix `Omega` and a matrix `k`: the collapsed range step, the coefficient and loading draw, the trait regressions and the factor scores. The binary model supplies `k = z - 0.5` and Polya-Gamma weights. Any likelihood that can be written in that form reuses all of it.
- The Polya-Gamma sampler `rpg()` in `src/jsdm.cpp` takes a whole-number shape and sums that many unit draws. A count model needs shapes that are large and not whole numbers.
- The occupancy state is drawn in `sample_z_cpp_parallel()` from the product of the occupancy probability and the probability of the site's sample outcomes under presence and under absence. The collapsed intercept shift in `R/b0_shift.R` uses the same two per-site log-probabilities, `l1` and `l0`.
- `computeSiteWAIC()` accepts only the binary, occupancy and two-stage models (`R/site-waic.R` line 121). The covariate-response plots choose between an identity and a logit link (`R/output.R` lines 573 and 624).

## The model

Sites are i, samples within a site are m, species are s. A site may have one sample.

```
z[i,s]            ~ Bernoulli(psi[i,s])            occupancy
y[i,m,s] | z = 1  ~ NegBin(mean = E[i,m] * mu[i,s], size = r[s])
y[i,m,s] | z = 0  = 0

logit(psi[i,s]) = eta[i,s]
log(mu[i,s])    = tau_0[s] + tau[s] * eta[i,s]        tau positive
```

- `eta` is the JSDM predictor: the species intercept, the covariates, the latent factors and the spatial field, with the trait regression on the coefficients. It is the occupancy predictor, exactly as in the package today. Alex confirmed on 23 November 2026 that the shared predictor is the whole of it (question 16).
- `tau_0` sets the level of abundance where a species is present. `tau` sets how steeply that abundance rises with the occupancy predictor, and is positive, so that a site more likely to be occupied also holds more individuals when it is. Both are per species to begin with (Alex, 23 November 2026; question 17).
- `E` is a known effort for each sample (trap-days, hours of observation, area searched), entering as an offset.
- `r` is the size: the variance of a count is its mean plus the squared mean divided by `r`. Large `r` is Poisson.

The **plain count JSDM**, built first, has no occupancy state: every count is negative binomial with `log(mu)` equal to `eta`, so there the JSDM predictor is the log of the expected count.

### What the shared predictor gives

- **Both uses from one JSDM.** Occupancy keeps the full JSDM, so every existing occupancy output keeps its meaning: occupancy maps, species associations in occurrence, trait effects, the spatial field. The counts feed the same predictor, so the field, the factors and the coefficients all learn from the numbers.
- **One set of coefficients.** There is no second block of coefficients, loadings and spatial coefficients, so the stored output, diagnostics and plots do not double.
- **The predictor stays on the logit scale.** `occCovariates` keeps its meaning, and the grouped embedding prior keeps its calibration.
- **A coherent reading of zeros with one sample per site.** The model adds two parameters per species, not a second regression. A zero in good habitat cannot be explained away separately from the rest of the model: it is an unlikely absence, an unlikely miss, or a reason for the shared predictor, including its spatial field, to be lower there.
- **Predictions.** For any location: the occupancy probability; the expected count where present; their product, the expected count; and the probability of counting at least one individual at a standard effort. All come with posterior draws.

### What it costs

- **Occupancy and abundance must move together.** Within a species, every covariate, the spatial field and the factors act on abundance in the same direction as on occupancy and in the same proportion. A species that is widespread but sparse in one habitat and local but dense in another cannot be represented, and neither can a covariate that changes numbers without changing occurrence. A positive relationship between occupancy and abundance is common, so this is a reasonable first model, but it is an assumption. Stage 4 checks what happens when it is false, and held-out predictions can be compared with the plain count JSDM.
- **The scale of the predictor may be weakly determined.** Abundance sees only the product of `tau` and `eta`. Occupancy alone fixes the scale of `eta`, and with one trap per plot the occupancy state is mostly inferred from the counts. Where the count distribution explains the zeros by itself, a larger `eta` with a smaller `tau` fits almost as well, and the priors on the coefficients then set the scale. See question 18.
- **Mixing.** `tau` multiplies the whole predictor, so it will be strongly correlated with the coefficients in the posterior. Alternating draws may move slowly along that direction. Alex has proposed a joint update for this; see *The interweaving step*.
- **The plain count JSDM is not a special case.** Occupancy near certainty needs a large `eta`, which also forces a large mean. The two are different models, which is why the plain count JSDM is kept as a comparison.

## How it is fitted

- **Negative binomial by Polya-Gamma augmentation** (Pillow and Scott 2012; Polson, Scott and Windle 2013). Write `phi` for the log mean plus an offset `o`, the log effort minus the log size. The likelihood of a count `y` is proportional to `exp(phi)^y / (1 + exp(phi))^(y + r)`. With a weight drawn from a Polya-Gamma distribution of shape `y + r` at `phi`, it becomes Gaussian in `phi`.
- **The plain count JSDM.** Here `phi = eta + o`, so `Omega` is the weight and `k` is `(y - r)/2` minus the weight times the offset. That is the form the binary branch already consumes.
- **The shared-predictor model.** Each cell has an occupancy term and, where occupied, a count term. Given `tau_0` and `tau`, both are Gaussian in `eta` once their weights are drawn, as Alex stated; Claude checked the algebra and confirmed it numerically.

```
occupancy term, every cell:    kappa1 = z - 1/2      omega1 from PG(1, eta)
count term, cells with z = 1:  kappa2 = (y - r)/2    omega2 from PG(y + r, phi)
                               phi = tau_0 + tau * eta + o

k     = kappa1 + z * tau * (kappa2 - omega2 * (tau_0 + o))
Omega = omega1 + z * tau^2 * omega2
```

- **The JSDM block is unchanged.** It receives `k` and `Omega` as it does now, so the range step, the coefficient draw, the trait regressions and the factor scores need no new algebra.
- **`tau_0` and `tau`.** Given `eta`, the pair has a bivariate normal full conditional for each species: a weighted regression of the count term's working response on `eta` over the occupied cells, with `tau` truncated to positive values.
- **Replicate samples** add up: the count term is summed over a site's samples, because all samples at a site share the site's predictor and differ only by the known offset. The JSDM block therefore always sees a site-by-species matrix, as now.
- **The size** is updated per species by a random walk on the log scale with an adaptive step during burn-in and a proper prior, replacing `sample_rnb()`.
- **Poisson** is offered as the negative binomial with the size fixed at a large value. This is an approximation: at a size of 1,000 a mean of 1,000 has twice the Poisson variance. See question 4.
- **The Polya-Gamma sampler** needs a version for any positive shape. Proposed: the existing exact sampler for the whole-number part up to a cut-off, a truncated sum of gammas for the fractional part, and above the cut-off a normal draw with the exact Polya-Gamma mean and variance, as `BayesLogit` does above a shape of 170. Written in C++ on the package's own random stream, to keep fits reproducible across thread counts.
- **The occupancy state.** A site with any positive count is occupied. A site whose counts are all zero is occupied with probability proportional to the occupancy probability times the product over its samples of the probability of a zero count.
- **The interweaving step and the collapsed intercept shift** are two further updates, each with its own section below.

### The interweaving step

Alex proposed this on 23 November 2026, in answer to question 18.

- **The problem it addresses.** For a species present at every site, occupancy says nothing about the scale of `eta`, so `eta` and `tau` cannot be told apart. That is the extreme case. With one sample per site a milder form applies to every species, because the occupancy state is unknown wherever the count is zero.
- **Alex's sketch.** (a) Update the coefficients, meaning the whole predictor, given `tau`. (b) Update `tau` given the predictor. (c) Change variables to `beta_tilde`, the coefficients multiplied by `tau`, and update `tau` again with `beta_tilde` held fixed.
- **What step (c) involves.** This is Claude's working, which Alex confirmed on 25 November 2026. Holding `beta_tilde` fixed holds the product of `tau` and `eta` fixed, so the expected counts do not change and the count term drops out. `tau` is then informed by the occupancy term, whose predictor is that fixed product divided by `tau`, and by the priors on the coefficients, evaluated at `beta_tilde` divided by `tau`, with the Jacobian of the change of variables. The result is not a standard distribution, so step (c) is a one-dimensional Metropolis update for each species; Alex agrees that one is needed.
- **What is rescaled.** The factor scores are shared between species and cannot be rescaled for one of them. The rescaling applies to a species' own parts of the predictor: its intercept, coefficients, loadings and spatial coefficients. The prior on each enters step (c), including the community prior on the intercept, the trait regression and, with an embedding block, the grouped prior.
- **A link to the intercept shift.** Step (c) changes the occupancy probabilities and leaves the expected counts alone, as option 3 of the intercept shift does. Both could therefore be done with the occupancy state summed out, at the same low cost.
- **Stage 3 builds it and measures it**, reporting the effective sample size of `tau` with and without step (c).

### The collapsed intercept shift in count fits

This is question 19, set out in full because Alex asked for a clearer account on 23 November 2026.

- **What the step does today.** With the community prior on the occupancy intercepts, which is the default, every species' intercept `B0` and their community mean `mu_b0` drift together slowly. That common direction is informed only through the occupancy states `z`, and the sampler changes those one site at a time. The step in `R/b0_shift.R` proposes adding one number, `delta`, to every `B0` and to `mu_b0`, and accepts or rejects it by Metropolis with `z` summed out. It was added after `mu_b0` failed its convergence check in fits of 100 species.
- **What it relies on.** With `z` summed out, each site and species contributes `psi * L1 + (1 - psi) * L0`. `L1` is the probability of the site's data if the species is present, and `L0` if it is absent. In an occupancy fit both come from the detection part of the model and do not involve the intercepts. A shift therefore changes `psi` only, and `L1` and `L0` are computed once per iteration and reused for every proposal.
- **Why that fails under the shared predictor.** Here `L1` is the negative-binomial probability of the site's counts, and their expected value is `E * exp(tau_0 + tau * eta)`. The intercept is inside `eta`, so adding `delta` to it multiplies the expected count by `exp(tau * delta)`. `L1` is no longer fixed, and the step as written would compute the wrong acceptance ratio.
- **Option 1: recompute.** Evaluate `L1` at the shifted predictor for every proposal. This is exact, and the cost stays of the same order. But the proposal then changes every species' expected count at once, which the counts will mostly reject, so the step may rarely be accepted.
- **Option 2: switch the step off** in count fits and check whether `mu_b0` still mixes. This is the simplest, but the slow drift may return, because `z` is still unknown wherever a count is zero.
- **Option 3: shift `tau_0` to compensate.** This is Claude's suggestion and is untested. Propose `delta` on every `B0` and on `mu_b0`, and at the same time subtract `tau * delta` from each species' `tau_0`. Then `tau_0 + tau * eta` does not change, so `L1` is fixed again and the step is the existing one, with one addition to the acceptance ratio: the prior on `tau_0` at its new value. The move changes the occupancy probabilities and leaves the expected counts where a species is present alone. That is the direction the counts cannot inform and only `z` can, which is the direction the step was built for.
- **Alex's view (25 November 2026).** His guess is that the counts now inform `B0` more directly, so the slow mixing may not appear in count fits at all. He would try option 3 in any case, and asked whether options 1 and 3 can be tried together.
- **Options 1 and 3 together.** This is possible. Each is a valid Metropolis step by itself, so both can run in the same iteration. They move in different directions: option 3 changes the occupancy probabilities and holds the expected counts, and option 1 changes both. Option 1 costs a little more, because `L1` is recomputed for each proposal, and each move adapts its own step size.
- **What stage 3 builds and compares.** Both moves, each with its own switch. Three arms: no shift, which tests Alex's guess; option 3 alone; and options 1 and 3 together. The arms are judged on the convergence and effective sample size of `mu_b0` and on the acceptance rate of each move. The default for count fits is chosen from that comparison.

### Approaches considered for the likelihood

- **Polya-Gamma negative binomial (chosen).** Reuses the whole Gaussian block. Costs: a general-shape sampler, a Metropolis step for the size, and Poisson only as a limit.
- **Poisson-lognormal with a latent log intensity.** Each cell gets a latent log intensity, the JSDM runs in its existing continuous mode on it, and the residual standard deviation plays the part of overdispersion. Exact Poisson. Costs: one Metropolis update per cell per iteration, slow mixing when overdispersion is small, and a site-level likelihood that needs numerical integration. It does not combine with the shared predictor as simply. Its structure is the same as the latent log biomass in the eDITH draft. See question 3.
- **Conjugate size update** (Zhou and others 2012). Puts the regression on the logit of the negative-binomial probability, so the size has a gamma full conditional. Rejected for now because the intercept then stops being the log mean and trades off against the size.
- **Metropolis or Hamiltonian updates of the coefficients.** Rejected: it discards the conjugate structure the package is built on, including the collapsed range step.

## Configurations considered

The package has one block of JSDM structure, and a count model with an occupancy state has two predictors. Four ways of connecting them were considered.

- **A, the JSDM on abundance.** Occupancy is a per-species logistic regression on covariates, a switch that excuses zeros. Doug's choice on 12 November, because the spatial field sees the counts. Set aside for D: A's occupancy has no traits, associations or spatial field, its switch competes with the spatial field for zeros in good habitat when there is one sample per site, and it moves the embeddings to the log scale.
- **B, the JSDM on occupancy.** The counts replace the binary detections as the observation layer, through a per-species count regression. Set aside because the spatial field would still see only occupancy.
- **C, a JSDM on each.** Two blocks of coefficients, loadings and spatial coefficients. Set aside because it roughly doubles the stored output, diagnostics and plots, and, as Alex noted, repeats every coefficient.
- **D, one JSDM with a shared predictor (adopted).** Alex's proposal of 15 November. It gives what B gives for occupancy and what A gives the spatial field, at the cost of tying the two responses together.

## Cautions

- **Occupancy is weakly identified without replicates.** With one sample per site, absence and overdispersion both produce extra zeros, and only the assumed shape of the count distribution separates them. A negative binomial with a small size can absorb most zeros. With replicate samples the two separate properly, because an absent site is zero in every sample while misses are independent between samples. The shared predictor reduces the problem, because occupancy is not free to vary apart from abundance, but it does not remove it. Stage 3 measures it.
- **Latent factors also create overdispersion**, so the factors, the size and the occupancy state compete for the same variation. Recovery has to be checked with all three present.
- **Embeddings and the spatial field compete.** Embeddings are numerous and vary smoothly in space, so they can absorb the pattern the spatial field is meant to find, and the field can take variation that belongs to the embeddings. A source-sink reading needs the environmental part of the model to be the estimate of habitat quality and the field to be what habitat does not explain. How well the two separate has to be measured by simulation before any field is interpreted. This is the count-model case of the `TODO.md` item *Keep environmental effects stable when spatial or latent-factor terms are added*.
- **A spatial residual is evidence of a pattern, not of demography.** A positive field in poor habitat is what a sink would produce, and also what a missing, spatially smooth covariate would produce. Traps that intercept flying insects also catch individuals in transit. The model can show where occurrence and numbers exceed what habitat predicts; calling that a sink needs an argument from outside the model.
- **Many rare species.** With several hundred species, many seen a handful of times, a separate size, `tau_0` and `tau` per species will be poorly determined without pooling. See questions 5 and 17.
- **Sequence read counts must not be fitted with this model.** They are not counts of organisms: they are compositional within a sequencing run and carry species-specific amplification bias, so their size says nothing absolute about abundance. The hazard that remains is practical: a table of read counts is also a table of whole numbers, so nothing in the data stops a user from fitting one. The help page should say so, and the model should never be chosen automatically from the data's type (question 10).

### A check from species-richness estimators

Doug asked for this on 16 November 2026. Estimators of sample completeness, such as Chao1 and the sample coverage reported by the iNEXT package, estimate from a single sample how many species at a plot were missed, using the numbers of species caught once and twice. They are a check on the model and on the data. They are not an input to the model and not a substitute for replicate samples.

- **A data check before the first fit.** Sample coverage per plot shows which plots were sampled thinly. Where trap-days differ between plots, its relation to trap-days shows whether the effort offset has something to explain.
- **A model check after fitting.** The shared-predictor model gives, for each plot, the expected number of species present but not caught: the sum, over the species with a zero count there, of `P(z = 1 | y = 0)`. Compare it plot by plot with the Chao1 estimate of undetected species. A large disagreement points at the count distribution or at the occupancy part.
- **Simulation first.** How closely the two should agree is not known. Stage 3 makes the comparison on simulated data, where the true number of missed species is known, before it is read on real data.
- **It names no species.** The estimators give a total per plot. The model needs the chance of a miss for each species at each plot, and the total cannot supply it.
- **It is not independent evidence.** The estimators and the count model both work from how many individuals of each species were caught in the one sample. Agreement shows that the two are consistent, not that detection has been measured. Replicate samples would measure it without assuming a shape for the counts.
- **The two totals differ in scope.** Chao1 is a lower bound and includes species never caught anywhere in the survey, which the model cannot see. Sample coverage is weighted by individuals, so a plot can have high coverage and still lack many of its rare species.
- **Longer trapping is already in the model.** More trap-days raise the expected count through the offset, so a zero from a long exposure is stronger evidence of absence. Nothing extra is needed for that.

## Proposed interface

Open to change; see question 10.

- A new argument to `runOccJSDM()`, `countModel`, `NULL` by default. When `NULL`, behaviour is as now, except that the "not supported" error tells the user how to ask for a count model.
- When set, a list: `family` ("negbin" or "poisson"), `occupancy` (`TRUE` for the shared-predictor model, `FALSE` for the plain count JSDM) and `offset` (the name of an effort column in `data$info`, on the natural scale).
- With `countModel` set, repeated rows per site are replicate counts and `threshold` is not applied.
- New stored output: the draws of the size, `tau_0` and `tau`, with their convergence diagnostics; the expected count where present. The occupancy state and probability are stored as the occupancy model stores them now.
- `simulateOccJSDMData()` gains both count models, and its `rnb` input, documented today as having no effect, becomes live.

### Example calls

A presence and absence or eDNA fit is written exactly as it is today:

```r
fit <- runOccJSDM(data, occCovariates = c("elevation", "canopy"))
```

The shared-predictor count model:

```r
fit <- runOccJSDM(
  data,
  occCovariates = c("elevation", "canopy"),
  countModel = list(family = "negbin", occupancy = TRUE, offset = "trap_days")
)
```

The plain count JSDM, with no occupancy layer:

```r
fit <- runOccJSDM(
  data,
  occCovariates = c("elevation", "canopy"),
  countModel = list(family = "negbin", occupancy = FALSE, offset = "trap_days")
)
```

### What the covariate arguments mean

- In the shared-predictor model `occCovariates` act on occupancy, as now, and through `tau` on abundance. No new covariate argument is needed.
- In the plain count JSDM there is no occupancy, so `occCovariates` act on the log of the expected count, despite the name. The help page says so.
- `collCovariates` are not used by either count model in the first version: differences between samples enter through the offset only.

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
- **The shared-predictor model leaves embeddings where they are**: occupancy covariates on the logit scale, under a prior calibrated for that scale. The combination is supported from the start and tested in stage 3. One thing to watch: the embedding's contribution now also moves the expected count, multiplied by `tau`, so a prior draw that is large on the logit scale is larger still in its effect on counts.
- **The plain count JSDM puts embeddings on the log scale**, where the grouped prior's default has not been checked. The first use does not need that combination, so `countModel` with `occupancy = FALSE` stops with a clear error when `embeddingVars` is set. Lifting it is left for later; see question 15.
- **Taxonomy depends on that branch.** Count fits get taxonomic borrowing only from part 2 of the grouped-prior work (see *Decided*). Until part 2 exists, rare species in a count fit borrow strength through the trait regression on their coefficients, through the community prior on the occupancy intercepts, and through whatever pooling is decided for the size, `tau_0` and `tau`. When part 2 is built, it has to be tested on a count fit as well as on the binary and two-stage fits.
- **A guard that already exists.** The embedding branch ships a reference fit that must be reproduced when its prior is off. The count work must keep that test passing.

## Questions for Alex and Doug

Questions keep their numbers from earlier versions of this document. Those marked resolved are kept for the record.

1. **(Doug) Which data come first? Resolved 12 November 2026:** counts of individuals. Sequence read counts are excluded; see *Decided*.
2. **Where should the JSDM structure sit? Resolved 15 November 2026:** one JSDM with a shared predictor, configuration D, proposed by Alex and adopted by Doug. See *Configurations considered*.
3. **Polya-Gamma negative binomial or Poisson-lognormal?** The shared predictor relies on the Polya-Gamma form, so this is settled in its favour unless Alex wants the Poisson-lognormal route kept open for the eDITH fork.
4. **Is Poisson as a large fixed size acceptable**, or should Poisson be left out of the first version?
5. **The size.** One per species, as the simulator has it, or a shared or partially pooled size? Which prior: a gamma on the size, or a half-normal on its inverse square root? The first use has several hundred species, many of them rare, which favours partial pooling.
6. **Counts at unoccupied sites.** Exactly zero, as written above, so that any positive count fixes the site as occupied? Or a small false-positive count rate per species, the analogue of `theta0`? With counts of individuals a positive count at an unoccupied site can only be a misidentification, so the recommendation is exactly zero in the first version.
7. **Occupancy with one sample per site.** Now part of question 18.
8. **The intercepts.** In the shared-predictor model the JSDM intercepts are occupancy intercepts and keep the community prior they have by default. For the plain count JSDM, where the intercepts are log expected counts: fixed prior or community prior?
9. **The Polya-Gamma approximation.** Is a moment-matched normal draw above a shape cut-off acceptable in a sampler otherwise exact, and is 170 the right cut-off?
10. **(Doug and Alex) Interface.** One `countModel` list, as proposed, or separate arguments? Should `occupancy` default to `TRUE`? Should integer data without site replication fit a count model automatically, as the current inference half-intends? The recommendation on the last is no: a table of sequence reads is also whole numbers, and automatic inference would fit it as counts of individuals.
11. **Site-level WAIC.** Leave counts out of `computeSiteWAIC()` at first, with a clear error, and keep only the complete-data WAIC?
12. **Variation partitioning.** In the shared-predictor model, partition on the logit scale of the predictor, as now? In the plain count JSDM the stub takes standard deviations of `exp(eta)`, which the sites with the largest means will dominate; partition on the log scale instead?
13. **(Doug) Where and when.** Pull requests to `main` from this branch, one per stage, or in Qin's fork? Before or after the beta work in `TODO.md`?
14. **Covariate argument names. Resolved 15 November 2026** by the shared predictor: `occCovariates` keeps its meaning. See *What the covariate arguments mean*.
15. **Embeddings on the log scale.** Needed only for the plain count JSDM with an embedding block, which the first use does not need. Deferred.
16. **What is in `X_i beta`? Resolved 23 November 2026 (Alex):** the whole JSDM predictor, including the latent factors and the spatial field.
17. **Are `tau_0` and `tau` per species? Answered 23 November 2026 (Alex):** both per species to begin with, to be changed later if they prove unidentifiable. So there is no pooling of `tau` across species in the first version. The priors are not yet chosen; Claude will propose them in the stage 3 plan.
18. **How is the scale of `eta` kept identified** when there is one sample per site and most of the information is in the counts? **Answered 23 and 25 November 2026 (Alex):** by an interweaving step that updates `tau` and the predictor jointly; see *The interweaving step*. Alex confirmed Claude's working of step (c) on 25 November 2026, and that it needs a Metropolis update. Whether the step is enough at one sample per site is measured in stage 3.
19. **The collapsed intercept shift in count fits. Answered 25 November 2026 (Alex):** try option 3, shifting `tau_0` to compensate, and possibly option 1 with it. He expects the slow mixing may not appear in count fits. Stage 3 builds both moves and compares them with no shift; see *The collapsed intercept shift in count fits*.

## Staged plan

Each stage is one pull request to Alex and leaves `main` working. Stage 2 starts from `main` after the embedding branch has merged; see *Interaction with the embedding covariates*. A task-by-task plan for each stage is written in `dev/superpowers/plans/` once its questions are answered: stage 1 needs questions 5 and 9, stage 2 needs 4, 8, 10, 11 and 12, and stage 3 needs 6.

- **Stage 1: likelihood pieces, tested alone.** The general-shape Polya-Gamma sampler, checked against the exact mean and variance across shapes and against the existing sampler at whole-number shapes. A negative-binomial log-likelihood with offsets. The size update with its prior, checked on simulated counts with a known predictor. One name for the model, "counts", in R and C++. No change to any fit.
- **Stage 2: the plain count JSDM.** Wire the count branch into `update_jSDMcoef()`; accept `countModel`; offsets and replicate samples; size draws in the output, the posterior-draw helpers and the diagnostics; a log link in the covariate-response plots; the complete-data WAIC with the negative-binomial density; `simulateOccJSDMData()` with a count model; a posterior predictive check that compares observed and predicted numbers of zeros per species. Done when a simulate-and-recover test passes for coefficients, size and residual correlations, existing fits of the other four models are unchanged draw for draw, and a comparison with an independent negative-binomial JSDM (Hmsc or gllvm) on one simulated dataset agrees within Monte Carlo error.
- **Stage 3: the shared-predictor model.** The combined `k` and `Omega`; the update of `tau_0` and `tau`, with the truncation and the priors proposed in this stage's plan (question 17); the interweaving step (question 18); the count version of the occupancy draw; the two intercept shift moves and their comparison with no shift (question 19); draws of `tau_0` and `tau` in the output and the diagnostics; the predictions listed under *What the shared predictor gives*; the simulator's shared-predictor case; a fit with an embedding block. Done when recovery and interval coverage of the coefficients, `tau_0`, `tau` and the size are reported at 1, 3 and 5 samples per site, with the effective sample size of `tau` with and without the interweaving step, so that question 18 is answered by measurement, when the model's expected number of missed species per plot has been compared with the truth and with Chao1 on the same simulated data (see *A check from species-richness estimators*), and when existing fits are still unchanged draw for draw.
- **Stage 4: is the first use's spatial aim achievable?** A simulation, not a feature, with three arms on the same simulated sites: presence and absence, the plain count JSDM and the shared-predictor model. Part one: on communities from the spatial recoverability pilot's design, compare recovered amplitude, range and held-out field. Part two, the `TODO.md` item *Source-sink inference scenario*: simulate species whose environmental and spatial effects oppose each other, at the first use's number of plots, with and without an embedding block, and ask whether the fitted field is positive where the environmental prediction is low. Part three: simulate species whose abundance and occupancy respond differently to a covariate, and report what the shared-predictor model does to the coefficients, the field and the predictions. Run detached, at most five workers, power mode recorded.
- **Stage 5: documentation.** The count help page, `?runOccJSDM`, `NEWS.md`, the README line "Count responses are not supported", and a lesson section once Doug has a dataset to teach from.

## Non-goals

- No model of sequence read counts, with or without spike-ins; that remains the eDNAPlus item, and it estimates changes in abundance, not absolute abundance.
- No configuration A, B or C.
- No sample-level covariates on the counts beyond the offset.
- No embedding block in the plain count JSDM.
- No taxonomic or phylogenetic structure; that belongs to part 2 of the grouped-prior work.
- No change to the binary, continuous, occupancy or two-stage models, their defaults or their stored output.

## References

- Chao (1984), Nonparametric estimation of the number of classes in a population, *Scandinavian Journal of Statistics* 11.
- Chao and Jost (2012), Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size, *Ecology* 93.
- Hsieh, Ma and Chao (2016), iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers), *Methods in Ecology and Evolution* 7.
- Pillow and Scott (2012), Fully Bayesian inference for neural models with negative-binomial spiking, *Advances in Neural Information Processing Systems* 25.
- Polson, Scott and Windle (2013), Bayesian inference for logistic models using Polya-Gamma latent variables, *Journal of the American Statistical Association* 108.
- Windle, Polson and Scott (2014), Sampling Polya-Gamma random variates: alternate and approximate techniques, arXiv 1405.0506.
- Zhou, Li, Dunson and Carin (2012), Lognormal and gamma mixed negative binomial regression, *Proceedings of the 29th International Conference on Machine Learning*.
- Diana and others (2024), eDNAPlus, *Journal of the American Statistical Association*, as cited in `TODO.md`.

These citations are from memory and have not been checked against Zotero.
