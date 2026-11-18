# Grouped covariate prior for occupancy coefficients

Date: 11 November 2026. Status: design written from a conversation with Doug on 11 November 2026. As of 12 November 2026 part 1's package code is implemented and unit-tested on branch `claude/grouped-covariate-prior`; the real embedding extract, the validation study and part 2 (taxonomy) are not started. A prior check on 12 November 2026 found the first default `axis_scale` far too wide; Doug chose `0.7 / K` the same day and it is implemented. Two weaknesses seen under the new default, and a numerical warning, are recorded under *What the prior check found*. A whole-branch review on 12 November 2026 found the sampler correct and fits under the fixed prior unchanged; its fixes, among them the starting state of the blend weights, are in. Amended on 12 November 2026 with the section *What the management coefficient estimates*. This is the `TODO.md` item "Regularised covariate selection" under *Other model extensions*. It is project A of five; projects B (running with thousands of species), C (prediction to planning units), D (a baseline analysis of the first real dataset) and E (occupancy from the beetle counts) are described at the end and are not designed here.

## Purpose

Let `runOccJSDM()` take a geospatial embedding (AlphaEarth, Tessera, OlmoEarth or similar: tens to a hundred or so uninterpretable numeric columns per site) as occupancy covariates alongside a few measured covariates, and let species borrow strength from their taxonomic relatives, in addition to the existing borrowing from species with similar traits, which is kept. Taxonomy matters most where there are no traits. Today every occupancy coefficient shares one prior scale that the data cannot move. This design replaces that single scale, as an option, with scales per group of covariates, per embedding axis and per taxonomic rank, all learned from the data.

## Why

Doug wants to predict continuous species distributions and run systematic conservation prioritisation on the predicted maps. The embedding goes in together with measured covariates such as elevation and lidar metrics, and with management covariates such as protected-area boundaries or time since clearcut. For a management covariate he wants one conclusion: whether it has an effect beyond what the satellite sees. He wants this in occJSDM, not in sjSDM or Hmsc, for three reasons: full posterior draws of the whole community to carry uncertainty into the prioritisation, the trait regression on coefficients, and one model for unreplicated and replicated datasets.

The first datasets are BioDivAbove (Norway, west and north of Lake Hurdal): soil fungi, beetles, deadwood fungi and malaise-trap arthropods on one shared set of about 220 to 270 plots, with no sample or PCR replicates, so they are fitted with `model = "binary"`. Only the beetles have traits; the other three have taxonomy. Later the same feature will run on the Gaoligongshan data, which is two-stage. So the feature must work in both modes.

The embedding therefore has two jobs. For the map it is a predictor, where strong shrinkage is safe. For the management question it is a control that makes sites comparable in habitat, where too much shrinkage leaves habitat differences unexplained and the management coefficient absorbs them. The prior has to serve both, and the validation has to test both.

## What the management coefficient estimates

Added on 12 November 2026, after Doug asked whether "an effect beyond what the satellite sees" is a kind of covariate matching for causal analysis. It is, with one complication. This section fixes what the coefficient means and what has to be reported with it. It changes no code: every fit it calls for can be run with the interface below.

**It is covariate adjustment.** Matching compares managed and unmanaged sites that are alike in their covariates. The model does the same by regression: it estimates the management effect with the other covariates held fixed. Both rest on the same two assumptions.

- *Nothing important is missing.* Everything that influences both where management happened and where species occur is among the covariates. Causes with a spatial pattern, such as remoteness, count; a spatial fit is one way to check for them.
- *Overlap.* Managed sites have unmanaged counterparts that look alike, and the reverse. Matching shows a failure of overlap as unmatched sites. A regression extrapolates and gives no warning, so overlap has to be checked separately.

The model adds a third assumption that matching does not need: the embedding acts linearly on the logit scale.

**The embedding is measured after the management.** In ordinary matching the covariates describe a site before it was treated. The embedding describes the vegetation now, and the vegetation now is partly a result of the management: a logged stand looks logged from space. So the embedding is partly a pathway through which management acts, and holding it fixed removes that pathway. That is what "beyond what the satellite sees" asks for. The management coefficient from a fit with the embedding is therefore a *direct effect*: the effect of management that does not pass through habitat the satellite can see. It is not the *total effect* of management.

A direct effect needs one further assumption: nothing unmeasured affects both the vegetation and the species. An example of how it fails: an old stand and a once-logged stand look the same from space. For them to look the same despite their histories, they probably differ in something else, such as soil fertility. If fertility also affects the species, the comparison is biased, and it is including the embedding that caused the bias. This cannot be tested from the data and has to be argued for each study.

**Covariates are classified twice.** The group of a covariate (measured or embedding) sets its prior. Its timing (fixed before the management, or able to change after it) sets what the management coefficient means. The two do not line up. Lidar terrain metrics are measured covariates fixed before management. Lidar canopy metrics and field measurements of stand structure or deadwood are measured covariates that management can change, so for this purpose they behave like the embedding.

**What an analysis with a management covariate reports.**

- *A total-effect fit:* the management covariate with only the covariates fixed before management. No embedding and no other covariate that management can change.
- *A direct-effect fit:* the same plus the embedding. Whether the other covariates that management can change also go in is decided for each study, and it changes the question to "beyond what the satellite and the lidar see".
- *Both, side by side, on the probability scale:* the change in predicted occupancy probability between levels of management, averaged over the sites. The two logit coefficients are not compared directly, because a logit coefficient grows when any covariate that predicts occupancy is added, even one unrelated to management. The gap between the two fits is then read, on the assumptions above, as the part of the effect that acts through visible habitat.
- *An overlap check:* how much of the variation in the management covariate is left after regressing it on the embedding and the other covariates. This works for a continuous covariate such as time since clearcut. If little is left, the direct effect is reported as not estimable, whatever the fit returns.
- *A matched check, for a two-level contrast:* match managed and unmanaged plots on the leading principal components of the embedding and the covariates fixed before management, refit on the matched plots, and see whether the answer holds.
- *An embedding from before the management, where one exists:* AlphaEarth starts in 2017. For a management change since then, the embedding of an earlier year describes the site before it was treated and can go in the total-effect fit as an ordinary control. This is not available for logging history back to 1936.
- *A statement about catchability:* a `binary` fit cannot separate an effect of management on occurrence from an effect on how well a trap catches, for example a window trap in an open clearcut. In `occupancy` and `two_stage` fits the management covariate also goes in the detection stage.

These are steps in the analysis of real data, first in project D. Within project A, the documentation of `embeddingVars` states that a coefficient on a management covariate fitted alongside the embedding is a direct effect and points to the total-effect fit, and the validation notes what truth 4 does and does not test.

## What the code does today

Verified on 11 November 2026 against `main` at `4f2ee6c`. Line numbers refer to that commit.

- The occupancy coefficients `B` are a covariates by species matrix. Each entry has a Normal prior whose mean is the species' row of `Tr %*% G + A %*% C` and whose standard deviation is the single number `sigma_b`. The joint draw of the intercept, `B`, the factor loadings and the spatial coefficients fills its prior diagonal from that scalar (`src/jsdm.cpp` lines 880 to 890 and 1100 to 1112).
- The variance of `sigma_b` has an InvGamma(10, 1) prior with no option (`R/runOccJSDM.R` line 1062), updated by `sample_sigmab()` (`R/jsdmfun.R` line 876). The simulation studies of 15 to 17 October 2026 found it sits at its prior mean of about 0.33 in every fit.
- `Tr` is the species by traits matrix. It has no intercept column (`remove_intercept = TRUE`, `R/runOccJSDM.R` lines 952 to 956) and has zero columns when no traits are supplied (line 962). So there is no explicit community-average response to a covariate. Categorical traits are accepted and dummy coded with a reference level.
- `A` is a species by `gt` matrix of latent traits and `C` is their `gt` by covariates weights. `gt` defaults to `floor(sqrt(min(S, ncov_psi)))` (line 992) and is set by `listParams$n_lattrait`.
- `sample_GC()` (`R/jsdmfun.R` line 1016) updates one covariate at a time: it regresses that covariate's row of `B` on `cbind(Tr, A)` with noise `sigma_b` and an independent Normal(0, 1) prior on every entry. `sampleBuniv()` (`src/jsdm.cpp` line 517) does not scale that prior by the noise.
- `sample_A()` (line 1064) updates one species at a time: it regresses that species' coefficients, less the trait part, on `t(C)` with noise `sigma_b` and a Normal(0, 1) prior.
- The spatial range step integrates the same coefficient block, and builds its prior precision from the scalar `sigma_b` in `spatial_range_logweights()` (line 1131) and its grouped variant (line 1192).
- The same `Tr` also drives the trait regression of the spatial coefficients (`Gs`, `As`, `Cs`, lines 1413 to 1426).
- Within an iteration the order is: range, joint coefficient draw, `mu_b0`, `sigma_b`, `sample_GC()`, `sample_A()` (lines 1347 to 1411).
- A half-Cauchy update for a shared coefficient scale already exists as `sample_spatial_sd_half_cauchy()` (line 897), used for `sigma_bs`. The tag `archive/sigma-b-prior` applied it to `sigma_b` under `sigma_b_prior = "half_cauchy"`; that code is not on `main`.
- `nthin` thins before storage: the loop runs `nburn + niter * nthin` iterations and keeps `niter` (`R/chains.R` line 645). The output arrays are allocated in full (`R/runOccJSDM.R` lines 1131 to 1143).
- The internal simulator draws the coefficient leftovers with one `sigma_b` and builds `C` upper triangular with a unit diagonal (`R/jsdmfun.R` lines 669 to 675).

## Goals

- An opt-in prior in which measured covariates and embedding axes are shrunk differently, with every scale learned.
- The same for taxonomy: a species' response is centred on effects of its order, family and genus (or whatever ranks are supplied), with a learned scale per rank.
- An explicit community-average response to each covariate.
- Works for `binary`, `occupancy` and `two_stage` fits, spatial or not, serial or parallel chains.
- With the option off, every existing fit is reproduced bit for bit.
- A simulation study that shows, before anyone relies on it, that the prior passes the three criteria under *Validation*.

## Non-goals

- No scale per individual coefficient (no species-by-axis horseshoe) and no inclusion probabilities.
- No change to the spatial coefficients, the factor model, the detection stages or the intercept prior.
- No change to any default for a call that does not name `embeddingVars`: there `b_prior` stays `"fixed"` until Doug decides otherwise after validation. Naming `embeddingVars` switches the grouped prior on, as set out under *Interface*.
- No nonlinear terms on embedding axes. Naming an embedding column in `splineVars` is an error.
- No exported package function to download or summarise embeddings. The validation study downloads a real extract with study scripts under `dev/simstudy/`, described under *Real embedding extract*.
- No work on memory or speed at thousands of species (project B) and no change to `predictNewSites()` (project C).

## Anti-goals

- A prior that improves held-out prediction while biasing the management coefficient is a failure, even if every test passes otherwise.
- Any change to results when the option is off.
- Presenting an individual embedding axis or blend as an ecological finding. The blends are identified only up to rotation.
- Presenting the management coefficient from a fit with the embedding as the total effect of management. It is the direct effect, as set out under *What the management coefficient estimates*.

## The modelling strategy in brief, for Alex

Added on 12 November 2026, from a description Doug asked for. It restates the two model sections below in a statistician's terms and ends with the points on which Alex's view is wanted.

**Aim.** Fit K = 64 to 128 embedding columns as occupancy covariates beside a few measured ones, with about 250 sites and 100 or more species. The prior is opt-in, and fits under the existing prior are reproduced bit for bit.

**Model.** Only the prior on `B` changes. For covariate k and species s:

```
B[k,s]  = (Tr_b G)[s,k] + (A C)[s,k] + e[k,s]
Tr_b    = cbind(1, Tr)             row 1 of G is a community-average response
e[k,s]  ~ N(0, sd_k^2)             sd_k = sigma_b on a measured covariate,
                                   sigma_emb * lambda_k on an embedding column
G[,k]   ~ N(0, lambda_k^2)         lambda_k = 1 on measured covariates
C[,k]   ~ N(0, lambda_k^2 / gt)    on embedding columns; N(0, 1) on measured ones
A[s,]   ~ N(0, I)
sigma_b ~ C+(1),  sigma_emb ~ C+(1),  lambda_k ~ C+(0.7 / K)
```

So `lambda_k` rescales the whole of row k, centre and leftover alike. It is a horseshoe-type scale per column, not per coefficient, and there are no inclusion indicators. `A C` over the embedding columns is a reduced-rank regression of rank `n_lattrait`, and the latent traits are shared with the measured covariates.

**Sampler.** Everything stays conjugate. The three kinds of scale use the inverse-gamma auxiliary form of the half-Cauchy. The joint draw of the coefficient block and the spatial range step take a prior standard deviation per covariate in place of the scalar. `sample_GC()` uses noise `sd_k` per covariate, and `sample_A()` becomes a weighted regression. The scales are drawn after the joint coefficient draw, where `sample_sigmab()` ran.

**Known problems.** The first four are set out under *What the prior check found*.

- *Prior scale.* A scale of 1 over the square root of K on `lambda_k` gave the embedding a prior contribution with a median of 16 logits at K = 64, because the heavy tail makes it grow as K times the scale. At `0.7 / K` the median is about 1.5, but the 95% point is about 28.
- *Mixing.* In a small test the `lambda_k` of inactive columns had not converged across chains at 1,500 burn-in and 1,500 kept iterations.
- *Starting values.* `C` started at ones, which threw every `lambda_k` far above its prior scale in the first sweep, with a transient longer than 2,000 sweeps that R-hat cannot see because all chains share it. The embedding columns of `C` now start at their prior scale. Whether a tight or an overdispersed start is better is open.
- *Numerics.* `chol` in `sampleBuniv()` warned about symmetry when `sd_k` was very small. Its argument is now made exactly symmetric, which changes no draw.
- *Replicated data.* In `occupancy` and `two_stage` fits the community-average row and the scales reach the data only through the latent occupancy states, which is the structure that made `mu_b0` fail to converge before the collapsed shift. Short fits run, but convergence there has not been examined.

**Where Alex's view would help.**

- Half-Cauchy, half-t or a regularised horseshoe on `lambda_k`, given binary data with many rare species?
- Is rescaling the whole row, `G` and `C` included, better than scaling the leftovers only?
- Is there a blocked or collapsed move for `lambda_k` together with row k?
- Is a collapsed move like the one for `mu_b0` likely to be needed for the community-average row in replicated fits?

## The model, part 1: covariate side

Every occupancy covariate column belongs to one of two groups: *measured* (the default) or *embedding* (the columns the user names). Factor and spline columns are always measured.

A species' coefficient on a covariate is a centre plus a leftover.

- The centre is, as now, the species' row of the trait matrix times `G`, plus its latent traits times `C`. The trait matrix used on the covariate side gains a leading column of ones, so the first row of `G` is the community-average response to each covariate. The spatial side keeps the old matrix.
- The leftover is Normal with mean 0 and standard deviation equal to the group scale times the axis scale of that covariate.

The scales are:

- **Measured-group scale**, one number. Half-Cauchy prior with scale `sigma_b_scale`, default 1.
- **Embedding leftover scale**, one number. Half-Cauchy prior with scale `sigma_emb_scale`, default 1.
- **Axis scales**, one per embedding column, fixed at 1 for measured columns. Half-Cauchy prior with scale `axis_scale`, default `0.7 / K` for K embedding columns. The first default, 1 over the square root of K, was changed on 12 November 2026; see *What the prior check found* below.

The axis scale of a covariate also multiplies the prior standard deviation of that covariate's column of `G` and of `C`. So an axis scale multiplies everything in that covariate's row of `B`, centre and leftover alike; it is equivalent to the model rescaling that embedding column. An axis no species responds to has its scale shrink towards 0 and drops out of every blend.

The defaults are meant to keep the embedding's prior contribution to the linear predictor from growing with the number of embedding columns or the number of blends. The size of the embedding is carried by the axis scales, which act on the whole row. On embedding columns the base prior standard deviation of the entries of `C` is 1 over the square root of `gt`, where today it is 1, so the blends together contribute about as much prior variance as the community average or the leftover. The target, checked by prior simulation, is a prior median of 1 to 2 for the standard deviation across sites of the embedding's contribution to a species' linear predictor.

**What the prior check found.** Added on 12 November 2026. The first default for `axis_scale` was 1 over the square root of K, for K embedding columns. Claude set it by treating each axis scale as if it sat at its prior scale, which gives the whole embedding a contribution of about 1.7 logits. That reasoning ignored the tail of the half-Cauchy: among K draws a few are very large and they dominate the total, so the contribution grows in proportion to K times the scale. The prior check in `dev/simstudy/grouped-prior-checks-20261111/`, with 8 latent traits, gave a prior median of about 16 logits at 64 columns and 24 at 128. Three independent simulations agreed. One of them, a reviewer's scratch run that is not kept in the repository, also ran 512 columns and confirmed that the scale must fall as 1 over K.

- *The default now:* `0.7 / K`, chosen by Doug on 12 November 2026 on Claude's recommendation, keeping the half-Cauchy. The prior median is 1.4 at 64 columns and 1.5 at 128.
- *Known weakness, the upper tail:* the 95% point is about 28 logits at both column counts, so one prior draw in twenty gives occupancy probabilities that are almost all 0 or 1. Part of the tail does not come from the axis scales: with them held fixed the 95% point is still about 12, because the embedding leftover scale is also half-Cauchy. The alternative Doug did not choose, a half-t with 3 degrees of freedom on the axis scales at 0.6 over the square root of K, gives a median of 1.6 to 1.7 and a 95% point of 11 to 14, at the cost of weaker selection among columns. It is the fallback if the validation study shows the tail does harm, alongside the regularised horseshoe named below.
- *Known weakness, mixing of the axis scales:* in the small shrinkage test (150 sites, 30 species, 10 embedding columns of which 3 are active, 2 chains of 1,500 burn-in and 1,500 kept) the two strongly active columns are unmistakable: their pooled median axis scales are 33 and 9.1 times the largest inactive one. The weakly active column is 2.4 times it. The axis scales of the inactive columns had not converged at that chain length: R-hat reached 1.49 at one seed and 4 to 7 at three others. The validation's convergence rule covers the axis scales, and this result warns that its 5,000 burn-in and 5,000 kept iterations may not be enough for them.
- *The starting state, found by the whole-branch review:* as first built, the blend weights `C` started as a block of ones, so the first sweep threw every axis scale far above its prior scale, to about `gt` over the square root of the number of species plus trait rows plus `gt`, whatever `axis_scale` was. The transient outlasted 2,000 sweeps, and because every chain shared it, R-hat could not see it. Under that start the weakly active column above was separated from the inactive ones by a factor of only 1.10. Claude ruled that under the grouped prior the embedding columns of `C` start at their prior scale, so that the first sweep agrees with the stated starting axis scale. There is a trade-off for the validation study to settle, with chains started at different axis scales: in the reviewer's experiment with 32 embedding columns, 3 of them active, the inflated start overestimated the embedding's contribution (3.74 against a truth of 3.02, and 3.13 with a tight start) but found a weakly active column that the tight start missed.
- *A numerical warning, now silenced:* one fit printed `chol(): given matrix is not symmetric` three times in one burn-in sweep, from the draw of `G` and `C`, for the three embedding columns whose axis scales were smallest at that moment. The cause is a rounding-level asymmetry, a relative mismatch of 2.4e-11, in a matrix that is symmetric by construction, magnified by division by a very small variance. Armadillo's check trips at a tolerance of 1e4 times machine epsilon and its factorisation reads only the upper triangle, so the draw was not affected. The argument of the factorisation is now made exactly symmetric: rerunning that fit gives no warning and identical draws, and fits under the fixed prior are unchanged.

The "blends" Doug and Claude discussed are the existing latent traits restricted to the embedding columns: a blend's weights are a row of `C` over those columns, and a species' response to the blend is its entry in `A`. No new structure is added. The number of blends is the existing `n_lattrait`, with its existing default. The latent traits stay shared between the two groups, so a species' latent traits predict its responses to measured covariates and embedding axes alike.

**What a blend is, in plain terms.** Added on 12 November 2026 from Doug's questions. Two figures go with it: [a PCA axis against a blend on the same sites](2026-11-11-grouped-covariate-prior-figures/blend-versus-pca.svg), and [species' vectors lying in the sheet that two blends span](2026-11-11-grouped-covariate-prior-figures/species-vectors-in-blend-sheet.svg). Both are schematic drawings with made-up points, not results.

- *A blend is a direction in the embedding space.* It is a set of weights, one per embedding column: 64 weights for AlphaEarth.
- *Each species has its own vector.* Its 64 coefficients point in the direction along which its chance of occurring changes fastest, which for one species is the direction that best separates its presences from its absences. With 100 species there are 100 such vectors.
- *The model does not find those vectors first and then compress them.* It assumes from the start that every species' vector is built from the same small set of shared vectors, the blends: a species' vector is the community average, plus its weight on blend 1 times blend 1, plus its weight on blend 2 times blend 2, and so on, plus a small leftover of its own. The blends are common to all species, the weights differ between species, and both are estimated together in one fit.
- *The geometric picture.* Draw the species' vectors as arrows from one origin in the 64-dimensional space. The assumption is that they do not point in unrelated directions but lie close to a flat sheet of low dimension: a plane for two blends, an eight-dimensional sheet for eight. The blends are the axes of that sheet, and a species' weights are its coordinates on it. The leftover is how far a species sits off the sheet, and because the leftover scale is learned, species are held close to the sheet when the blends explain the responses well.
- *How this differs from principal components of the embedding.* PCA sees only the embedding table of sites by columns, and picks the combinations of columns along which the sites differ most, whether or not any species responds to them. A blend is picked because species respond to it, however little of the embedding's own variation it carries. Suppose sites differ a great deal along column 1 and only a little along column 2, and species respond only to column 2: one principal component keeps column 1 and is useless, while one blend puts its weight on column 2. In ordination terms, regressing species on principal components of the embedding is an unconstrained ordination of the environment, and the blends are like RDA or CCA axes. The blends are also close to what a PCA across species of their separately fitted coefficient vectors would give, which is a different PCA from one across sites.
- *Why one fit and not two steps.* A species with five detections cannot determine 64 numbers, but it can determine eight weights on blends that the other species have located. Separate vectors for 100 species need 6,400 numbers; eight blends need 8 times 64 for the blends and 100 times 8 for the weights, about 1,300.
- *How many blends.* As many as there are latent traits, set by `n_lattrait`. Two latent traits give exactly two blends. The reporting function can return no more blends than there are embedding columns or species.
- *What the axis scales do to the picture.* An embedding column that no species responds to has its axis scale shrink towards zero, which flattens the sheet along that column.
- *The blends are not unique.* Any rotation of the axes spans the same sheet and gives the same fit, so the sheet is meaningful and an individual blend is not. This is the reason for the anti-goal about presenting a blend as an ecological finding.
- *The cost.* Principal components are fixed before the fit and cost nothing. Blend weights are estimated from the same data as everything else, so with few sites or weak responses they are noisy. That is why the validation has a *components* arm: it tests whether learning the blends beats fixing them in advance at about 250 sites.

All updates are conjugate given the others, using the auxiliary-variable form of the half-Cauchy that `sample_spatial_sd_half_cauchy()` already uses:

- Each group scale: inverse-gamma, from the sum over that group of squared leftovers each divided by its squared axis scale.
- Each axis scale: inverse-gamma, from that covariate's squared leftovers divided by the squared group scale, plus the squared entries of its columns of `G` and `C` divided by their base prior variances.
- `sample_GC()`: noise and prior standard deviations become specific to the covariate.
- `sample_A()`: the regression becomes weighted, because the covariates now have different noise.
- The joint coefficient draw and both range functions take a vector of prior standard deviations, one per covariate, in place of the scalar. The vector is the same for all species, so the shared prior precision matrix stays shared.

If the simulations show axis scales running away (possible with a half-Cauchy and weakly informative binary data), the fallback is the regularised horseshoe, which caps large scales with a slab. The builder reports to Doug before adopting it.

## The model, part 2: taxonomy

The user supplies a table with one row per species and one column per rank, broadest first. A species may be unassigned below some rank. The covariate-side trait matrix gains one indicator column per taxon at every rank, with no reference level dropped. Columns that contain a single species, or exactly the same species as a column of a broader rank, are removed, since they cannot be told apart from the species' leftover or from the broader taxon.

The effect of a taxon on a covariate is Normal with mean 0 and standard deviation equal to a rank scale times the covariate's axis scale. There is one rank scale per rank and per covariate group, each with a half-Cauchy prior of scale `taxon_scale`, default 1.

A species' centre is then: community average, plus its order's effect, plus its family's, plus its genus's, plus the parts from any real traits and from the latent traits. A species unassigned below family takes the order and family effects only.

Taxa within a rank do not overlap, so each rank is updated as a block of independent Normal draws, one per taxon and covariate, from the mean of the residuals of the species in that taxon. This keeps the cost linear in the number of species; a joint draw over hundreds of taxon columns would cost a large Cholesky factorisation per covariate per iteration. The rank scales have inverse-gamma updates of the same form as the group scales.

## Interface

- There is one fitting function. Doug asked on 12 November 2026 whether the feature should have a second function, so that it cannot complicate standard use. It stays in `runOccJSDM()`, with no new argument: a second function would run the same sampler, would have to be tested against every accessor, and would become redundant if the prior is later made a default. Standard use is protected by the default, by the reproduction tests, and by the new outputs existing only under the grouped prior.
- `listParams$embeddingVars`: character vector of names from `occCovariates`. Each must be numeric and absent from `splineVars`. Naming it is enough to turn the grouped prior on (Doug's decision, 12 November 2026). This also changes the prior on the measured covariates, to the learned measured-group scale, and adds the community-average response; the help page says so and the fitted object records which prior was used.
- `listPriors$b_prior`: `"fixed"` (today's prior) or `"grouped"`. When it is absent it is `"grouped"` if `embeddingVars` is named and `"fixed"` otherwise. Setting it to `"grouped"` is needed only to use the prior without an embedding. `"fixed"` together with `embeddingVars` is an error.
- Embedding columns listed in `occCovariates` but not in `embeddingVars` are ordinary measured covariates under whichever prior applies. This is how the *raw* arm of the validation is fitted.
- `data$taxonomy`: data frame with species as row names and ranks as columns, broadest first. Requires `b_prior = "grouped"`; whether supplying a taxonomy should switch the prior on by itself, as `embeddingVars` does, is settled in part 2's plan. Missing values are allowed only from some rank downwards. Species absent from the table are an error.
- `listPriors$sigma_b_scale`, `sigma_emb_scale`, `axis_scale`, `taxon_scale`: positive numbers, with the defaults above. Supplying one that does not apply is an error, following `read_b0_prior()`.
- `"grouped"` with no `embeddingVars` and no taxonomy is valid: it gives the learned measured-group scale and the community-average response.
- The feature has its own help page (Doug's decision, 12 November 2026), opened by `?groupedPrior` or `?embeddingVars`. It holds the settings, the outputs, the accessors and the reading of a management covariate. `?runOccJSDM` gains only two sentences that point to it, and the package site lists the page and the accessors in a section of their own, marked experimental.
- The resolved settings, the group of every covariate column and the taxon of every indicator column are saved in `infos`.

## Outputs

- `sigmab_output` keeps its shape and holds the measured-group scale.
- New, present only under `"grouped"`: the embedding leftover scale, the axis scales, the rank scales, and the community-average row of `G`.
- Taxon effects are not stored by default, because they can run to hundreds of rows per covariate. `listParams$store_taxon_effects = TRUE` adds them to `G_output`. Predictions need only `B`.
- `thinOutput()`, the chain-combination code and the diagnostics handle the new arrays.
- `returnCovariateScales()`: posterior summaries of the group, axis and rank scales.
- `returnEmbeddingBlends()`: for each draw, take the low-rank part of the embedding coefficients and reduce it by singular value decomposition to ordered blends with a fixed sign convention, as `reparamFactorModel()` does for the factors. Returns blend weights and site scores, with a documented warning that blends are defined up to rotation.
- `returnTaxonomicSignal()`: for each covariate group, the share of among-species variance in responses attributed to each rank, the real traits, the latent traits and the leftover. It is computed per draw from quantities that are always stored, so it needs no taxon effects: the squared rank scales, the squared group scale, the variance across species of the real-trait part from the trait rows of `G`, and the sum of squared entries of the covariate's column of `C`. The axis scale multiplies every component of a covariate equally, so shares are formed per covariate and then averaged within the group, weighted by the squared axis scale so that unused axes do not count. The community average is a mean, not a variance, and is excluded.
- Variation partitioning reports measured covariates and the embedding as separate shares.
- Convergence is assessed on quantities that do not rotate: the scales on the log scale, `B`, and the embedding's contribution to the linear predictor. `A` and `C` are excluded.
- Two of the outputs above were built differently in part 1, as Task 7 of the plan records. Variation partitioning is left as it is, because its three-way partition cannot hold a fourth share, and the split between measured covariates and the embedding is given by `returnEmbeddingShare()`. The package's diagnostics report R-hat on the scales as stored; the log-scale version is left to the validation study.

## Simulator

`simulateOccJSDMData()` gains optional inputs to generate truths for the study: an embedding block, a chosen number of true blends over a chosen number of active axes, separate leftover scales for the two groups, and a nested taxonomy with a scale per rank. The embedding block is either supplied by the caller as a matrix, which is how the study passes rows of the real extract, or simulated with correlated columns from a factor model with a decaying spectrum, which is what the unit tests use so that they need no network. Existing calls are unaffected. A management covariate correlated with a true blend is built by the study harness, not the simulator.

## Real embedding extract

The study's embedding design matrix is real, not simulated, because real embedding columns have a correlation structure that can only be guessed at. Doug decided on 11 November 2026 that the embeddings are to be downloaded as part of this work, and that the download is a step of this design, not something done ahead of it. Nothing has been downloaded.

What is extracted:

- **AlphaEarth Foundations**, the annual Satellite Embedding V1: 64 columns named `A00` to `A63`, 10 m pixels, years 2017 to 2025. This is the main design matrix.
- **Tessera**: 128 columns, 10 m pixels, read at the same points. It is used for one confirmatory set, because 128 columns on about 250 sites is the harder case. One dataset version and variant is chosen and recorded, since versions are separately learned feature spaces and must never be mixed.
- **OlmoEarth**, added on 12 November 2026 after Doug pointed to it: the v1.3 annual embeddings from the Allen Institute for AI, released on 22 October 2026. 128 columns, of which the first 64 work on their own as a smaller embedding, 10 m pixels, built from Sentinel-1, Sentinel-2 and Landsat. The years 2023 to 2025 are global and 2017 to 2022 cover "select countries"; whether Norway is among them is not known. It is read at the same points and saved with the other two, but it is not used in the simulation study, which needs only a realistic correlation structure among columns and gets that from AlphaEarth. Its use is in project D, where the three embeddings are compared on the real data. Its nested widths allow a like-for-like comparison with AlphaEarth at 64 columns and with Tessera at 128. One version is chosen and recorded, on the assumption that v1, v1.2 and v1.3 are separate feature spaces.
- **Points:** about 6,000 random points in the study landscape, of which the forest points are kept, using ESA WorldCover 2021 class 10 (tree cover). Plots are in forest, and water or farmland points would add contrasts to the embedding that the plots do not have. Each replicate draws its fitted and held-out sites from the forest points without replacement.
- **Footprint:** both the point's own pixel and the average of the 3 by 3 pixels around it are saved. The average is closer to a 250 square metre plot and is the one the study uses.
- **Extent and year:** the box used so far, longitude 10.80 to 11.10 and latitude 60.30 to 60.55, and the year 2024 are Claude's guesses at "west and north of Lake Hurdal". They are replaced by the real extent of the BioDivAbove plots and their sampling year when Doug has them.

How it is fetched, checked on 11 November 2026 against the providers' documentation:

- AlphaEarth is public on Source Cooperative as cloud-optimised GeoTIFFs under `https://data.source.coop/tge-labs/aef/v1/annual/`, one folder per year and UTM zone, with an index of file footprints in `aef_index.parquet` (78 MB). No Earth Engine account is needed. Only the window covering the box is read. The files are stored bottom-up and are opened through their companion `.vrt` files. Values are signed 8-bit integers with -128 as missing, converted with `((v / 127.5)^2) * sign(v)`, after which each pixel's vector has unit length.
- Tessera is fetched with the `geotessera` Python package from the same host, as 0.1 degree tiles of about 80 MB each at this latitude, roughly 1 GB for the box. An unofficial R port, `GeoTessera`, exists on GitHub and has not been tried.
- OlmoEarth v1.3 is on the same host, at `https://source.coop/olmoearth/olmoearth-v1-3-embeddings`. This was checked on 12 November 2026 against the provider's [release post](https://allenai.org/blog/olmoearth-v1-3) only. The post does not give the file format, the data type or any quantisation, so the builder reads the dataset's own documentation before writing a reader. A browser tool, the Embeddings Explorer, exports GeoTIFFs and is a fallback for a small box. The licence is the "OlmoEarth Artifact License"; its terms have not been read, and the builder reads them and records any attribution it requires before the extract is committed. The provider reports that its 128-column version scores highest on 20 of 27 of its own evaluations against AlphaEarth and Tessera v2; whether any of those tasks resembles species distribution modelling has not been checked.
- Doug asked whether this can be done from R. `terra` should be able to read the AlphaEarth files directly from their web addresses; this is untested. The builder tries that first and falls back to Python.
- Two draft Python scripts are in `dev/simstudy/embedding-extract-20261111/`. They compile and have never completed a run. They are a starting point, not a tested tool.
- The extract needs network access to `data.source.coop`, to the WorldCover bucket on Amazon S3 and, for the Python route, to `pypi.org`. On 11 November 2026 a Claude session could not reach the first and last of these at first and could later. Raw downloads go to a scratch directory outside the repository and are not committed. The extracted tables, the points and a provenance file are committed with the study.
- The AlphaEarth licence (CC-BY 4.0) requires the attribution "The AlphaEarth Foundations Satellite Embedding dataset is produced by Google and Google DeepMind", which goes in the provenance file and in anything published from the study.

Before the extract is used, the builder reports its basic properties to Doug: the number of forest points, the spectrum of the column correlation matrix, and how many principal components carry 90% of the variance. Those numbers also set the simulated matrix used in the unit tests.

## Validation

A simulation study under `dev/simstudy/`, following that directory's conventions, sized to BioDivAbove: about 250 sites, 64 embedding columns, 4 measured covariates, 100 species with a rare-dominated spread of prevalence, `binary`. A confirmatory set uses 300 species. A smaller `two_stage` set (about 100 sites, 30 species, 3 field samples, 3 PCRs) covers the Gaoligongshan case.

Fitted arms:

- *measured only*: the fixed prior on the measured covariates, with no embedding;
- *raw*: the fixed prior on the measured covariates and all raw embedding columns;
- *components*: the fixed prior on the measured covariates plus the first 5 or 10 principal components of the embedding;
- *grouped*: the grouped prior, covariate side;
- *grouped with taxonomy*, in truth 5 and the re-check of truth 4.

Every arm runs 2 chains of 5,000 burn-in and 5,000 kept iterations, the length at which `B` reached nominal coverage in the October studies. Every arm uses the same number of latent traits.

Truths, with starting strengths. A short pre-study tunes them so that the *measured only* arm has a realistic held-out discrimination (an AUC of about 0.7 to 0.8), and they are then fixed and recorded before any main run.

1. Low rank: species respond to 3 blends of 10 active axes. The embedding and the measured covariates each contribute a variance of about 1 on the logit scale, and leftovers carry about 10% of the embedding's share.
2. No blend structure: each species responds to 3 axes of its own, with coefficients of standard deviation 0.5. This is the unfavourable case.
3. Null: the embedding has no effect.
4. Confounded management: truth 1 plus a binary management covariate correlated at about 0.5 with a true blend. In half the replicates its true effect is zero for every species. In the other half the community-average effect is 0.5 and it increases by 0.5 per standard deviation of one trait.
5. Taxonomy: truth 1 with the measured-covariate responses built from a nested taxonomy of about 5 orders, 20 families and 50 genera. In the conserved version the rank scales are 0.2, 0.5 and 0.2 with a species leftover of 0.2. In the unconserved version all rank scales are 0 and the leftover carries the same total variance.

The embedding design matrix is the real AlphaEarth extract described under *Real embedding extract*, with columns standardised as `runOccJSDM()` standardises any covariate. A confirmatory set repeats truths 1 and 4 with the 128-column Tessera extract, for the *components* and *grouped* arms only.

How comparisons are made. Arms are fitted to the same simulated datasets, so every comparison is a paired difference per replicate, summarised by its mean and 95% interval across replicates. "Better" means the interval lies wholly on the better side of zero. "Not worse" means the interval excludes a loss larger than the tolerance, which is 1% of the *measured only* arm's mean log loss per site-by-species cell. The lowest prevalence band is species with true prevalence under 5%.

Criteria. Every table is reported overall and by band of true prevalence. A criterion must hold both overall and in the lowest band.

- **Prediction:** log score at 250 held-out sites drawn from the same truth. *Grouped* must be not worse than *components* in truths 1 to 4 and better in truth 1. In truth 3 it must be not worse than *measured only*.
- **Rarity:** the signed and absolute error of predicted prevalence at the held-out sites. In the lowest band the absolute error under *grouped* must be no more than 10% larger than under *components*, in truths 1 to 4.
- **Management effect:** assessed in truth 4 at three levels: each species' coefficient, the community-average row of `G`, and the trait slope. The true community average is the mean of the species' true effects. For each level, coverage of the 95% interval must be at least 0.93 and no more than 0.02 below the *components* arm. Under a true effect of zero, species' intervals must exclude zero in no more than 7.5% of species-replicates. If the *components* arm itself falls below 0.93, that is reported as a finding about the fixed prior and *grouped* is judged on the relative rule alone.
- **What truth 4 does and does not test:** in truth 4 the management covariate is generated from a blend of the real embedding, so there the embedding is a confounder and the total and direct effects are the same number. Had management changed the blend instead, the model for the species given both would be the same, and a fit cannot tell the two apart. So truth 4 tests whether the prior recovers the management coefficient with the embedding held fixed, which is what either reading needs. It does not test the causal assumptions, which no choice of prior can repair. The *measured only* arm's management coefficient in truth 4 is reported alongside, as the size of the error from leaving the control out.
- **Taxonomy (part 2):** in the conserved version of truth 5, *grouped with taxonomy* must be better than *grouped* on held-out log score in the lowest band. In the unconserved version it must be not worse. The management criterion must still pass on truth 4 with a taxonomy supplied.

Convergence. A fit fails if R-hat exceeds 1.1 on any group, axis or rank scale (on the log scale) or on the management coefficients. Criteria are computed on all fits and again on converged fits only, and both are reported. If more than 10% of *grouped* fits fail in any cell, the prior fails on mixing grounds whatever the other criteria say.

Also recorded: coverage and bias of the other measured coefficients, effective sample sizes, and wall time per fit against the fixed prior.

Replicates and budget. A first pass uses 30 replicates per cell, and the cells that decide a criterion are confirmed at 100. The pre-study times one fit per arm. If the first pass is projected to take more than 48 wall-hours on five workers, the builder stops and reports to Doug before launching it.

The `two_stage` set uses truths 1 and 4 and the *components* and *grouped* arms, keeps all 64 embedding columns, and applies the prediction and management criteria. Its numbers of sites and species above are placeholders, to be replaced by those of the Gaoligongshan dataset.

What happens on failure. If part 1 fails a criterion at 100 replicates, the builder reports to Doug and part 2 is not started. In `occupancy` and `two_stage` fits the community-average row and the taxon effects reach the data only through the latent occupancy states, the structure that made `mu_b0` fail R-hat until the collapsed shift in `R/b0_shift.R` was added. If they fail to converge for that reason, a new sampler move is not part of this project: the builder reports and stops.

## Staging and tests

Part 1 (covariate side) is built, tested and validated on truths 1 to 4 before part 2 (taxonomy) is started. Part 2 starts only when part 1 has passed the prediction, rarity and management criteria at 100 replicates. Part 2 is validated on truth 5 and re-checked on truth 4.

Unit tests cover: bit-identical fits with `b_prior = "fixed"` against a fit saved before the change, for all four model types; each new conditional update against a direct calculation; input validation; chain combination and thinning of the new arrays; one short spatial fit and one parallel-chain fit under `"grouped"`.

## Constraints

- Work happens in the worktree `.worktrees/grouped-covariate-prior` on branch `claude/grouped-covariate-prior`, including every run directory, log and result.
- At most five fitting workers across all concurrent runs on the laptop. Long runs are launched detached. Runs use R 4.6.1, and the power mode is recorded with any timing.
- Paired refits outside the study runner use the runner's RNG stream.
- Markdown follows `CLAUDE.md`.
- How the work reaches `main`, and whether Alex reviews it as a pull request, is Doug's decision once part 1 is validated.

## Who decided what

Doug decided: the purpose and the three reasons for doing this in occJSDM; that the embedding sits alongside measured and management covariates; that the management conclusion wanted is an effect beyond what the satellite sees, and (12 November 2026) that he reads this as covariate adjustment for causal analysis and wants the total and direct effects, the extra assumption and the with-and-without comparison written into the design; that BioDivAbove is first and is a pure JSDM; that Gaoligongshan, two-stage, follows; that taxonomy is part of this design; that the study uses downloaded AlphaEarth or Tessera embeddings, and that the download belongs in this design and is not done ahead of it; and (12 November 2026) that OlmoEarth is added as a third source in the extract; and (12 November 2026) that naming `embeddingVars` is enough to turn the grouped prior on, after asking how to keep the feature from complicating standard use, and that it has its own help page; and (12 November 2026) that the default `axis_scale` is `0.7 / K`, after the prior check showed the first default to be far too wide, with a half-t on the axis scales the alternative he did not choose.

Claude's calls, for Doug to check:

- one switch, `b_prior`, in place of several, and one fitting function in place of two;
- half-Cauchy priors, with the size of the embedding carried by the axis scales and the blend weights scaled by the number of blends;
- the first default of 1 over the square root of K for the axis scale, which the prior check showed to be wrong;
- after the whole-branch review: starting the embedding columns of `C` at their prior scale under the grouped prior, having the shrinkage test assert only on the strongly active columns, and adding one pointer sentence for `embeddingVars` under `listParams` in `?runOccJSDM`;
- the axis scale multiplying the whole row, including `G` and `C`;
- latent traits shared between the two groups;
- the community-average row added only under `"grouped"` and only on the covariate side;
- one rank scale per rank and covariate group, updated rank by rank;
- taxon effects not stored by default, and the taxonomic signal computed from the scales;
- the argument and function names;
- the study's sizes, truths and strengths, chain length, tolerance and pass thresholds, including that a criterion must hold in the lowest prevalence band as well as overall;
- judging the management effect at three levels, absolutely and relative to the *components* arm;
- everything under *What the management coefficient estimates* beyond Doug's reading of it as covariate adjustment: the direct and total effects, the two fits compared on the probability scale, the overlap and matched checks, and the classification of covariates by timing;
- not simulating a truth in which an unmeasured cause affects both the embedding and the species. Its bias falls on every arm alike, so it would measure the assumption and not the prior. It can be added as an illustration if Doug wants one;
- the 48-hour rule for the first pass, taken from the rule Doug set for the Stage 2 comparison;
- reporting and stopping, not building a new sampler move, if the new community-level parameters fail to mix in replicated fits;
- building taxonomy after the covariate side;
- for the extract: AlphaEarth as the main design matrix with Tessera as a confirmatory set, OlmoEarth extracted but kept out of the simulation study and compared on the real data in project D, the guessed box and year, keeping forest points only, using the 3 by 3 average, and trying R before Python.

## Left to the builder

How the vector of prior standard deviations is threaded through the R and C++ signatures; how embedding names map to design-matrix columns; the layout of the new output arrays; how the simulator generates correlated embedding columns and a nested taxonomy; how variation partitioning is split by group; the study harness; test file layout.

## The other four projects

- **B, thousands of species.** BioDivAbove has up to about 4,500 OTUs. Needs optional storage of the largest arrays, a rule for OTUs with too few detections (taxonomy may let rare OTUs stay in), and time and memory measurements. Sized once the real files show how many OTUs are modellable.
- **C, prediction to planning units.** Joint posterior draws across species and units, in chunks, plus a flag for units unlike any sampled site. Needed when a study reaches the prioritisation step.
- **D, BioDivAbove baseline.** The principal-components arm on the real data with the current release, as a benchmark, run for each of the three embeddings (AlphaEarth, Tessera, OlmoEarth) so that they are compared on held-out prediction for real species. Where a management covariate is included, it follows the reporting list under *What the management coefficient estimates*. Waits for the data.
- **E, occupancy from the beetle counts.** Recorded here so that it is not lost; it is not part of project A, and nothing in A depends on it.
  - *The idea, Doug's:* the BioDivAbove beetles are counts of individuals per window trap (609 species, about 90,000 individuals), with no replicate samples. Divide each trap's catch at random into several subsamples and treat them as replicate samples, so that the occupancy model can be fitted and false negatives estimated, where the other three datasets can only be fitted as a pure JSDM.
  - *Provenance:* this entry is written from Doug's one-line summary on 11 November 2026. His earlier, fuller message did not reach the conversation, and the entry should be corrected against it.
  - *When it is sound:* if individuals arrive at a trap independently of one another, a random division of the catch into K parts gives K independent catches, each equivalent to running the trap for 1/K of the time. An occupancy model fitted to the parts then estimates the chance that a whole trap misses a species that is present. This holds however the catch is spread over the season, because the division is random and not by date. Doug's judgement on 12 November 2026 is that independent arrival is a reasonable assumption for these traps, and Claude agrees. An earlier version of this entry said beetle catches are "usually clumped", which was imprecise: what is well established is that counts differ between traps more than chance alone would give, and that is variation in abundance between plots, covered in the next item, not beetles arriving in groups. The exceptions are species that swarm.
  - *What it corrects, and where it falls short:* where a species is scarce, a trap often misses it. Estimating that miss rate is the purpose of the method, not a fault in it. The limit is that the occupancy model gives a species one detection probability, adjusted only by detection covariates, whereas the real one rises with local abundance and so differs from plot to plot. The model learns detection from the plots where the species was caught, which are the plots where it is commoner, and applies it to the plots with empty traps, where the species, if present, is scarcer and easier to miss. So the correction goes in the right direction but is too small, and occupancy stays underestimated, by more the more abundance varies among occupied plots. Entering the environmental covariates in the detection stage absorbs the part of that variation they explain. If abundance grades continuously down to zero, no method can separate "very scarce" from "absent", and the answer then rests on the model's assumptions.
  - *Total catch as a collection covariate, Doug's addition on 12 November 2026:* the total number of beetles in a trap, of all species, on the log scale, entered as a collection covariate. It measures how well that trap caught, whatever the cause (how open the stand is, exposure, how long the trap ran), so it absorbs the part of the plot-to-plot variation in detection that all species share. It also answers the catchability point under *What the management coefficient estimates*. It does not absorb a single species being locally commoner or scarcer than the total would suggest, so the shortfall above gets smaller but does not vanish. Two details to settle: for a dominant species the total is partly its own count, which is circular, so the total may need to leave out the few dominant taxa; and the slope should be allowed to differ by species.
  - *What it estimates:* a species with one individual in a trap appears in one part of K, which the model reads as low detection, and so it infers misses at plots with empty traps. "Occupied" then means "would be caught by a longer trapping period", which has to be stated wherever the maps are used.
  - *How it relates to other models:* it is not N-mixture modelling, which needs real repeat counts and learns detection from how counts vary between independent visits. A single catch only identifies the product of abundance and detection. What subsampling can do is tell an absent species from a present but uncaught one, using how large the catches are where the species was caught. That is the logic of a zero-inflated count model, and such a model fitted directly to the counts uses all of the information and needs no choice of K. The Royle-Nichols model is the closest relative that infers abundance from detections, and it too needs real replicates.
  - *What is in its favour:* occJSDM cannot fit counts, but it can fit occupancy with replicate samples. Subsampling is therefore a way to get inference of the zero-inflated kind, with covariates, traits and the joint species structure, from the sampler that already exists. The zero-inflated model it stands in for is one model fitted to the whole catch, not one per subsample. The information lost by subsampling shrinks as the number of parts grows, because with many small parts the number of non-empty parts approaches the count itself. So it should be done with many parts, not two or three.
  - *Recommendation, Claude's, 12 November 2026:* fit the occupancy model on the subsamples; do not replace it with a zero-inflated count model. With many parts and a constant detection probability, the subsampled occupancy model is close to a zero-inflated Poisson model, so little is given up. What a count model could add is an allowance for abundance varying among occupied plots beyond what covariates explain (the negative binomial version), which is the shortfall described above. Against that, no count model Claude knows of offers the traits, the taxonomy, the embedding prior and the joint posterior draws that are the reasons for using occJSDM. The count model's place is as a check: in the simulation below, and on a subset of the real species, to measure how large the shortfall is.
  - *If occJSDM gains count data:* Doug said on 12 November 2026 that Alex plans to add it. The recommendation then depends on what is built. A zero-inflated model (an occupancy state, plus a count where the species is present) that carries the same traits, taxonomy and embedding prior would replace subsampling for the beetles: it uses the whole catch, needs no choice of K and no random division, and in its negative binomial form allows abundance to vary among occupied plots. A plain count model with no occupancy state would not replace it, because it estimates expected catch and has no false negatives. The count branch sketched in the code today is of that plain kind (`R/jsdmfun.R` line 775, a negative binomial with mean `exp(eta)`). One caution carries over: zero inflation and negative binomial dispersion both produce zeros and are hard to tell apart, so the occupancy estimate would depend on the count distribution assumed. Until a count model exists, subsampling works with the current release, and afterwards it is the check on the new model, since the two should agree when the parts are many and the count model is Poisson.
  - *What "an occupancy state" means:* a hidden yes or no for each species at each plot: is the species there at all? occJSDM already has one, the latent presence `z` that the `occupancy` and `two_stage` models estimate.
    - *A count model with an occupancy state* has two parts. First, is the species present at the plot, with a probability that depends on the occupancy covariates. Second, if it is present, how many the trap catches, with an expected value that can depend on its own covariates, such as total catch.
    - *A zero then has two possible causes:* the species is absent and the trap could never have caught it, or the species is present and the trap caught none, which is a false negative. The model works out how likely each is for every empty trap, from the size of the catches where the species was caught. With an expected catch of 5 where present, a Poisson catch is empty under 1% of the time, so an empty trap almost certainly means absent. With an expected catch of 0.5, it is empty about 60% of the time, so many empty traps are misses. This is the reasoning of the subsampling idea, done directly on the counts.
    - *A count model without one* has only the second part. Every plot gets an expected catch, which is always above zero and is just smaller in poor habitat. There is no "absent", and a zero is never a false negative: it is what a low expected catch often produces.
    - *Why it matters here:* the two give different maps. With an occupancy state the map is the probability that the species is present; without one it is the expected number caught. Only the first corrects for false negatives, which was the reason for the subsampling idea.
    - *The question for Alex:* whether the count model keeps `z` and models the count only where `z` is 1, or models the counts directly.
  - *Published precedents, searched on 12 November 2026:* eight web searches found no paper that divides one count at random into parts and fits an occupancy model to them. That is a limited search, not proof that none exists; Doug's Zotero library was not searched. The close relatives are:
    - *Subsamples of one collection as replicates:* [Reitan, Ergon and Liow 2023](https://doi.org/10.1017/pab.2022.17), *Paleobiology* 49: 23 to 38. Fossil shells from one site are the replicates in a multispecies occupancy model used to estimate relative abundance. The parts are physical units, not a random division of individuals. Abundance varying between sites is handled with a beta-binomial for the number of detections, which is a remedy for the shortfall described above. Read at source.
    - *PCR replicates in eDNA occupancy models:* [Schmidt et al. 2013](https://doi.org/10.1111/2041-210X.12052), *Methods in Ecology and Evolution* 4: 646 to 653. PCR replicates are random subsamples of one DNA extract, treated as replicates. occJSDM's own `two_stage` model follows this convention. Title and DOI checked; the paper was not read.
    - *One continuous record cut into occasions:* [Pautrel et al. 2024](https://doi.org/10.1111/2041-210x.14314), *Methods in Ecology and Evolution*. Camera-trap records are routinely cut into occasions, by time and not at random. This comparison of discrete-time and continuous-time occupancy models reports that the length of the interval had minimal effect on occupancy estimates, which is the analogue of the choice of K. Abstract and preprint excerpts only.
    - *Count splitting in statistics:* [Neufeld et al. 2022](https://doi.org/10.1093/biostatistics/kxac047), *Biostatistics*, for Poisson counts, and [Neufeld et al. 2023](https://arxiv.org/abs/2307.12985), a preprint, for negative binomial counts. They give the theory: parts of a Poisson count are independent, and parts of a negative binomial count are correlated. The method is used to validate models, not to make replicates. Preprint excerpts only.
    - *Single-visit models, which attack the same problem with covariates and no replicates:* Lele, Moreno and Bayne 2012, *Journal of Plant Ecology* 5: 22 to 31, for occupancy, and Sólymos, Lele and Bayne 2012, *Environmetrics*, for zero-inflated counts with detection error. Both are in the R package `detect`. Titles checked; the papers were not read.
  - *What to settle first:* whether the window traps were emptied more than once with the collections kept separate, which would give real replicates and make subsampling unnecessary; and whether each plot had more than one trap.
  - *If it goes ahead:* a simulation with known occupancy, independent arrivals, and abundance that varies among occupied plots by a controlled amount (from none to strong, partly explained by covariates and partly not), comparing subsampled occupancy fits at several K against a pure JSDM on presence and absence and against a zero-inflated count model, scored by occupancy error by species rarity. occJSDM does not fit counts today (the `TODO.md` item *Count-data models*), so the count arm needs either that work or another package.
