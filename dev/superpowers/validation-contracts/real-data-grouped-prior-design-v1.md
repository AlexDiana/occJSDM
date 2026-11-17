# Real-data check of the grouped covariate prior: design contract, version 1

Issued by occJSDM on 17 November 2026, on draft PR 41 (grouped covariate prior, package 0.3.2). Derived from the AlphaEco shakedown design, which is the worked instance. The package's own design is `dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md`, and the simulation validation study described there is separate and unchanged.

## Purpose

To learn whether the grouped prior behaves sensibly on a real community, where there is no known truth. This is a shakedown, not a validation, so the study looks for warning signs and answers a set of questions about prediction. It does not produce pass marks. The simulation study is what tests recovery of a known truth.

The expected result on a dataset with a weak habitat signal is that the prior shrinks the embedding towards nothing and predicts no worse than the measured covariates alone. That shows the prior does no harm when there is little to find. It cannot show that the prior finds a strong signal.

## Who does what

- The dataset section runs everything and keeps all data, fits and per-unit results in its own repository.
- occJSDM reads the evaluation that comes back (see the output contract) and decides what it means for the package.
- The dataset section changes no package code. If a fit shows what looks like a package defect, it stops and reports to Doug, who passes it to occJSDM.

## What the dataset section needs

- Detections for one community at many sites: presence or absence with one row per site, or replicated detections for a two-stage fit. The contract is written for one row per site and the optional modules say what changes otherwise.
- A small set of measured covariates (for example five), chosen by the dataset owner and recorded.
- At least one embedding source (AlphaEarth, TESSERA or OlmoEarth) extracted at two scales: a pond scale or site scale, the mean within a small radius of the sampling location (for AlphaEco, 30 m), and a landscape scale, the mean over a much wider area (500 m). Both together make a third set. The radii are the owner's choice and are recorded.
- A grouping of species into a few groups, if the data have one, used as a trait and to break results down.
- Site coordinates, used only to build spatial folds.

## Embedding rules

- A pixel that is zero in every band is missing data written as a number. It is left out of every mean. Check for zeros as well as NA in any extract.
- Check the extract before relying on it: the vector length of a single pixel (about 1 for AlphaEarth), whether any site has little valid area at a scale, and that the extraction year matches the sampling year.
- Embedding columns enter the model through `embeddingVars`, which switches the grouped prior on. The package defaults are used unchanged, including the axis scale of 0.7 over K.

## Questions

Answered first in the result, before the warning signs. Each embedding set is a source at a scale, or two scales together.

1. Absolute performance. How well do the embeddings alone, and the measured covariates alone, predict? Mean per-species AUC is called poor under 0.7, decent from 0.7 to 0.8 and good above 0.8, with the number of species in each band. Skill against the frequency-only null is one minus the arm's held-out log loss over the null's, so 0 is no better than knowing how common each species is.
2. Relative performance. Are the embeddings better or worse than the measured covariates? Embedding only against grouped with no embedding, as a paired difference in held-out log loss and in per-species AUC.
3. Combining. Does adding the embedding to the measured covariates raise performance? Grouped against grouped with no embedding, and against embedding only. If the combination beats the measured covariates alone, sites still need measuring and the embedding is free extra information. If embedding only matches the combination, the measurement adds nothing.
4. Scale. Which scale carries the signal: site scale, landscape scale or both, overall and by species group.
5. Source. If more than one source is available, which is better at the same scale.

## Arms

All arms use the same model type, the same latent traits, the same MCMC length and the same trait table (see Fixed settings). They differ only in which covariates they are given and which prior they use.

Required, per embedding set:

- grouped: the measured covariates and all embedding columns under the grouped prior.
- embedding only: the embedding columns under the grouped prior with no measured covariate. This is what could be predicted for a site nobody has visited.

Required, once per study (they use no embedding):

- measured only: the measured covariates under the fixed prior.
- grouped, no embedding: the measured covariates under the grouped prior with no `embeddingVars`. It separates the effect of the embedding from the effect of the looser prior on the measured covariates, which grouped differs from measured only in two ways.
- the frequency-only null (needs no fit, see below).

Required, for one embedding set only (the site scale of the first source), because their job is to test the prior and not to compare scales or sources:

- raw: the measured covariates plus all embedding columns as ordinary covariates under the fixed prior.
- components: the measured covariates plus the first 5 and, separately, the first 10 principal components of the embedding, under the fixed prior. Components are computed from the training sites of each fold and projected onto the held-out sites, and the scaling is stated.

Optional: grouped with no trait, for one embedding set, which isolates the trait's own effect.

## Fixed settings

- Model type is the one the data need. Version 1 covers `binary` with one row per site and no spatial random effect. A spatial fit is a separate contract.
- Two chains of 5,000 burn-in and 5,000 kept iterations. If the full-data grouped fit fails the convergence sign at this length, that fit alone is rerun once at 20,000 and 20,000, and the rerun counts against the budget.
- The same number of latent traits in every arm, below the number of measured covariates so that the measured-only arm is not over-parameterised. Recorded.
- Two hidden environmental covariates (`n_factors` = 2) in every arm, so that any competition between them and the embedding for shared variation shows up.
- The species group enters every arm as one categorical trait, with the largest group as the reference level under the fixed prior. Under the grouped prior the trait matrix has a column of ones, so every group has a mean response.
- Strongly skewed measured covariates are transformed before standardising, and the transformation is recorded. AlphaEco used `log1p` of pond area.
- Columns are read by name, never by position.

## Held-out evaluation

- Spatially blocked 5-fold cross-validation. Sites are grouped into square blocks on the coordinate grid and blocks are assigned to folds at random with a recorded seed. The block side is chosen so that few held-out sites lie within 1 km of a training site, and the achieved distances are reported. AlphaEco used 10 km blocks and folds of 64 ponds.
- Every site is held out exactly once. Each arm gives one out-of-fold probability for every site-by-species cell, and all held-out measures come from that one table. Each arm is also fitted once to all sites, for the in-sample diagnostics.
- The frequency-only null gives each species at every held-out site the probability (detections + 0.5) over (training sites + 1), computed per fold. It is reported beside every held-out measure.
- Primary measure: log loss per site-by-species cell. It is a proper score and it is not distorted by pooling common and rare species.
- Per-species AUC, for species with at least 10 detections, counted over pairs of sites within the same fold and summed over folds. Pooling out-of-fold predictions across folds is wrong: a fold with more detections of a species has a lower training frequency, so the null scores below 0.5 when pooled.
- Every difference is paired by site-by-species cell and carries an interval from resampling blocks. Results are reported overall, by species group and by prevalence band (under 5%, 5 to 15%, over 15% of sites).

## Warning signs

Each is reported for every arm, overall and by group and band. Signs 3 to 5 compare grouped with both measured only and grouped with no embedding: the first difference is the whole change, the second is the embedding's part, and the gap between the two baselines is the prior's part.

1. Convergence. R-hat over 1.1 on the log scale for the measured-group scale, the embedding leftover scale, any axis scale, or `B`; and effective sample sizes under 400 for the same quantities. A scale with an R-hat of NA (a chain that never moved) counts as failing. Reported as the share that fail, per arm.
2. No separation among axes. The axis scales show no clear gap between a few large ones and the rest, reported as the ratio of each of the five largest to the median. For a fit with both scales, how many of the ten largest belong to each scale.
3. Extreme predictions. The share of in-sample and held-out probabilities below 0.001 or above 0.999, against measured only.
4. Worse held-out prediction. Held-out log loss worse under grouped than under measured only, or worse under any fitted arm than under the null. If no fitted arm beats the null, that is the headline, and signs 2 and 5 are read as behaviour with nothing to find.
5. Unstable measured coefficients. Community-average responses to the measured covariates that change sign or move by more than their posterior standard deviation between measured only and grouped. A move is not wrong in itself, so it is reported with the embedding share from `returnEmbeddingShare()`.
6. Implausible prevalence. Predicted prevalence departing from observed by more than a factor of 2, for species with at least 10 detections.
7. Group responses that contradict known ecology. The mean response of each species group to each measured covariate, with intervals, from the main arms, set out for the dataset owner to judge against what is known of these animals. There is no automatic rule. Responses to embedding axes are not tabulated, because an axis has no meaning of its own.

Blends are defined only up to rotation, so a blend is never read as an ecological finding.

## Things that have gone wrong before

- Site order. Before package 0.3.1, a one-row-per-site fit paired each site's detections with another site's covariates when site names did not sort in file order, which made every real fit show no signal. Pass integer site ids in file order, and stop after every fit unless the fit's design matrix matches the training covariates row for row.
- Cost of wide fits. Before 0.3.2 the end-of-fit convergence summary used memory growing with the square of the number of coefficients. It is fixed, but a full 256-column fit still stores about 1 GB of draws.
- Never time or run fits under `devtools::load_all()`: it compiles without optimisation and was about 9 times slower. Install a pinned commit from a git export into a study library.
- Pin the package commit and record it. A summary must stop if extracts mix commits, data checksums or MCMC settings.
- Check power mode before timing: Low Power Mode on battery halved sampler speed. Record it with the timings.
- Long fits started from a tool session die when that session ends. Start them detached.
- At most 5 fitting workers at once on a 32 GB laptop.

## Budget and stop rules

- Time one grouped fit at each width of embedding (for example 64, 128 and 256 columns) before the full run, then project the whole set from the measured times. AlphaEco measured 29, 87 and 265 minutes for 64, 128 and 256 columns, each step costing about 3 times the last.
- If the projection is over 48 wall-hours, drop the widest arms first and say so in the run record. If it is still over, stop and report to Doug before launching.
- Nothing is launched at full scale without the dataset owner's go-ahead.

## Optional modules

- Replicated detections (`two_stage`): the same arms, with the same sign checks on the occupancy coefficients. In two-stage and occupancy fits the community-average row reaches the data only through the latent occupancy states, so non-convergence there is a lead about `mu_b0` and the collapsed shift in `R/b0_shift.R`.
- A management or protection covariate: report a total-effect fit (pre-management covariates only) and a direct-effect fit side by side on the probability scale, with an overlap check and a matched check. The embedding is measured after management, so its coefficient is a direct effect. The package spec's section on the management coefficient says how to read it.
- Taxonomy: not part of version 1. Part 2 of the package design is not built.

## Versioning

A new version is issued by occJSDM. The dataset section writes the version it followed, and any deviation with its reason, in the run record.
