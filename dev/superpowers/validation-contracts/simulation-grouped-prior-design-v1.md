# Simulation study of the grouped covariate prior: design contract, version 1

Issued by occJSDM on 17 November 2026, on draft PR 41 (grouped covariate prior, package 0.3.2). This is the package spec's Validation section (`dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md`) turned into a directive, adapted to run on a dataset that has real sites, real measured covariates and a real embedding, with species simulated on top. Version 1 runs on AlphaEco only. It is to be run before the real-data check (`real-data-grouped-prior-design-v1.md`), because a known truth is the only way to show that the prior recovers a signal.

## Purpose

To show, before anyone relies on the grouped prior, that it passes three criteria against known truths: prediction (it predicts no worse than fixing the blends in advance, and better when the truth is low rank), rarity (it does not do worse for rare species), and the management effect (a coefficient fitted beside the embedding is recovered with honest intervals). A fourth criterion, convergence, applies throughout. The real-data check can then be read as a test of behaviour, with this study as the evidence that the prior works when there is something to find.

## Who does what

- The dataset section builds the harness, runs it and keeps everything, in its own repository. BioDivAbove's `simulation/` folder is a template for the harness, not a place to run it.
- occJSDM reads the output (see the output contract) and decides what a failure means for the package.
- The dataset section changes no package code. A criterion that fails is reported to Doug, who passes it to occJSDM, and part 1 is not extended.

## The design matrix is real

- Sites: the real survey sites (for AlphaEco, 320 ponds), with their real coordinates.
- Measured covariates: the real measured covariates (for AlphaEco, the five pond variables, transformed and standardised as in the real-data contract). They are real so that their correlation with the embedding is real, which is what the confounding truth depends on.
- Embedding: OlmoEarth v1.3 at 64 bands, at the site scale, 64 columns, standardised on the training sites. A confirmatory set uses the two scales together, 128 columns, for the components and grouped arms only. It is the harder case, as Tessera's 128 columns were in the package spec.
- Species are the only simulated part: a new set for every replicate, from a truth, on this fixed design.
- A consequence to state in the result: every replicate shares the same sites and covariates, so variation across replicates comes only from the species and the split. The findings are conditional on this design, and they say how the prior behaves on embeddings of this kind, not on every landscape.

## Truths

Every truth has S species with a rare-dominated spread of prevalence (for example Beta(0.5, 8) floored at 0.01, recorded with the realised spread) and two simulated traits per species. Species responses to the measured covariates depend on the traits plus a species leftover. Each replicate has the same number of hidden environmental covariates as the real-data contract (two), with species loadings of standard deviation 0.5.

1. Low rank. Species respond to 3 blends of 10 active axes chosen at random from the 64 columns. The embedding and the measured covariates each contribute a variance of about 1 on the logit scale to the linear predictor. Leftovers on the active axes carry about 10% of the embedding's share.
2. No blend structure. Each species responds to 3 axes of its own, with coefficients of standard deviation 0.5. This is the unfavourable case for the prior.
3. Null. The embedding has no effect.
4. Confounded management. Truth 1 plus a binary management covariate correlated at about 0.5 with a true blend (built from a score on the real embedding, dichotomised at the training median). In half the replicates its true effect is zero for every species. In the other half the community-average effect is 0.5 and it increases by 0.5 per standard deviation of the first trait. The harness builds this covariate, because the dataset has none. In truth 4 the embedding is a confounder, so the total and direct effects coincide there, and the truth tests recovery of the coefficient with the embedding held fixed. It does not test causal assumptions.

A pre-study tunes the strengths so that the measured-only arm has a realistic held-out discrimination (an AUC of about 0.7 to 0.8). The strengths are then fixed and recorded before any main run. Taxonomy truths belong to part 2 of the package design and are not in version 1.

## Split and replicates

- Each replicate draws one blocked split on the real coordinates: blocks of the same side as the real-data contract, five folds, one fold held out at random with a recorded seed. For AlphaEco that fits 256 ponds and holds out 64. (The package spec asked for 250 held-out sites, which this dataset cannot supply alongside 250 fitting sites.)
- All arms are fitted to the same simulated dataset and the same split. Every comparison is a paired difference per replicate, summarised by its mean and 95% interval across replicates.
- A species seen fewer than 3 times in the fitting sites is dropped from every arm, since a fit cannot learn from none. The number dropped is recorded.

## Arms

- measured only: the fixed prior on the measured covariates (and the management covariate in truth 4), no embedding.
- raw: the fixed prior on the measured covariates and all embedding columns.
- components: the fixed prior on the measured covariates plus the first 5 and, separately, the first 10 principal components of the embedding, computed on the training sites and projected onto the held-out sites.
- grouped: the grouped prior, covariate side, with the package defaults, including the axis scale of 0.7 over K.

Every arm uses 2 chains of 5,000 burn-in and 5,000 kept iterations, the length at which `B` reached nominal coverage in the October studies, and the same number of latent traits. The same fixed settings as the real-data contract apply otherwise.

## Criteria

Every table is reported overall and by band of true prevalence. The lowest band is species with true prevalence under 5%. A criterion must hold both overall and in the lowest band. "Better" means the interval lies wholly on the better side of zero. "Not worse" means the interval excludes a loss larger than the tolerance, which is 1% of the measured-only arm's mean log loss per site-by-species cell.

- Prediction. The log score at the held-out sites. Grouped must be not worse than components in truths 1 to 4 and better in truth 1. In truth 3 it must be not worse than measured only.
- Rarity. The signed and absolute error of predicted prevalence at the held-out sites. In the lowest band the absolute error under grouped must be no more than 10% larger than under components, in truths 1 to 4.
- Management effect, in truth 4, at three levels: each species' coefficient, the community-average row of `G`, and the trait slope. The true community average is the mean of the species' true effects. For each level, coverage of the 95% interval must be at least 0.93 and no more than 0.02 below the components arm. Under a true effect of zero, species' intervals must exclude zero in no more than 7.5% of species-replicates. If the components arm itself falls below 0.93, that is reported as a finding about the fixed prior and grouped is judged on the relative rule alone.
- Convergence. A fit fails if R-hat exceeds 1.1 on any group, axis or rank scale (on the log scale) or on the management coefficients. Criteria are computed on all fits and again on converged fits only, and both are reported. If more than 10% of grouped fits fail in any cell, the prior fails on mixing grounds whatever the other criteria say.

Also recorded: coverage and bias of the other measured coefficients, effective sample sizes, and wall time per fit against the fixed prior.

## Pre-study, replicates and budget

- The pre-study times one fit per arm on truth 1 and checks the measured-only AUC under every truth. It tunes the truth strengths, then fixes them. It also fixes the number of species S. Using the real data's 48 species costs about half as much per fit as the package spec's 100.
- A first pass uses 30 replicates per cell, and the cells that decide a criterion are confirmed at 100. A cell is a truth (truth 4 split into its two effect halves) by arm.
- The pre-study's timings project the first pass. At 30 replicates it is about 600 fits, one per arm per replicate per truth with truth 4 doubled, and at a measured cost per fit that is likely many days on five workers. If the projection is over 48 wall-hours on five workers, the dataset owner stops and reports to Doug before launching. The usual cuts, in order, are fewer species, fewer replicates in the first pass (10), truths 1 and 3 first, and dropping the 10-component arm. Any cut is recorded as a deviation.
- At most 5 fitting workers at once, single-threaded fits, and nothing is launched at full scale without the dataset owner's go-ahead.

## Things to know before building

- The fits use the installed, pinned package, never `devtools::load_all()`. See the real-data contract's list of earlier mistakes, which all apply.
- The package's simulator `simulateOccJSDMData()` (`R/simulateData.R`) accepts a supplied embedding matrix and builds low-rank truths. BioDivAbove's `simulation/helpers.R` builds all four truths directly and is a working template.
- Under the grouped prior the package starts the blend weights at their prior scale, a tight start. Whether overdispersed starts behave better is open and cannot be tested without a package change, so version 1 uses the default. If the pre-study shows axis-scale R-hat failures, say so in the result.

## Optional module

- A two-stage set (about 100 sites, 30 species, 3 field samples and 3 PCRs) for the Gaoligongshan case, on truths 1 and 4 and the components and grouped arms. Not in version 1, because AlphaEco has no replicates.

## What happens on failure

If part 1 fails a criterion at 100 replicates, the dataset owner reports to Doug and nothing further is built on the prior. In `occupancy` and `two_stage` fits the community-average row and any taxon effects reach the data only through the latent occupancy states, which is what made `mu_b0` fail R-hat until the collapsed shift in `R/b0_shift.R` was added. If a failure traces to that, a new sampler move is not part of this study: report and stop.

## Versioning

A new version is issued by occJSDM. The dataset section writes the version it followed, and any deviation with its reason, in the run record.
