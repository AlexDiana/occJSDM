# Simulation study of the grouped covariate prior: design contract, version 2

Status: written by occJSDM (Claude, for Doug) on 21 November 2026. Not issued. It comes into force when Doug says so, and until then version 1 (`simulation-grouped-prior-design-v1.md`) is the contract in force.

Version 2 is pinned to the head of draft PR 45 (commit a2a64a1, package 0.3.2). That commit holds the grouped covariate prior (PR 41), the sampler moves, staged start and convergence check of PR 44, and the gradient functions of PR 45. Like version 1, it is the package spec's Validation section (`dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md`) turned into a directive, adapted to run on a dataset that has real sites, real measured covariates and a real embedding, with species simulated on top. It runs on AlphaEco only, and it is to be run before the real-data check.

Version 2 changes the convergence rule, adds two criteria (gradients, and the shared part), and raises the number of chains from 2 to 4. It also fits the raw arm in truths 1 and 4 only, keeps the strengths of the version 1 run in place of a new pre-study, moves the pinned commit, lets a fit's chains run on separate cores, and shortens the list of usual cuts. The rest is version 1's. It was written after the results of a version 1 run had been seen. The section "Disclosure" says what had been seen. Read it before weighing any result made under this contract.

## Purpose

To show, before anyone relies on the grouped prior, that it passes four criteria against known truths: prediction (it predicts no worse than fixing the blends in advance, and better when the truth is low rank), rarity (it does not do worse for rare species), the management effect (a coefficient fitted beside the embedding is recovered with honest intervals), and gradients (the habitat gradients are recovered when they exist). A convergence rule applies throughout. A fifth criterion, the shared part, is scored and reported and does not stop the work. The real-data check can then be read as a test of behaviour, with this study as the evidence that the prior works when there is something to find.

## What version 2 judges

- Version 2 judges only fits made after it is issued, on seeds that no earlier run used. The run record lists the seeds, and the dataset section confirms that none was used before.
- The 150 grouped fits of the version 1 run (20 November 2026, package 2b9b890) keep version 1's verdict: the prior failed on mixing grounds in every cell. The convergence draft describes how they would score under version 2. That description is not a verdict.
- The verdict of the version 1 rule is computed and reported for every fit made under version 2, beside the version 2 verdict.

## Who does what

- The dataset section builds the harness, runs it and keeps everything, in its own repository. BioDivAbove's `simulation/` folder is a template for the harness, not a place to run it.
- occJSDM reads the output (see the output contract) and decides what a failure means for the package.
- The dataset section changes no package code. A criterion that fails is reported to Doug, who passes it to occJSDM, and part 1 is not extended.

## The design matrix is real

- Sites: the real survey sites (for AlphaEco, 320 ponds), with their real coordinates.
- Measured covariates: the real measured covariates (for AlphaEco, the five pond variables, transformed and standardised as in the real-data contract). They are real so that their correlation with the embedding is real, which is what the confounding truth depends on.
- Embedding: OlmoEarth v1.3 at 64 bands, at the site scale, 64 columns, standardised on the training sites. A confirmatory set uses the two scales together, 128 columns, for the components and grouped arms only. It is the harder case, as Tessera's 128 columns were in the package spec. It is not part of the first pass counted below, and version 2 does not yet say how the gradients and shared-part criteria apply to it.
- Species are the only simulated part: a new set for every replicate, from a truth, on this fixed design.
- A consequence to state in the result: every replicate shares the same sites and covariates, so variation across replicates comes only from the species and the split. The findings are conditional on this design, and they say how the prior behaves on embeddings of this kind, not on every landscape.

## Truths

Every truth has S species with a rare-dominated spread of prevalence (for example Beta(0.5, 8) floored at 0.01, recorded with the realised spread) and two simulated traits per species. Species responses to the measured covariates depend on the traits plus a species leftover. Each replicate has the same number of hidden environmental covariates as the real-data contract (two), with species loadings of standard deviation 0.5.

1. Low rank. Species respond to 3 blends of 10 active axes chosen at random from the 64 columns. The embedding and the measured covariates each contribute a variance of about 1 on the logit scale to the linear predictor. Leftovers on the active axes carry about 10% of the embedding's share.
2. No blend structure. Each species responds to 3 axes of its own, with coefficients of standard deviation 0.5. This is the unfavourable case for the prior.
3. Null. The embedding has no effect.
4. Confounded management. Truth 1 plus a binary management covariate correlated at about 0.5 with a true blend (built from a score on the real embedding, dichotomised at the training median). In half the replicates its true effect is zero for every species. In the other half the community-average effect is 0.5 and it increases by 0.5 per standard deviation of the first trait. The harness builds this covariate, because the dataset has none. In truth 4 the embedding is a confounder, so the total and direct effects coincide there, and the truth tests recovery of the coefficient with the embedding held fixed. It does not test causal assumptions.

The strengths and the number of species are those fixed and recorded in the version 1 run. They are not tuned again. Version 1's pre-study tuned them so that the measured-only arm had a realistic held-out discrimination (an AUC of about 0.7 to 0.8). That AUC is computed again for every truth in this run and recorded, and the result says so if it has left that range. Taxonomy truths belong to part 2 of the package design and are not in version 2.

A cell is a truth, with truth 4 split into its two halves (true effect zero, true effect not zero). So there are five cells.

## Split and replicates

- Each replicate draws one blocked split on the real coordinates: blocks of the same side as the real-data contract, five folds, one fold held out at random with a recorded seed. For AlphaEco that fits 256 ponds and holds out 64. (The package spec asked for 250 held-out sites, which this dataset cannot supply alongside 250 fitting sites.)
- All arms are fitted to the same simulated dataset and the same split. Every comparison is a paired difference per replicate, grouped minus the other arm, summarised by its mean and 95% interval across replicates.
- A species seen fewer than 3 times in the fitting sites is dropped from every arm, since a fit cannot learn from none. The number dropped is recorded.

## Arms

- measured only: the fixed prior on the measured covariates (and the management covariate in truth 4), no embedding.
- raw: the fixed prior on the measured covariates and all embedding columns. It is fitted in truths 1 and 4 only, where the gradients criterion uses it.
- components: the fixed prior on the measured covariates plus the first 5 and, separately, the first 10 principal components of the embedding, computed on the training sites and projected onto the held-out sites.
- grouped: the grouped prior, covariate side, with the package defaults at the pinned commit, including the axis scale of 0.7 over K and the staged start (binary fits, from 5,000 burn-in).

Every arm uses 4 chains of 5,000 burn-in and 5,000 kept iterations, the length at which `B` reached nominal coverage in the October studies, and the same number of latent traits. Version 1 used 2 chains. The same fixed settings as the real-data contract apply otherwise.

## Criteria

Every table is reported overall and by band of true prevalence. The lowest band is species with true prevalence under 5%. A criterion must hold both overall and in the lowest band. The gradients and shared-part criteria have no bands, since their measures belong to sites and to fits. "Better" means the interval lies wholly on the better side of zero. "Not worse" means the interval excludes a loss larger than the tolerance. The tolerance is 1% of the measured-only arm's mean log loss per site-by-species cell for the log score, and 0.05 of R-squared for the gradients criterion.

- Prediction. The log score at the held-out sites. Grouped must be not worse than components in truths 1 to 4 and better in truth 1. In truth 3 it must be not worse than measured only.
- Rarity. The signed and absolute error of predicted prevalence at the held-out sites. In the lowest band the absolute error under grouped must be no more than 10% larger than under components, in truths 1 to 4.
- Management effect, in truth 4, at three levels: each species' coefficient, the community-average row of `G`, and the trait slope. The true community average is the mean of the species' true effects. For each level, coverage of the 95% interval must be at least 0.93 and no more than 0.02 below the components arm. Under a true effect of zero, species' intervals must exclude zero in no more than 7.5% of species-replicates. If the components arm itself falls below 0.93, that is reported as a finding about the fixed prior and grouped is judged on the relative rule alone.
- Gradients, in truths 1 and 4 (both halves of truth 4). The true blends must be recovered from the gradients of the fitted embedding effect, at the held-out sites. The measure is the adjusted blend R-squared defined below, with each arm at its own number of gradients. Grouped must be better than each components arm, and not worse than the raw arm, with a tolerance of 0.05 of R-squared. Truths 2 and 3 have no blends and are not scored. A fit is scored whatever its state.
- Shared part. For every cell the result gives, in its first lines on that cell, the number of grouped fits in the shared state, the mixed state (the chains of one fit in different states) and the leftover state, beside what the truth holds. A grouped fit is in the right state if every chain is in the state "shared structure" where the truth has shared gradients (truths 1 and 4), or every chain is in the state "leftovers dominate" where it has none (truth 2). Truth 3 has no effect and is not scored. If more than 10% of the grouped fits in a cell are not in the right state, the cell fails this criterion. The result says so, in the word "fails", in its first lines on that cell, and makes no claim about the shared part for that cell. A failure of this criterion does not stop the work (see "What happens on failure").

Also recorded: coverage and bias of the other measured coefficients, effective sample sizes, and wall time per fit against the fixed prior.

## Convergence

A species' response to the embedding is the sum of a shared part and a leftover. The data pin down the sum. They do not pin down how it is divided, which is the split. So a fit is judged on the targets, the quantities that a criterion judges or a user reads, and not on the scales that set the split. Every grouped fit gets two labels: whether it converged on the targets, and whether its split can be read.

- The targets. A fit fails if any of these holds. (a) R-hat is over 1.1 on a management coefficient at any of three levels: a species' coefficient, the community average, or the trait slope. The last two are the quantities whose coverage the management criterion judges. In the version 1 run the harness formed them per draw from the species' coefficients, by a ruling of that run, and the same is done here. (b) R-hat is over 1.1 on any species' coefficient for a measured covariate. (c) The chains' whole embedding effects at the fitted sites differ by more than 0.30 of the effect's size. (d) The chains' gradients agree less than 0.90. (e) A draw behind (a) to (d) is missing or not finite.
- For a grouped fit, (b), (c) and (d) are the rows `rhat_measured`, `disagreement` and `gradients` of `returnGroupedConvergence()`. For the other arms, (a) and (b) apply, from the package's convergence table, where the measured covariates are the pond variables and the management covariate.
- Small effects. When the embedding's effect is small the package does not assess (c) and (d). Such a fit is counted as "not assessed" on them, in a column of its own, and is judged on (a), (b) and (e). It is never counted as having passed (c) and (d). For every cell the result says in words how many fits this left unchecked: "the chains' agreement on the embedding was not checked in n of N fits", where N is the number of fits in the cell. For every fit the largest difference between two chains' whole embedding effects, on the scale of the linear predictor, is saved and reported, and decides nothing.
- The split. A grouped fit's split can be read only if every chain is in the state "shared structure" and R-hat on the logs of `sigma_emb` and `sigma_b` is at most 1.1 (the rows `state` and `rhat_scales` are both unflagged). Otherwise nothing is claimed from that fit about the shared part, the blends or the scales.
- A cell. If more than 10% of the grouped fits in a cell fail on the targets, the prior fails on mixing grounds in that cell, whatever the other criteria say.
- The wording for a cell. A cell that passes on the targets and fails the shared-part criterion is never reported as plainly "converged". It is reported as "converged on the targets; fails the shared-part criterion (n of N fits in the right state)".
- Fit sets. Every criterion is computed on all fits, on the fits converged under this rule, and on the fits converged under the version 1 rule. All three are reported. A pair of fits belongs to a converged set if its grouped fit converged, and the other arm's own failures are counted and reported beside the result. Where the other arm is the raw arm, the result says in words that the comparison is against fits that failed their own convergence check.
- The version 1 rule, for the third fit set and for the verdict reported beside this one: a fit fails if R-hat is over 1.1 on the log of `sigma_emb`, of `sigma_b` or of any axis scale, or on a species' management coefficient. If more than 10% of grouped fits fail in a cell, the prior fails on mixing grounds under version 1. This is the rule as the version 1 run computed it. Version 1's text also names a rank scale, and no rank scale exists at this commit.
- Reported, and deciding nothing. R-hat of every log axis scale (how many are over 1.1, and the largest). R-hat of log `sigma_emb` and of log `sigma_b`, each on its own. The share of embedding coefficients over 1.1. The row `disagreement_beyond`. Each chain's median `sigma_emb`. The size of the embedding's effect. Effective sample sizes of the targets. Each chain's held-out log loss. The R-hat of each species' embedding contribution. The version 1 verdict of every fit, with what failed it.
- The two rules side by side. For every cell the result gives the version 1 verdict first, in version 1's own words, and then the version 2 verdict, which is never given without it. Per cell, the four counts of version 1 pass or fail against version 2 pass or fail are reported, and every fit that fails version 1 and passes version 2 is listed by cell and replicate number, with what failed and its largest R-hat. The package's own verdict ("ok" or "review") is reported as a third line.

## The gradients measure

These definitions are the ones used in the two comparisons of 21 November 2026, which reproduced the package's own functions on the grouped arm.

- The effect. The arm's embedding columns as the fit has them (the 64 columns for grouped and raw, the 5 or 10 components for the components arms) times the posterior-mean coefficients, pooled over chains, centred at the fitted sites. It is a table of sites by species.
- The gradients. The leading directions of the effect across the fitted sites, from its singular value decomposition.
- Own number of gradients. The smallest number whose directions hold 90% of the effect at the fitted sites. This is the default of `returnEmbeddingGradients()`. A components arm cannot have more than its 5 or 10.
- Held-out scores. The gradients carried to the held-out sites through their weights on the embedding columns, as the `weights` of `returnEmbeddingGradients()` do. No held-out data enter the gradients.
- Blend R-squared. Each true blend's score at the held-out sites is regressed on the held-out scores of the k gradients, with an intercept. The value is the mean R-squared over the 3 blends.
- Adjusted. 1 minus (1 minus R-squared) times (m minus 1) over (m minus k minus 1), where m is the number of held-out sites and k the number of gradients.
- For the grouped arm the gradients and weights come from `returnEmbeddingGradients()`. That function refuses a fit without `embeddingVars`, so for the other arms the study's own code makes the same calculation. The study must show, on the grouped arm, that its code reproduces the package.
- Reported, and deciding nothing. Every arm at 3 gradients, and at each number from 1 to 10, with the unadjusted R-squared beside the adjusted one. Each arm's own number of gradients per fit. The reference with no fit: the true blends on the first 5 and the first 10 principal components of the embedding. The value at the fitted sites. The grouped arm by state. How much of each arm's effect its leading 3 gradients hold.

## Replicates and budget

- A first pass uses 30 replicates per cell, and the cells that decide a criterion are confirmed at 100.
- The first pass is 690 fits: 30 replicates in 5 cells for the grouped, measured-only and two components arms (600 fits), and 30 replicates in 3 cells for the raw arm (90 fits).
- No strengths are tuned. The dataset owner times one fit per arm with 4 chains, and projects the first pass from those times.
- A projection made by occJSDM from earlier runs, to be replaced by the dataset owner's own. With 2 chains a fit took 19 to 26 minutes (grouped), about 25 (raw), about 4 (10 components), about 3 (5 components) and under 4 (measured only). Doubled for 4 chains, the first pass is about 215 to 260 core-hours, which is 43 to 52 wall-hours on five workers. That sits at the line below.
- If the projection is over 48 wall-hours on five workers, the dataset owner stops and reports to Doug before launching. The usual cuts, in order, are fewer replicates in the first pass (10), truths 1 and 3 first, and dropping the 10-component arm. (Version 1's first cut, fewer species, is not open here, since the number of species is fixed.) A cut that removes an arm or a cell that a criterion needs leaves that criterion undecided by the run, and the result says so. The gradients criterion needs both components arms, the raw arm and truths 1 and 4. Any cut is recorded as a deviation.
- At most 5 fitting workers at once, and nothing is launched at full scale without the dataset owner's go-ahead. Every chain is single-threaded. A fit's chains may run one after another in one worker, or on separate cores (`ncores` in `MCMCparams`), provided no more than 5 cores are fitting at any time. On a fixed number of cores the wall time of the study is the same either way.

## Things to know before building

- The fits use the installed, pinned package, never `devtools::load_all()`. See the real-data contract's list of earlier mistakes, which all apply.
- The package's simulator `simulateOccJSDMData()` (`R/simulateData.R`) accepts a supplied embedding matrix and builds low-rank truths. BioDivAbove's `simulation/helpers.R` builds all four truths directly and is a working template.
- At the pinned commit the package has scale moves and a staged start for the grouped prior. Version 1's note about a tight start and axis-scale R-hat no longer applies.
- The harness of the version 1 run is the starting point. It needs new fields in its per-fit extract, which the output contract lists, and the gradients calculation, which exists in the dataset repository as two diagnostic scripts.
- With 4 chains the rows `disagreement` and `gradients` are the worst value over six pairs of chains. Their limits were set on fits with 4 chains.

## Optional module

- A two-stage set (about 100 sites, 30 species, 3 field samples and 3 PCRs) for the Gaoligongshan case, on truths 1 and 4 and the components and grouped arms. Not in version 2, because AlphaEco has no replicates.

## What happens on failure

- If part 1 fails on mixing grounds in a cell, or fails the prediction, rarity, management or gradients criterion, at 100 replicates, the dataset owner reports to Doug and nothing further is built on the prior.
- A failure of the shared-part criterion is reported to Doug in the same way and does not stop the work. What it does: every claim that the shared part finds the gradients is withdrawn for that cell, in the package's help and in any paper, and `returnEmbeddingBlends()` keeps its warning that it is to be relied on only where the split can be read.
- In `occupancy` and `two_stage` fits the community-average row and any taxon effects reach the data only through the latent occupancy states, which is what made `mu_b0` fail R-hat until the collapsed shift in `R/b0_shift.R` was added. If a failure traces to that, a new sampler move is not part of this study: report and stop.

## Disclosure

Written to be read by Alex. The full account is in two drafts beside this file: `simulation-grouped-prior-convergence-v2-draft.md` (section 6) and `simulation-grouped-prior-gradients-criterion-draft.md` (section 4). They also record each decision Doug made and what he was told before making it.

- Version 1 was written on 17 November 2026, before any fit. By its rule the prior failed on mixing grounds in every cell of the run of 20 November 2026 (150 grouped fits, 2 chains).
- Version 2 was written on 21 November 2026, after that run and after the diagnostics that followed it had been seen. Both versions, the package's convergence check and the gradient functions were written by Claude sessions working for Doug. Nobody independent of the results chose the quantities.
- The convergence rule stops judging the scales, which are the quantities that failed, and the two checks it adds on the embedding side were already known to flag none of the 150 fits. The reason given is that the data pin down the sum and not the split. The evidence for that reason was gathered after the failure.
- How the existing fits would score under the new rule was known exactly when the rule was decided: on the targets, 0, 5, 0, 9 and 4 of 30 fits would fail in the five cells, so two cells would pass and three would still fail. These counts are lower bounds, since the management check could be scored at the species level only. The failures are on the measured coefficients, which version 1 did not judge, and which the arms with few columns never fail.
- Two parts of the rule were decided without knowing how the existing fits score on them: the management check at the community-average and trait-slope levels, and the size of the difference between chains in fits whose embedding effect is small. How many fits the small-effect rule leaves unchecked was known: 0, 25, 30, 3 and 0 of 30 fits, in truths 1, 2 and 3 and the two halves of truth 4.
- The shared-part criterion was decided knowing that both halves of truth 4 would fail it on the existing fits (26 of 30 fits in the wrong state in each).
- The gradients criterion was written after two comparisons of the same kind had been seen on 10 replicates of two of its three cells. The form chosen would pass on those replicates. The lead over the raw arm in the zero-effect half of truth 4, where the signal is weak, was thin: a mean of +0.03, with an interval that reached just below zero.
- The limits of checks (c) and (d) have almost no support. 0.30 was set on 18 simulated fits whose answer was known, with next to no margin, and 0.90 is tested by no saved fit. Neither has ever flagged a fit of the validation study.
- The tolerance of 0.05 in the gradients criterion has no strong argument behind it.
- In the zero-effect half of truth 4 the grouped arm was not far above a reference that uses no species data: the true blends regressed on the first 10 principal components of the embedding gave 0.77, against the grouped arm's 0.81.
- For these reasons version 2 is a test only on fits made after it is issued, on new seeds.

## Versioning

A new version is issued by occJSDM. The dataset section writes the version it followed, and any deviation with its reason, in the run record.
