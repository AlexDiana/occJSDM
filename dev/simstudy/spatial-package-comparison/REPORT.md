# Spatial package comparison

Draft execution and diagnostic report, 10 October 2026. Independent numerical reconstruction passed all 288 jobs. Truth scores are saved; independent compact arithmetic, recovery conclusions and final teaching exports remain pending. Process completion does not establish diagnostic acceptance or ecological recovery.

## Experiment and provenance

The comparison uses three seeded binary communities, eight species, four arrangements of 100 training sites, and 1,600 held-out lattice sites per community. Regional occupancy targets are 5%, 25% and 75%. The generating spatial field has a squared-exponential covariance with length scale 0.03 in unit-square coordinates and logit-scale SD 1. Unstructured residual species association and imperfect detection are disabled. The two native spatial factors in Hmsc and gllvm remain, so those arms impose a different structure on eight independently generated species fields.

The original input archive was unavailable. The original generator, settings and seeds were used to regenerate the inputs, with 61 recorded checks against surviving evidence and a replay of occupancy uniforms. These are audited reconstructions of the same seeded communities, not restored historical files. All three complete landscapes and binary observations are retained. The original GP and oracle summaries remain separately labelled historical references; their schedules and field-summary convention differ from the new comparison.

All fitting used the frozen scientific snapshot and compatible R 4.5.0 launcher, isolated package library and pinned CPU Python environment. E12 corrects only the administrative Hmsc loaded-version spelling check. Its eight pilot extensions were administrative and must not be counted as evidence of poor mixing. E14's bounded window, E15 resumption and E16 eight-worker handoff are preserved in the archive. E16 finished naturally with 354 completed tasks, zero interrupted or deferred tasks, and successful controller/supervisor exits. A fresh receipt confirmed all 357 recorded fitting/control PIDs and groups exited. A separate integrity audit verified every completed request identity and artifact checksum.

## Finite budgets and diagnostic acceptance

Bayesian comparison fits use four chains, an initial schedule of 3,000 burn-in and 4,000 retained iterations per chain, and at most one prescribed extension of 8,000 burn-in and 16,000 retained iterations. The screen requires rank-normalised Rhat at most 1.01 and bulk and tail ESS at least 400 for the monitored parameters and training targets. Point fits use three starts, with three additional starts only when an acceptable best solution is not replicated. Numerical stationarity and objective/prediction agreement are checked without recovery truth. No diagnostic threshold or attempt ceiling was relaxed.

The immutable records account for 192 comparison cells and 96 species-specific oracle references. There were 462 comparison attempts, within the ceiling of 768, and 165 oracle attempts, within the ceiling of 192. No resource-interrupted attempt was retried. All 96 oracle references passed their diagnostic screen. Comparison statuses are:

| Package | Arm | Passed selection | Unresolved MCMC | Unstable point fit | No acceptable fit | Attempts |
|---|---|---:|---:|---:|---:|---:|
| occJSDM | Environment | 12 | 0 | 0 | 0 | 12 |
| occJSDM | Quadratic trend | 12 | 0 | 0 | 0 | 12 |
| occJSDM | Moran basis | 12 | 0 | 0 | 0 | 12 |
| occJSDM | GP | 0 | 12 | 0 | 0 | 24 |
| sjSDM | Environment | 12 | 0 | 0 | 0 | 36 |
| sjSDM | Quadratic trend | 12 | 0 | 0 | 0 | 39 |
| sjSDM | Moran basis | 12 | 0 | 0 | 0 | 39 |
| sjSDM | Coordinate network | 0 | 0 | 1 | 11 | 39 |
| Hmsc | Environment | 12 | 0 | 0 | 0 | 15 |
| Hmsc | Quadratic trend | 11 | 1 | 0 | 0 | 20 |
| Hmsc | Moran basis | 8 | 4 | 0 | 0 | 22 |
| Hmsc | Two spatial factors | 0 | 12 | 0 | 0 | 24 |
| gllvm | Environment | 12 | 0 | 0 | 0 | 36 |
| gllvm | Quadratic trend | 12 | 0 | 0 | 0 | 36 |
| gllvm | Moran basis | 12 | 0 | 0 | 0 | 36 |
| gllvm | Two spatial factors | 7 | 0 | 5 | 0 | 60 |
| **Comparison total** | | **146** | **29** | **6** | **11** | **462** |

An unresolved Bayesian or unstable point result retains its best prescribed numerical output for descriptive scoring, with that status attached. A cell with no acceptable solution has no fabricated prediction. In particular, the coordinate network supplies only one unstable cell; it cannot support an aggregate ranking against the other representations. All failed cells remain in denominators.

## Prediction targets and verification

The common target is marginal latent occurrence probability, with four versions of each selected fit: full, environment variation removed, spatial variation removed and both removed. Each draw's training mean spatial effect is moved into the baseline, preserving full probabilities. The same training centring is carried to the lattice. The spatial probability contrast is full minus spatially constant prediction. Means are obtained by averaging draw-level probabilities, not by applying an inverse link to a mean linear predictor.

Native logit-field recovery is reported only for matching logit units. Hmsc and the pinned sjSDM use probit-family likelihoods, so their native field amplitudes are not scored against logit truth. The pinned sjSDM fork's internal diagonal option removes off-diagonal residual association but retains independent normal integration in training. Its ordinary native prediction differs from that integrated marginal target; the study's saved-parameter reconstruction and certified quadrature preserve the fitted target. This finding is scoped to the pinned implementation.

Bayesian intervals use retained draws. occJSDM's subset-of-regressors prediction has coefficient uncertainty without added full-GP residual conditional variance. Hmsc uses certified conditional latent fields. gllvm's spatial intervals include fitted variational and conditional spatial uncertainty, holding coefficients, loadings and ranges fixed. sjSDM and gllvm fixed-basis arms have no invented spatial-probability intervals. Truth inclusion is descriptive for these landscapes, not a repeated-simulation coverage guarantee.

The archive-required numerical audit reconstructs all supported training summaries, complete point-fit lattice summaries and Bayesian predictions at the prespecified 50 lattice audit sites. It checks unthinned draw IDs, fixed RNG states, means, field medians, contrasts, intervals and Monte Carlo standard errors at absolute tolerance 1e-10. Truth is not opened in this audit. Scoring is gated on its successful completion and verified absence of all owned audit processes. The audit completed with 554 passed checks and 11 unavailable failed-fit cells, with maximum absolute discrepancy 1.9429e-16. Hash-bound completion records and a fresh process audit confirmed that all audit groups and descendants exited before truth scoring.

## Recovery, attribution and development recommendation

Pending verified truth scoring. The final report will distinguish spatial patch placement from amplitude, compare held-out prediction with both the selected fit's spatial ablation and a separately fitted environment arm, retain equal-community denominators and negative oracle gaps, and evaluate rare-species attribution against held-out benefit. Native hyperparameters will retain their kernel and distance definitions; Moran graph distances and basis amplitudes are not fitted ranges.

## Deliverables and limits

### Native ranges and sampling distances

The staged archive-free reporting audit independently checks occJSDM's raw-axis half-correlation distances against the training-coordinate SDs and the native squared-exponential length scale. Different axis scales produce different distances along the two raw axes; a single isotropic raw length-scale recovery claim is therefore excluded. Hmsc's two exponential factor ranges are ordered within each draw, with zero denoting spatial independence. gllvm's exponential factor ranges are point estimates without invented intervals. These quantities are retained with diagnostic status and kernel definitions. The spatial coefficient prior SD is distinct from realised field amplitude, and Moran graph distances do not estimate a fitted kernel range. The audit passed 193 checks on the diagnostic cache; final archive-bound export acceptance remains pending.

The generating squared-exponential kernel reaches correlation one half at raw distance 0.0353. The held-out lattice's median nearest-training distance is 0.0395 for the grid, 0.0470--0.0489 for spread sites, 0.0514--0.0539 for close-pair designs and 0.1168--0.1297 for clustered sites. Only 25% of grid lattice locations, 31.0--32.6% for spread sites, 26.9--28.1% for pairs and 6.9--7.8% for clusters are within that half-correlation distance of a training site. These are geometric descriptions, not thresholds guaranteeing recovery. Close pairs supply short-distance information among sampled sites while leaving fewer distinct locations to cover the landscape.

The completed deliverables will include an archive-free verified bundle, six sets of scientific figures, Lesson 5's spatial section and reciprocal Lesson 7 links, rendered HTML and unwrapped Markdown, and a final execution handover. Long fitting chunks stay disabled during knitting. Lessons remain drafts for Doug's line-by-line review; no publication or production sampler/API change is part of this work.

This experiment covers three communities, one generating spatial scale and a finite fitting budget. It cannot establish universal package superiority, universal unidentifiability of 5% species, an exact finite-sample oracle bound, or an ecological mechanism from spatial association. The original non-spatial experiment, this controlled spatial comparison and the historical detection/design study remain distinct.

Poor recovery here does not establish that published spatial submodels are generally underpowered. Statistical power was not estimated in this experiment. Hmsc and gllvm were constrained to two shared spatial factors, whereas truth contains eight independently generated species fields. Two shared patterns with species-specific loadings cannot generally represent eight independent maps. Their fitted exponential factor kernels also differ from the squared-exponential generating kernel. These restrictions describe the tested configurations, not the full flexibility of either package.

Information and inference must also be separated. Each species has one perfectly observed binary occurrence per training site; absence of detection error does not remove Bernoulli variation. At 5% regional occupancy, 100 sites provide roughly five presences on average. Short-range fields and gaps between sampled locations further limit information about unsampled patches. However, unresolved MCMC, unstable point solutions and failures to meet the numerical acceptance screen cannot be counted as reliable fits demonstrating a lack of spatial signal. Increasing site count alone would not remove the two-factor structural restriction.

The sjSDM coordinate network used a fixed two-layer 30-by-30 architecture, fixed regularization and diagonal residual covariance, without outcome-based tuning. The pinned fork's training and ordinary prediction targets differ, as documented above. Its acceptance failures concern this configuration and implementation. They are not a general assessment of published sjSDM performance or a reproduction of optimally tuned published applications.

A follow-up intended to distinguish information limits from inference should vary sampling density, shared versus independent generating fields and fitted factor dimension, while retaining a reliably converged reference. The present comparison can identify weaknesses under its declared stress test; it cannot justify abandoning a published spatial representation or selecting a universal winner.
