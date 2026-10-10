# Spatial package comparison

Completed scientific report, 10 October 2026. Unpublished teaching drafts are retained as records; Doug chose external-only evidence storage and cancelled repository bundle delivery. Independent numerical reconstruction passed all 288 jobs. The teaching export passed full archive checks and 449,110 independent arithmetic/accounting assertions, preserving the saved truth scores. Process completion does not establish diagnostic acceptance or ecological recovery.

## Experiment and provenance

The comparison uses three seeded binary communities, eight species, four arrangements of 100 training sites, and 1,600 held-out lattice sites per community. Regional occupancy targets are 5%, 25% and 75%. The generating spatial field has a squared-exponential covariance with length scale 0.03 in unit-square coordinates and logit-scale SD 1. Unstructured residual species association and imperfect detection are disabled. The two native spatial factors in Hmsc and gllvm remain, so those arms impose a different structure on eight independently generated species fields.

The original input archive was unavailable. The original generator, settings and seeds were used to regenerate the inputs, with 61 recorded checks against surviving evidence and a replay of occupancy uniforms. These are audited reconstructions of the same seeded communities, not restored historical files. All three complete landscapes and binary observations are retained. The original GP and oracle summaries remain separately labelled historical references; their schedules and field-summary convention differ from the new comparison.

All fitting used the frozen scientific snapshot and compatible R 4.5.0 launcher, isolated package library and pinned CPU Python environment. E11 authorised up to six single-threaded jobs after the four-worker pilot, retaining the conservative four-worker projection and 48-hour gate. E12 corrects only the administrative Hmsc loaded-version spelling check. Its eight pilot extensions were administrative and must not be counted as evidence of poor mixing. Original error/completion bytes and recovery receipts are preserved, with independent finish-time upper bounds for the four workers held during recovery, so idle hold time does not inflate the runtime projection. E14's bounded window, E15 resumption and E16 eight-worker handoff are preserved in the archive. E16 finished naturally with 354 completed tasks, zero interrupted or deferred tasks, and successful controller/supervisor exits. A fresh receipt confirmed all 357 recorded fitting/control PIDs and groups exited. A separate integrity audit verified every completed request identity and artifact checksum.

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

The verified scores and figures distinguish spatial patch placement from amplitude and compare held-out prediction with both the selected fit's spatial ablation and a separately fitted environment arm. They retain equal-community denominators and negative oracle gaps. Rare-species training attribution does not translate into a consistent held-out benefit: the attribution figure retains zero and negative benefits instead of interpreting a fitted spatial contribution as recovered spatial structure. Native hyperparameters retain their kernel and distance definitions; Moran graph distances and basis amplitudes are not fitted ranges.

The design did not include a positive-control scenario in which the oracle demonstrably recovered the target spatial structure at held-out sites. That is a serious limitation for the stated development question: deciding which representation is worth developing. A spatial SD of one in the generator establishes strong latent variation, but does not establish that sparse binary observations contain enough information to recover it. The historical sweep did find informative common-species field recovery at sampled clustered sites; that is not a demonstrated positive control for held-out recovery across the landscape. The historical oracle evidence should have prompted a scientific feasibility check before the full package comparison, alongside the runtime and integrity checks. Preserving the original scenarios took priority without adequately challenging whether they could discriminate useful methods.

Consequently, this experiment is chiefly a sparse, short-range stress test, rather than a decisive model-development comparison. It provides failure and diagnostic evidence but cannot cleanly separate insufficient recoverable information from inference and representation limitations. A replacement comparison should first demonstrate substantial oracle recovery in a prespecified positive-control regime, using the actual held-out probability and field targets. Only after that gate passes should the package arms be compared with representations capable of expressing the generating structure. Sparse, rare-species and structurally mismatched regimes can then serve as harder secondary cases. The gate should be calibrated on separate pilot simulations and fixed before evaluation, without selecting favourable seeds or changing criteria after package results are seen.

### Descriptive recovery from the saved scores

The following summaries passed the authoritative export's independent compact arithmetic checks. They average species within prevalence group, then give the three communities equal weight within each arrangement. The four arrangements have equal weight when all twelve community/arrangement cells are available. Unresolved and unstable selected numerical outputs remain descriptive results with their diagnostic labels; they are not promoted to accepted fits.

The oracle's held-out spatial probability contrast correlates with truth by 0.074, 0.200 and 0.207 for the 5%, 25% and 75% groups respectively. Its contrast RMS is 0.0183, 0.0510 and 0.0502, against true values 0.0531, 0.1722 and 0.1721. Matching-logit field attenuation slopes are 0.0114, 0.0419 and 0.0430. Thus even the reference knowing the nuisance parameters does not reliably recover patch placement or native field amplitude in this design. Probability contrast amplitude and logit-field attenuation measure different targets; their ratios must not be substituted for one another.

Held-out full-probability RMSE benefit relative to a separately fitted environment-only arm is shown below in probability percentage points. Positive means improvement. These small benefits do not establish spatial recovery: the occJSDM GP's mean improvements are about 0.043, 0.115 and 0.106 percentage points, and all twelve GP fits remain diagnostically unresolved. The Hmsc spatial-factor improvements are still smaller, with every fit unresolved. gllvm spatial factors slightly worsen mean prediction and include five unstable cells. Fixed trends and Moran bases worsen mean held-out prediction in all four packages under this short-range design.

| Package | Spatial arm | 5% benefit | 25% benefit | 75% benefit | Available / expected cells |
|---|---|---:|---:|---:|---:|
| occJSDM | GP | 0.043 | 0.115 | 0.106 | 12 / 12 |
| occJSDM | Quadratic trend | -0.298 | -0.948 | -0.626 | 12 / 12 |
| occJSDM | Moran basis | -0.438 | -1.100 | -1.228 | 12 / 12 |
| Hmsc | Two spatial factors | 0.003 | 0.029 | 0.007 | 12 / 12 |
| Hmsc | Quadratic trend | -1.567 | -2.164 | -1.888 | 12 / 12 |
| Hmsc | Moran basis | -2.618 | -3.143 | -3.400 | 12 / 12 |
| gllvm | Two spatial factors | -0.026 | -0.063 | -0.166 | 12 / 12 |
| gllvm | Quadratic trend | -8.317 | -3.543 | -3.154 | 12 / 12 |
| gllvm | Moran basis | -11.205 | -4.888 | -5.362 | 12 / 12 |
| sjSDM | Quadratic trend | -3.220 | -3.220 | -2.933 | 12 / 12 |
| sjSDM | Moran basis | -4.330 | -4.612 | -4.988 | 12 / 12 |
| sjSDM | Coordinate network | -3.035 | -7.058 | -3.645 | 1 / 12 |

The coordinate-network row describes its single unstable available cell, not a mean representative of all twelve datasets. Its eleven unavailable cells preclude ranking this arm. The environment-only comparison also differs from zeroing the selected fit's spatial term, because fitting a spatial term can change environmental coefficients. Both comparisons remain in the bundle and figures.

The development recommendation is to establish a discriminating positive-control design before choosing a replacement spatial representation. This run does not justify promoting the existing GP, replacing it with the fixed bases, or abandoning native spatial-factor or network models. A stable fit and a nonzero training attribution alone would not answer that question; held-out placement, strength, prediction benefit and uncertainty need to be recovered together.

## Deliverables and limits

### Native ranges and sampling distances

The archive-free reporting audit independently checks occJSDM's raw-axis half-correlation distances against the training-coordinate SDs and the native squared-exponential length scale. Different axis scales produce different distances along the two raw axes; a single isotropic raw length-scale recovery claim is therefore excluded. Hmsc's two exponential factor ranges are ordered within each draw, with zero denoting spatial independence. gllvm's exponential factor ranges are point estimates without invented intervals. These quantities are retained with diagnostic status and kernel definitions. The spatial coefficient prior SD is distinct from realised field amplitude, and Moran graph distances do not estimate a fitted kernel range. The audit passed all 193 checks on the final exported bundle; its receipt is `teaching-export-v3/native-distance-audit.rds` in the run archive.

The generating squared-exponential kernel reaches correlation one half at raw distance 0.0353. The held-out lattice's median nearest-training distance is 0.0395 for the grid, 0.0470--0.0489 for spread sites, 0.0514--0.0539 for close-pair designs and 0.1168--0.1297 for clustered sites. Only 25% of grid lattice locations, 31.0--32.6% for spread sites, 26.9--28.1% for pairs and 6.9--7.8% for clusters are within that half-correlation distance of a training site. These are geometric descriptions, not thresholds guaranteeing recovery. Close pairs supply short-distance information among sampled sites while leaving fewer distinct locations to cover the landscape.

The accepted export contains ordinary R objects and all 19 figures in six scientific figure sets. Its original bundle is 358,611,608 bytes, MD5 `9da703e5fa62b28442daaf06557edd69`. The successful export completed at 10:10:25 BST; a fresh receipt at 10:12 BST verified source/output/log hashes and the absence of all three owned process groups and descendants. The two earlier failed export stages are preserved with diagnoses of verifier defects; fitting, numerical verification and truth scoring were not repeated. Saved scores retain SHA256 `3b62b22cfe1e7c6d9775b62d5e8324f1fd98ecabdc65b661e7e3338d8715d6ef`.

Lesson 5's spatial section and reciprocal Lesson 7 links have rendered to HTML and unwrapped Markdown from the verified bundle. Every scientific figure was inspected for labels, clipping, units and missing-data treatment. Long fitting chunks stay disabled during knitting. Doug decided that this study should not be pursued for teaching publication and chose to keep the original bundle only on the external research drive (E17). Its copied bytes and sidecar were checked against the accepted export before replacing the laptop copies with ignored private links. No split, bundle commit or push is planned. Lessons remain unpublished records; no publication or production sampler/API change is part of this work. The whole-branch review found no new Critical or Important defect, and the final handover records the completed external-only disposition.

This experiment covers three communities, one generating spatial scale and a finite fitting budget. It cannot establish universal package superiority, universal unidentifiability of 5% species, an exact finite-sample oracle bound, or an ecological mechanism from spatial association. The original non-spatial experiment, this controlled spatial comparison and the historical detection/design study remain distinct.

Poor recovery here does not establish that published spatial submodels are generally underpowered. Statistical power was not estimated in this experiment. Hmsc and gllvm were constrained to two shared spatial factors, whereas truth contains eight independently generated species fields. Two shared patterns with species-specific loadings cannot generally represent eight independent maps. Their fitted exponential factor kernels also differ from the squared-exponential generating kernel. These restrictions describe the tested configurations, not the full flexibility of either package.

Information and inference must also be separated. Each species has one perfectly observed binary occurrence per training site; absence of detection error does not remove Bernoulli variation. At 5% regional occupancy, 100 sites provide roughly five presences on average. Short-range fields and gaps between sampled locations further limit information about unsampled patches. However, unresolved MCMC, unstable point solutions and failures to meet the numerical acceptance screen cannot be counted as reliable fits demonstrating a lack of spatial signal. Increasing site count alone would not remove the two-factor structural restriction.

The sjSDM coordinate network used a fixed two-layer 30-by-30 architecture, fixed regularization and diagonal residual covariance, without outcome-based tuning. The pinned fork's training and ordinary prediction targets differ, as documented above. Its acceptance failures concern this configuration and implementation. They are not a general assessment of published sjSDM performance or a reproduction of optimally tuned published applications.

A follow-up intended to distinguish information limits from inference should vary sampling density, shared versus independent generating fields and fitted factor dimension, while retaining a reliably converged reference. The present comparison can identify weaknesses under its declared stress test; it cannot justify abandoning a published spatial representation or selecting a universal winner.
