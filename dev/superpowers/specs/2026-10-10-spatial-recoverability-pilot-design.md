# Spatial recoverability pilot: common gradients and trait-linked dispersal fields

Date: 10 October 2026. Status: confirmed and in execution. Doug confirmed every open item with the review recommendations on 10 October and asked for the study to be run; the implementation plan is `dev/superpowers/plans/2026-10-10-spatial-recoverability-pilot.md`. A review earlier the same day found two mechanics that failed against the code this spec reuses and one outcome with no decision rule; the fixes are incorporated below and listed in "Revision record". The former open items are recorded as decisions in "Decisions recorded 10 October".

## Purpose

The spatial package comparison on branch `codex/spatial-comparison-plan` (report: `dev/simstudy/spatial-package-comparison/REPORT.md` on that branch) could not decide which spatial representation occJSDM should develop. Its design had no positive control: even an oracle that knew every nuisance parameter barely recovered held-out spatial structure (held-out contrast correlation 0.07, 0.20 and 0.21 for the 5%, 25% and 75% groups), so package differences could not be separated from a lack of information.

This study is that missing positive control, and nothing more. It asks, using oracle fits only: in which sampling regimes is spatial structure recoverable from single-season binary occupancy data, for two kinds of spatial truth?

1. Common gradients: a few shared latent fields, standing in for unmeasured environmental drivers that several species respond to.
2. Dispersal-like fields: an independent field per species whose spatial range is set by a dispersal trait, standing in for species-specific local dispersal.

The two truths differ in information content as well as in kind. Under truth 1, twelve species observe the same two fields, so each shared field is informed by roughly twelve times as many binary observations as any single species' field under truth 2. The expected result is therefore that truth 1 becomes recoverable at a smaller site count than truth 2, and every report from this study must say that this asymmetry follows from pooling, not from a difference in the fields themselves.

The eventual development question is which spatial representation to build into occJSDM. That question is deferred to a separate stage 2 specification, written only if this study finds recoverable regimes, and informed by its results. Doug's stated preference for the eventual model is a continuous trait controlling range; this study uses three dispersal classes only to establish whether the information exists.

## Non-goals

- No package fits (occJSDM, Hmsc, gllvm, sjSDM or others).
- No occJSDM prototype, including any trait-linked or per-species range sampler.
- No continuous-trait model.
- No mixed truth combining common gradients and dispersal fields in one community. Each truth is simulated separately.
- No dispersal dynamics (colonisation-extinction) generator and no multi-season data.
- No imperfect detection.
- No change to production code, priors, defaults or published lessons.

## Mechanism attribution is out of reach

A single-season snapshot cannot distinguish a species-specific spatial pattern caused by local dispersal from one caused by a species-specific response to an unmeasured local gradient. Both produce smooth, species-level patches. Attributing mechanism needs repeat seasons (colonisation depending on neighbours' earlier occupancy), a connectivity covariate, or dispersal traits that predict which species should be patchier. This study tests the statistical footprint that a snapshot can show, namely species-specific spatial structure at trait-linked scales, and every report from it must say that this footprint is not evidence of dispersal as a mechanism. The working assumption, stated by Doug, is that an ecologist can propose good candidate traits that correlate with dispersal; models in a later stage 2 would therefore receive the true trait.

## 1. Data-generating design (decided)

### Landscape and species

- Unit-square landscape. One measured environmental covariate, generated as in the spatial design sweep (`dev/simstudy/spatial-design-sweep/generator.R`): a smooth squared-exponential field with range 0.5 plus independent noise SD 0.3, with species slopes alternating +0.4 and -0.4.
- 12 species per community, three at each of four regional occupancy targets, 5, 25, 50 and 75%. Under truth 2 the three species at each occupancy are one from each dispersal class, so classes and occupancies are fully crossed. Under truth 1 there are no classes; the three species at each occupancy differ only in their loadings. Occupancy intercepts are solved so the mean occupancy over the held-out lattice equals the target, as in the sweep generator. (The four-level crossing was recommended in review; see open items.)
- 5% species are generated, fitted and reported, but are not used in the gate.
- No unstructured residual species association. Observations are true occupancy states (perfect detection), so spatial information is not confounded with detection error.

### Truth 1: common gradients

- Two shared latent fields, each a zero-mean squared-exponential Gaussian process with unit variance and range l in {0.03, 0.08, 0.15}. The range is a design factor (one value per cell).
- Each species s has loadings drawn independently from Normal(0, 1) on the two shared fields. Its spatial effect is `sqrt(0.9) * (lambda_s1 f1 + lambda_s2 f2) / norm(lambda_s) + sqrt(0.1) * g_s`, where g_s is a species-specific unit-variance field with the same range l. Every species therefore has total spatial SD 1 on the logit scale, with 90% of spatial variance shared and 10% independent. Because g_s shares the range of the shared fields, there is no scale separation between the two parts; the oracle separates them only through the known loadings. (The 10% independent part was recommended in review; see open items.)

### Truth 2: trait-linked dispersal fields

- Each species has a dispersal trait d_s, defined as the natural log of its spatial range, taking one of three values: log(0.03), log(0.08) and log(0.15), labelled short, medium and long. Four species per class.
- Each species has an independent unit-variance squared-exponential field with range `exp(d_s)`. Defining the trait as the log range makes range a genuine function of the trait, so a later study can draw d_s continuously and reuse the generator unchanged. (The three ranges are not log-linear in class labels 1, 2, 3, so a class-label parameterisation would have been a lookup table, not a function.)
- Spatial SD is 1 for every species.
- Every truth-2 community contains all three classes, so truth 2 has no range factor: recoverability is read out per class within the same communities.

### Half-correlation distances

For the squared-exponential kernel `exp(-d^2 / (2 * l^2))`, correlation falls to one half at d = l × sqrt(2 ln 2), about 1.18 l. The three ranges give half-correlation distances 0.035, 0.094 and 0.177. This is the distance convention used by the first comparison's report and distance audit, and all reports from this study use it.

### Sites and arrangements

- Training site counts n in {100, 300, 900}. 900 represents national-monitoring scale; if the short class fails at 900, single-season recovery of local dispersal at that scale is effectively out of reach, which is itself a reportable result.
- Two arrangements, both applied to the same community (the same latent fields and environment):
  - Random uniform: n points drawn uniformly in the unit square, with the sweep's rejection of draws whose x and y coordinate SDs differ too much. In the first study this arrangement was called "spread".
  - Clustered: k clusters of 10 sites, k = 10, 30 and 90 for n = 100, 300 and 900, each cluster drawn within the sweep's radius of 0.02 around its centre. The minimum centre separation scales with cluster count as `0.2 * sqrt(10 / k)`, giving 0.200, 0.115 and 0.067. The sweep's fixed separation of 0.2 cannot be kept: the review ran the sweep's own centre-placement function and found that the unit square holds at most about 20 centres at that separation, so 30 and 90 clusters fail to place. The scaled separations place successfully. Doubling effort therefore still means visiting more places rather than packing more sites into each, but the places are necessarily closer together; at k = 90 the separation is less than twice the cluster diameter, so neighbouring clusters nearly touch, and the report must note this when comparing clustered results across site counts.
- The gate (section 3) is judged on random uniform. Clustered is the realistic harder case and is reported alongside. It also shows whether close sites help short-range classes despite poorer coverage.
- Held-out recovery is scored on a regular 40 × 40 lattice of 1,600 points, disjoint from all training sites.

Expected sites within one half-correlation distance of a point, about n π r², for random uniform sites (coverage, about 1 - exp(-n π r²), in parentheses): short class 0.4, 1.2 and 3.5 sites at n = 100, 300, 900 (32%, 69%, 97%); medium class 2.8, 8.4 and 25; long class 9.8, 29 and 88. Edge effects make all of these somewhat lower. The n = 100 short-class coverage matches the 31 to 33% measured in the first study. These are geometric guides, not recovery predictions.

### Replication and cells

- 10 independent pilot communities per cell, each with its own seed. (10 was recommended in review; see open items.)
- Truth 1: 3 ranges × 3 site counts = 9 cells, 90 communities.
- Truth 2: 3 site counts = 3 cells, 30 communities.
- Each community is surveyed under both arrangements.
- Seeds use a new namespace. No community, seed or input from the first comparison or the sweep is reused, and no community is ever discarded or replaced after generation.

## 2. Oracles matched to each truth (proposed)

The oracle knows every generating parameter except the latent spatial fields: occupancy intercepts, environmental slopes, ranges, spatial SDs and, for truth 1, the loadings and the 90/10 variance split. It estimates only the latent fields from the binary training data and predicts held-out lattice values by Gaussian-process conditioning. It is therefore an upper reference for what a correctly specified model with perfectly known hyperparameters could recover, not a fitted competitor. Under truth 1 in particular, no fitted factor model will approach it, because a fitted model must also estimate the loadings and the number of factors; the report must label it as such.

- Truth 2 oracle: per species, the existing elliptical-slice sampler (`dev/simstudy/spatial-design-sweep/oracle.R`, which sources `dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp`) with that species' own range. The sampler itself needs no change. Its lattice prediction does: see "Numerical conditioning" below.
- Truth 1 oracle: a joint elliptical-slice sampler over the two shared fields and the twelve independent parts, with all twelve species' likelihoods. The prior is Gaussian with known block structure, so prior draws need only one factor of the n × n kernel per community, applied to each of the fourteen blocks. This is new code and must pass the verification in section 6 before any pilot fit. A per-species oracle would ignore the main source of information under this truth (twelve species observing the same patches) and would understate what is recoverable.
- Spatially uninformed reference: the best prediction without spatial information is the field-marginal probability `E[plogis(eta + f)]` with f drawn from the species' spatial prior, computed by Gauss-Hermite quadrature from the known intercept, slope and spatial SD. Using `plogis(eta)` instead would be biased and is not used.

### Numerical conditioning

The existing lattice predictor (`oracle_lattice`) forms the 1,600 × 1,600 conditional covariance of the lattice given the sites, adds jitter 1e-8 and takes a Cholesky factor. The review tested this on uniformly drawn sites in every cell of the grid. The Cholesky fails at range 0.15 for all three site counts and at range 0.08 for n = 300 and 900, five of the nine range × site-count combinations. The prior kernel itself remains factorisable but reaches a condition number near 1e10 at n = 900 and range 0.15, so the truth 1 joint sampler, which factorises the same kernel, inherits the hazard. The spec therefore requires:

- Lattice conditional covariance by symmetric eigendecomposition, keeping eigenvalues above a recorded relative tolerance of the largest (proposed 1e-10) and drawing from the truncated factor. The truncated covariance must reproduce the full conditional covariance to a recorded absolute tolerance, checked once per cell on the realised site layout. In the review's worst-cell test (n = 900, range 0.15) the truncation kept 546 of 1,600 eigenvalues and reproduced the conditional covariance to about 1e-6, with the conditional itself carrying a spurious negative eigenvalue near -4e-5 from forming it through an explicit inverse. The tolerance should therefore be of order 1e-6, not tighter, and the conditional should be formed with a triangular or eigen solve of the site kernel rather than an explicit inverse.
- The same eigendecomposition route for the prior kernel factor whenever the Cholesky with jitter 1e-8 fails or the condition number exceeds a recorded ceiling, with the chosen route and the condition number logged for every community and arrangement.
- No change to the Cholesky-based sampler path where it succeeds, so that the truth 2 oracle remains the sweep's sampler in the cells where the sweep's conditions hold.

Proposed MCMC budget, mirroring the first comparison's oracle references: four chains, 3,000 burn-in and 4,000 retained iterations, with one prescribed extension to 8,000 burn-in and 16,000 retained if the screen fails. Screen: rank-normalised Rhat at most 1.01 and bulk and tail ESS at least 400 for the field average, the field-truth projection and the probability average per species. The field-truth projection uses the true field, which is acceptable for an oracle that already knows every other parameter and follows the sweep's practice; the report states this. A cell whose oracle fails the screen after extension is "unresolved" and is neither pass nor fail.

## 3. Gate (structure discussed; thresholds proposed)

A cell passes only if both recoverability and headroom hold. Both are evaluated on held-out lattice points, random uniform arrangement, with the reference oracle matched to the truth. Under truth 1 they are evaluated per cell, pooling the 25%, 50% and 75% species. Under truth 2 they are evaluated per dispersal class within each site-count cell, pooling the same three occupancy levels.

### Recoverability: pattern and amplitude

- Pattern: correlation between the oracle's held-out spatial probability contrast (full prediction minus the spatially uninformed reference) and the true contrast. Proposed: median across the 10 communities at least 0.6, and at least 8 of 10 communities at least 0.4.
- Amplitude: attenuation slope of the oracle's held-out centred logit field regressed on the true centred field. Proposed: median between 0.5 and 1.5.

Both are required. A model that places patches correctly at a tenth of their true strength barely changes predictions, and a development decision needs both placement and strength. These thresholds are demanding against the first comparison's held-out attenuation slopes of about 0.04, which is appropriate for a positive control.

### Headroom: recovery must matter for prediction

- Held-out probability RMSE benefit of the oracle over the spatially uninformed reference, in percentage points. Proposed: mean across the 10 communities at least 1 percentage point, with a 95% interval across communities excluding zero.

Without headroom, every future package would land between the uninformed reference and the oracle in a band too narrow to separate them. One percentage point is a low floor relative to the ceiling that a species with logit-scale spatial SD 1 allows, so section 4 also reports each cell's benefit as a fraction of that ceiling, defined as the RMSE benefit of the true probability over the uninformed reference. The fraction is context for reading a bare pass; it is not part of the gate.

### Consistency across seeds

The median and count criteria above require recovery in most communities, not only on average. The full distribution across communities is reported for every cell, whether or not it passes.

### Outputs of the gate

- A recoverability map: pass, fail or unresolved for every truth × range or class × site count, with the clustered arrangement reported alongside each random-uniform result.
- A recommended positive-control setting and a stage 2 direction, read from the map by the rule below.
- Between-community SDs of every gate metric and of the oracle-minus-uninformed differences, to size a later package comparison (section 5).

### Reading the map: full, partial and no pass

The expected pattern (see Purpose) is that truth 1 passes at a smaller site count than truth 2. The rule therefore covers every combination, not only the case where both pass.

- Both pass: the recommended positive-control setting is the smallest site count at which truth 1 passes at range 0.08 and truth 2 passes for at least the medium class. A harder secondary setting is the next-smaller site count, or the short class at the recommended site count. Stage 2 would compare both shared-factor and per-species-range representations.
- Only truth 1 passes: the recommended setting is the smallest site count at which truth 1 passes at range 0.08. Stage 2 would be scoped to shared-factor representations, and the report states that species-specific ranges were not recoverable at any tested site count, or unresolved for any truth 2 cell whose screen failed (wording added 10 October, evening; see the revision record).
- Only truth 2 passes: the recommended setting is the smallest site count at which the medium class passes. Stage 2 would be scoped to per-species-range representations. This outcome is not expected and, if it occurs, the report must examine the truth 1 oracle for defects before any stage 2 is written.
- Neither passes: a stop result. If even n = 900 and range 0.15 fail under truth 1, that is reported as the finding and no stage 2 is proposed.

(This selection rule is a proposal.)

All thresholds are fixed and committed before any pilot community is generated, and are not changed after results are seen.

## 4. Scores and reported quantities (proposed)

For every species × community × arrangement × cell:

- Held-out contrast correlation and contrast RMS, against the true contrast RMS.
- Held-out centred logit field correlation, attenuation slope and RMSE, against a zero-field baseline.
- Held-out probability RMSE and MAE for the oracle, the uninformed reference and the true probability, so that the oracle's benefit can be expressed as a fraction of the attainable ceiling.
- Training-site versions of the same quantities, reported separately and never used in the gate.
- Per-species MCMC diagnostics and the screen status.
- Numerical conditioning records: kernel condition number, factorisation route and truncation tolerance used, per community and arrangement.
- Geometric descriptors: median nearest-training distance from lattice points and the fraction of lattice points within each half-correlation distance; for clustered layouts, the realised minimum centre separation.

Aggregation: average species within an occupancy level, then give communities equal weight. No species-level or community-level result is dropped; unresolved results are reported with their status.

## 5. Sizing a later package comparison (proposed)

This study does not run a package comparison. It supplies what is needed to size one. In a stage 2 comparison, arms would be compared within the same communities using paired differences d_i. With s the SD of d_i across communities and δ the smallest difference that would change a development decision, about 7.85 × (s / δ)² communities give 80% power for a two-sided 5% test, and about 15 × (s / δ)² give a 95% interval narrower than ±δ/2. Stage 1 cannot measure s between packages directly; the spread across communities of the oracle-minus-uninformed difference is reported as a conservative proxy. The meaningful difference δ is a development judgement for Doug and belongs in the stage 2 specification.

## 6. Verification before pilot fits (proposed)

- Generator: realised lattice occupancy equals each target; per-species spatial SD and the 90/10 shared split hold in expectation over seeds; empirical variograms reproduce each class's range; random-uniform and clustered layouts match the declared counts, cluster sizes and scaled separations, and centre placement succeeds for every community seed without exhausting the placement function's retry limit; the community seeds are distinct from earlier studies.
- Numerical conditioning: for every cell, on the realised site layout, the prior kernel's condition number is recorded; the lattice conditional covariance is formed by the eigendecomposition route and reproduces the full conditional covariance within the recorded tolerance; where the Cholesky route is used it agrees with the eigendecomposition route on prior draws with a fixed seed.
- Truth 2 oracle: the existing sampler's recorded deterministic checks still pass under the pinned runtime.
- Truth 1 joint oracle: reduces exactly to the existing single-field sampler when there is one species and one field with no independent part; recovers its prior when given no data; matches deterministic quadrature on a small logistic-normal case; reproduces fixed-seed draws.
- Uninformed reference: quadrature agrees with Monte Carlo integration to a stated tolerance.
- Scoring: synthetic cases with known correlation, slope and RMSE are reproduced exactly; an independent script recomputes every gate decision, including the map-reading rule in section 3, from saved per-community summaries.

## 7. Cost and execution (proposed)

Fit counts: truth 2 has 3 site counts × 10 communities × 2 arrangements × 12 species = 720 per-species oracle fits; truth 1 has 9 cells × 10 communities × 2 arrangements = 180 joint oracle fits, plus any extensions.

Before committing to this grid, time one 900-site oracle fit, per-species and joint, including lattice prediction, under the pinned runtime. Doug requested this timing on 10 October; it has not been done. The expected cost driver at 900 sites is the lattice prediction and the dense linear algebra described under numerical conditioning, not the slice sampling itself: the joint sampler draws fourteen fields but reuses one kernel factor, so its per-iteration cost is about fourteen single-species iterations. If the 900-site cost makes the full grid impractical on the laptop, reduce communities at n = 900 or move that level to another host, with the change recorded before any pilot fit.

Execution follows the existing conventions: the pinned R 4.5.0 launcher from the first comparison, single-threaded workers, raw draws on `/Volumes/occJSDM Research`, a run index under `dev/simstudy/results/`, and long jobs launched detached (`nohup` or `setsid`) so that an app restart does not kill them.

Provenance: an oracle-only pilot does not need the first comparison's full execution, recovery and audit infrastructure. It does need recorded seeds, hashes of every saved result, the conditioning log, and the independent gate recomputation listed in section 6. (Recommended in review; see open items.)

## Decisions recorded 10 October

Doug confirmed each former open item as recommended in review:

1. Four occupancy levels, 5, 25, 50 and 75%, crossed with the three classes under truth 2.
2. The 10% independent part in truth 1, sharing the shared fields' range.
3. Ten communities per cell. The 8-of-10 count criterion has coarse resolution at this size and is reported as such.
4. The truth 1 oracle receives the true loadings and is labelled an upper reference no fitted factor model will approach.
5. The gate thresholds and the map-reading rule in section 3, including the partial-pass branches, are fixed as written.
6. The lighter provenance set in section 7.
7. Eigendecomposition truncation at 1e-10 relative; reconstruction tolerance 1e-5 absolute; the prior factor switches from the Cholesky route to the eigen route when the Cholesky fails or the condition number exceeds 1e12.
8. The scaled cluster separation rule `0.2 * sqrt(10 / k)`.
9. The 900-site timing probe runs first. Fixed rule: the grid proceeds if the probe's upper-bound projection divided by eight workers is under 12 wall-clock hours; otherwise execution stops and Doug decides, and the grid is not reduced unilaterally.

Execution decisions taken when the run was set up, all within the spec's conventions:

- Runtime: the pinned private R 4.5.0 home built for the first comparison, reached through a study launcher that sets single-threaded numerics. System R is 4.6.1 and is not used. The run config records R and package versions and the Git revision.
- Storage: the raw per-iteration draws are not kept. Each job stores every tenth retained draw of the site fields (every fortieth under the extension, so that the stored count is the same), the lattice and site posterior means, the full convergence traces, screen results and conditioning records. Convergence diagnostics are computed inside the worker from the unthinned traces before thinning. Lattice prediction uses the stored thinned draws, about 1,600 per species.
- Posterior summary: the attenuation slope and field correlation use the posterior mean of the lattice field. The sweep used the median at sites; this study records the mean and says so.
- Job structure: one job is one community under one arrangement, so 60 truth 2 jobs of twelve per-species fits and 180 truth 1 joint fits, 240 jobs in all, run as forked single-threaded workers. An extension, if the screen fails, happens inside the job.
- Locations: run index and inputs at `dev/simstudy/results/spatial-recoverability-20261010/` in the ignored results directory; job files on `/Volumes/occJSDM Research/occJSDM/spatial-recoverability-20261010/jobs/` through a symlink; compact scores, the verification record and figures committed under `dev/simstudy/spatial-recoverability-pilot/results/`.
- Uninformed reference: 40-node Gauss-Hermite quadrature, checked against Monte Carlo integration.

## Revision record

10 October 2026, after review of the first draft:

- Clustered arrangement: replaced the fixed 0.2 centre separation, which cannot place 30 or 90 clusters, with a separation scaled by the square root of cluster count.
- Truth 2 oracle: withdrew the claim that no change is needed. The existing lattice predictor's Cholesky fails in five of nine range × site-count combinations; added the "Numerical conditioning" subsection and matching verification items.
- Truth 2 trait: defined d_s as the log range so that range is a function of the trait rather than a class lookup.
- Gate: added the expected truth 1 versus truth 2 asymmetry to Purpose, and replaced the single selection rule with a map-reading rule covering full, partial and no pass.
- Headroom: added the benefit as a fraction of the true-probability ceiling as a reported quantity, not a gate criterion.
- Section 2: labelled the truth 1 oracle as an upper reference no fitted factor model will approach, and noted that the convergence screen uses the true field.
- Section 7: recorded that lattice prediction and linear algebra, not slice sampling, are the expected cost driver, and set out the lighter provenance set.
- Wording: clarified that dispersal classes apply only to truth 2.

10 October 2026, evening, after the first run. Study directory `dev/simstudy/results/spatial-recoverability-20261010`; timing probe 476 s for the 900-site truth 2 job and 579 s for the 900-site truth 1 job, projected upper bound 36.9 CPU-hours, within the 12-hour rule; truth 2 fitted 15:13 to 16:22 BST, truth 1 16:23 to 20:19, scored 20:31, `verify.R` passing all 17 checks. All 180 truth 1 jobs failed the convergence screen after extension, so every truth 1 cell was unresolved and the map read stop. The diagnosis is `dev/simstudy/spatial-recoverability-pilot/DIAGNOSIS-truth1-mixing.md`, with the independent review and follow-up runs recorded in its last section. Doug decided the following after reading that review.

- Truth 1 sampler: the joint elliptical slice sampler (`joint-ellipse.cpp`, one ellipse angle for all fourteen blocks) is replaced by a blockwise one (`block-ellipse.cpp`, one angle per block) in commit ae13699. The prior, likelihood, loadings, variance split, budget in iterations, chains, thinning, extension rule, screen thresholds and seed labels are unchanged. This is a post hoc computational change made after results were seen, and it is recorded as such.
- What was known when the change was decided: the truth 1 screen outcome; the truth 1 per-species scores, which had been computed and written to `scores/species.csv` at 20:31 and which the review quoted (uniform arrangement, occupancy above 5 percent, range 0.08: median contrast correlation 0.61, 0.81 and 0.89 and median slope 0.40, 0.67 and 0.79 at n = 100, 300 and 900); and the finding that the first run's posterior means agree with the blockwise sampler's on the scored quantities to within 0.02 in slope and correlation on every job tested. The change is therefore expected to alter screen status, not scores, and a later pass of a truth 1 cell is not evidence that the sampler was tuned to pass.
- What differs in the re-run: the random-number stream, so draws are not comparable job by job; the per-iteration cost, about 1.36 times, so the compute budget rises although the iteration budget does not; and the source hashes, which supersede the `sources` field of `config.rds` for the files changed (block-ellipse.cpp 0f33cc9e766cfbd7820cea3d83900056, oracle.R 0045bd8d59651d371afbccbd6a1f58d6, run.R 46848c986819274d54a72e5c4aeda798, test-oracle.R 23140c3f4fbd5a0fae98334c9d375917, test-run.R 3870a6d02ce2ac8125fdaf4f4a255834; joint-ellipse.cpp 6e74865f7295a2f3ae675827e60c7756 is unchanged and still compiled and tested). The 180 first-run job files and their screen results are kept aside on the research volume and both result sets are reported. Truth 2 is not re-run.
- Expected ceiling, stated before the re-run: from the extension phase run on seven jobs with their own seeds (`diagnosis/extended-budget.log`), range 0.08 at n = 300 and range 0.15 at n = 900 are expected to resolve, range 0.08 at n = 900 and range 0.03 at n = 300 and above are not, and n = 100 is untested. Truth 2 at n = 900 is unresolved in all three classes under the frozen budget and the medium class never passes, so the map can reach at most the truth1_only branch. The report must say so and must not read the truth 2 n = 900 cells as failures.
- Section 3 wording: in the "Only truth 1 passes" branch the report says unresolved, not "not recoverable", for truth 2 cells whose screen failed. The decision rule itself is unchanged.
- Verification added to section 6 for the blockwise sampler: bit-for-bit reduction to the sweep's sampler for one species and one field, prior recovery with no data, agreement with three-dimensional quadrature on a one-site case with shared and independent blocks, rejection of invalid inputs, and a check that `run_joint_oracle` calls it. All six test files pass under the pinned runtime.

## Future work, out of scope here

- Stage 2 package comparison, including a continuous-trait occJSDM range prototype, written only if this study passes and informed by its results.
- A mixed truth combining shared gradients and trait-linked dispersal fields.
- A dynamic dispersal generator with multi-season data, which is needed for mechanism attribution.
- Imperfect detection and eDNA replicate structure.
