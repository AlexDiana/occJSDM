# Lesson 2 study-design sweep: when does site arrangement let occJSDM see a spatial field?

Protocol fixed on 1 October 2026, before any simulation or fit. Approved by Doug section by section in conversation on the same day; the decision log at the end records the choices made. Later changes go in dated amendments appended below, never by editing the frozen sections.

## Question and scope

The September spatial studies found that with 100 sites scattered at random over the study area and one occupancy state per site, the fitted spatial field was nearly flat and the range posterior uninformative, even when the spatial-amplitude prior was changed and even when an oracle sampler was handed every parameter except the field. The conditional diagnosis attributed this mainly to spacing: at the shortest generating range, 61 to 70% of sites had no other site correlated with them above 0.5. Doug's hypothesis is that the sites were spaced too far apart relative to the field's range.

This study tests that hypothesis under the design decision an ecologist actually controls. The survey budget is fixed at 100 sites in a fixed study area, and only the arrangement of those sites changes. It asks, for each arrangement, how much of the spatial field the data contain, how much of that occJSDM recovers from true occupancy states, how much is lost to eDNA detection error, and what the arrangement costs when predicting unsurveyed locations. The results become the worked part of Lesson 2.

The simulation is deliberately model-matched and controlled: every species shares one range and one amplitude, as the fitter assumes; the environment is a single broad gradient that every arrangement estimates equally well; nothing varies between cells except site placement. This isolates arrangement. Dispersal, species-specific ranges and same-scale environmental confounding are separate contrasts for later, and the lesson will say so. No package code, default or prior is changed by this study.

## Simulation, fixed before fitting

**Study area and coordinates.** The unit square. Raw coordinates lie in \[0, 1\] on each axis. The fitter standardises each axis by its own sites' mean and standard deviation; the generator records each arrangement's standardisation so that ranges and slopes can be compared on either scale.

**Environment.** One surface over the area: a zero-mean Gaussian random field with squared-exponential kernel, range 0.5 of the side and standard deviation 1, plus independent site noise with standard deviation 0.3. It is drawn once per community at every location used by any arrangement and at the test lattice. True slopes are defined on this raw surface scale. Each fit standardises the environment by its own survey sites, as a user would; the generator records each arrangement's mean and standard deviation so that fitted slopes convert to the raw scale for comparison with truth.

**Spatial field.** One independent Gaussian random field per species with squared-exponential kernel, correlation exp(-d^2 / (2 r^2)) at raw distance d, range r = 0.03 of the side and standard deviation 1 on the log-odds scale. Fields are drawn jointly at the union of all four arrangements' sites and the test lattice, so one realised landscape underlies every design within a community. The range is the shortest of the three used in September; on the fitter's standardised scale it is about 0.10 for a spread design, inside the grid of ten values from 0.01 to 0.30. Correlation falls below 0.5 at raw distance 0.035.

**Species.** Eight species with target region-wide occupancy 5, 5, 25, 25, 25, 75, 75 and 75%. Environmental slopes alternate +0.4 and -0.4 on the raw surface scale. Each intercept is solved by root finding so that the mean generating probability over the 1,600 test-lattice locations equals the target, with that species' field included. Prevalence is therefore a property of the landscape, the same for every arrangement, not of the surveyed sites. Zero-occupancy outcomes at the surveyed sites are retained if they occur.

**Test lattice.** A 40 by 40 lattice at the centres of equal cells, spacing 0.025, 1,600 locations. It receives true probabilities, environment and field values, and is never observed.

**Arrangements, 100 sites each.** Constructed from the community's data seed in this order, so that the spread design is the same set of points the pairs design extends.

1. Spread: 100 points uniform at random over the square.
2. Spread with close pairs: the first 80 points of the spread design; 20 of them chosen at random without replacement each receive a partner at raw distance 0.01 in a uniformly random direction, reflected back inside the square if necessary.
3. Clustered: ten cluster centres drawn uniformly at random subject to a minimum pairwise separation of 0.2, by rejection; each centre receives ten points uniform in a disc of radius 0.02 around it, clipped to the square. Within-cluster distances are mostly below 0.03.
4. Regular grid: the 10 by 10 lattice at coordinates 0.05, 0.15, ..., 0.95 on each axis, spacing 0.1. This is a control, not a realistic design.

For every arrangement the generator records: mean nearest-neighbour distance; the fraction of sites with at least one other site within raw distance 0.035, that is correlated above 0.5 under the true kernel; the effective rank tr(K)^2 / tr(K^2) of the true site covariance; the standardised range on each axis, which must lie within the fitter's grid of 0.01 to 0.30 or the design is flagged; and the range of environment values spanned by the sites as a fraction of the lattice range.

**Observations.** Two arms per community and arrangement.

- Binary control: the true occupied or unoccupied state at each surveyed site, one row per site.
- Two-stage survey: two field samples per site, two primers and six PCR replicates per primer, 200 field samples and 2,400 PCR rows, with direct Bernoulli detections exactly as in the September targeted recheck. Collection intercept probabilities are uniform on 0.2 to 0.5 with standardised collection slopes alternating +1 and -1 on a per-sample collection covariate drawn from a standard normal; false-positive collection probabilities uniform on 0.02 to 0.1; PCR true-detection probabilities uniform on 0.3 to 0.6; PCR false-positive probabilities uniform on 0.01 to 0.05, the low-contamination setting. Detection parameters are drawn once per community and shared by all arrangements; sample-level states and PCR outcomes are drawn per arrangement from that community's uniforms.

**Communities and seeds.** Three independent communities, replicates 1 to 3. Seeds follow the existing rule seed = sum of the label's byte values times 1,000 plus replicate, with labels spatial-design-data for the landscape, spatial-design-survey for sample and PCR draws, spatial-design-oracle for the oracle, and spatial-design-fit-ARM-ARRANGEMENT for each full fit, so every random stream is reproducible and distinct.

## The oracle, fixed before fitting

The elliptical-slice sampler from `dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp` is sourced with a recorded MD5 and not modified. For every community, arrangement and species, 96 posteriors, it receives the true occupied states at the surveyed sites, the true intercept, slope, range and amplitude, and the package-matched full-support covariance at the true range, and samples only the field. Schedule: four chains, 1,000 burn-in and 2,000 retained iterations, initialised at zero and at independent prior draws. Screens as in the diagnosis: rank-normalised Rhat at most 1.05 and bulk, tail and quantile ESS at least 100 for every field element and for the across-site average. A failing posterior is rerun once with both counts doubled; remaining flags are reported, never used to drop a cell. Oracle maps are pointwise posterior medians, centred across sites. The oracle is the perfect-observation ceiling for each cell; nothing in the full model can recover the field better from the same states.

## Full fits, fixed before fitting

All four arrangements receive full occJSDM fits in both arms for all three communities: 24 initial fits, prespecified so that no arrangement is chosen after seeing the oracle.

**Production code.** Frozen at the main revision current when the generator tests pass, installed into a dedicated library under the study directory with its hash recorded in `source-revision.txt` and in an amendment line appended below before the first fit. The fits use only that library.

**Settings.** `occCovariates = "environment"`, `spatCovariates = c("longitude", "latitude")`, `collCovariates = "collection"` for the two-stage arm, `listParams = list(n_factors = 0, n_lattrait = 0, n_supportpoints = 100)`. One support point per unique site is required: the default of 20 points chosen by clustering would place one point per cluster in the clustered design and discard exactly the within-cluster information that design exists to provide. All other arguments and every prior at their defaults.

**Schedule and selection.** Initial fits: two chains, 3,000 burn-in, 5,000 retained, thinning 1. A longer fit of four chains, 6,000 burn-in and 12,000 retained is run when any scored quantity has rank-normalised Rhat above 1.05, when a primary occupancy-group or species mean has bulk ESS below 100, when a scored non-range element has an unavailable Rhat, or when the fitter reports a convergence warning. Selection happens before any scientific outcome is compared. Both runs are kept; the longer one is reported where selected.

**Prediction.** For every fit, `predictNewSites()` at the 1,600 lattice locations with the lattice environment and coordinates, twice: with the spatial term, and with `useSpatial = FALSE` so that the field's contribution to prediction is visible. Latent factors are absent, so no averaging over hidden conditions is needed; the target at each lattice location is the true probability there, which includes the realised field.

**Execution.** Detached single-threaded worker processes launched with nohup, at most eight at once, each saving its full fit, warnings, RNG state and log before scoring. Existing archives are read-only. Raw outputs go to `dev/simstudy/results/spatial-design-20261001/`, which is gitignored and inside the repository tree so that the session can read it.

## Scoring, fixed before fitting

Every posterior probability is computed by transforming each linear-predictor draw and then averaging. Field draws at the surveyed sites are rebuilt from the saved spatial coefficients through the native basis, which is independently reconstructed from the kernel and checked against the package's matrices. Five outcomes per cell.

1. **Field recovery at surveyed sites.** RMSE and correlation of the centred posterior-median field against the centred true field, for the oracle and for both full arms, beside the zero-field baseline RMSE. The attenuation slope of median against truth is also recorded.
2. **Range and amplitude.** The posterior distribution of the range over the grid, with the probability within one grid step of the generating value and the probability at either boundary; posterior median and 95% interval of the spatial amplitude against 1.
3. **Occupancy at surveyed sites.** Signed error, mean absolute error and RMSE of posterior-mean probability against truth, separately below 20%, from 20% through 80% and above 80% truth, and separately for the 5%, 25% and 75% species groups. Report occupied-site counts and detections per group.
4. **Prediction at the lattice.** Signed error and mean absolute error of the posterior-mean predicted probability against the true probability at each of the 1,600 locations, overall and in four bins of distance to the nearest surveyed site: up to 0.02, 0.02 to 0.05, 0.05 to 0.1, and above 0.1. Reported with and without the spatial term, and for the oracle through its own field prediction at the lattice under the known kernel.
5. **Convergence.** Every Rhat, ESS and native warning for every scored quantity, reported beside the outcomes and never used to exclude a cell.

Communities receive equal weight. Each outcome is reported per community and as the mean of the three with their range; three communities cannot support confidence intervals and none are claimed. Paired contrasts are oracle against binary control, which isolates parameter estimation, and binary control against two-stage, which isolates detection, both within community and arrangement.

**Prespecified reading rules for the lesson.** An arrangement is described as informative for a species group when, in all three communities, the binary-control centred field RMSE is at least 20% below the zero-field baseline and the median-field correlation is at least 0.5. It is described as uninformative when the RMSE reduction is below 10% or the correlation below 0.3 in all three communities. Anything between is reported as intermediate. The range is described as recovered when at least half the posterior mass lies within one grid step of the generating value in all three communities. These rules fix the vocabulary of the lesson; the numbers themselves are reported regardless.

## Numerical validation

- `test-generator.R`: arrangements have exactly 100 sites inside the square, pairs are at distance 0.01, cluster centres are at least 0.2 apart, the grid is exact; fields reproduce from the seed; intercepts give the target lattice prevalence within 1e-6; survey rows map to the right sites, samples and primers; truth identities match.
- `test-oracle.R`: the sourced sampler's MD5 matches the diagnosis record; the covariance equals the package's full-support basis product to 1e-10; a known-answer check against logistic-normal quadrature for one site.
- `test-score.R`: scoring functions reproduce hand-computed values on a tiny synthetic fit; distance bins partition the lattice; the no-spatial prediction equals the environment-only calculation.
- `verify.R`: after fitting, independently rebuilds every selected fit's probabilities, field medians, range table and lattice predictions from the saved draws and compares them with the scored values to 1e-10; checks input, source and library hashes.

## Compute budget

Oracle: 96 posteriors at minutes each, under two hours in total. Full fits: 24 initial fits of roughly one to two hours each at 100 support points, eight at a time, about half a day of wall time; longer fits as selected, at most 24, about the same again. Prediction and scoring are minutes per fit. The study stops after the longer fits; it does not add arrangements, communities, ranges or priors without an amendment stating the question and budget.

## Lesson mapping

- 2A uses the generator's landscape maps and the arrangement table.
- 2B uses outcome 1 for the oracle and the oracle maps.
- 2C uses outcomes 1, 2, 3 and 5 for both full arms.
- 2D uses outcome 4.
- The compact bundle `vignettes/teaching-data/spatial-lesson.rds` carries truth, designs, oracle and fit summaries, lattice predictions and provenance hashes; `verify-teaching.R` for the lesson is written with the archive mode from the start.

## Reproduction

From the repository root, with `STUDY` an absolute path to a fresh directory under `dev/simstudy/results/`:

```sh
Rscript dev/simstudy/spatial-design-sweep/test-generator.R
Rscript dev/simstudy/spatial-design-sweep/test-oracle.R
Rscript dev/simstudy/spatial-design-sweep/test-score.R
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=STUDY --mode=prepare
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=STUDY --mode=oracle --workers=2
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=STUDY --mode=freeze
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=STUDY --mode=initial --workers=8
Rscript dev/simstudy/spatial-design-sweep/summarise.R --repo=. --study=STUDY --mode=select
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=STUDY --mode=long --workers=8
Rscript dev/simstudy/spatial-design-sweep/summarise.R --repo=. --study=STUDY --mode=final
Rscript dev/simstudy/spatial-design-sweep/verify.R --repo=. --study=STUDY
Rscript dev/simstudy/spatial-design-sweep/plot.R --repo=. --study=STUDY
Rscript dev/simstudy/spatial-design-sweep/export-teaching.R --repo=. --study=STUDY
```

Every script refuses to overwrite a completed output.

## Decision log

- **1 October 2026.** Doug chose to build Lesson 2 around a design sweep rather than the full 2A to 2D outline, with dispersal deferred. He chose site arrangement at a fixed budget of 100 sites as the primary axis, in place of site number at fixed extent, because effort is set by funding while placement is the ecologist's decision. He kept the regular grid as a control. He chose two-stage fits plus a binary control on the same communities, so that spacing and detection effects can be separated. He chose eight species with two at 5% and the rest common, dropping the 1% level already shown unrecoverable. He accepted a single broad environmental gradient after discussion of why a same-scale covariate would mix confounding into every cell; the confounding contrast is deferred to a separate bounded addition on the informative design. He approved the simulation, fitting and scoring, lesson structure and build-pipeline sections, and asked for this protocol to be written and committed before anything runs.

## Amendments

- **1 October 2026, amendment 1, before any simulation.** The fitter standardises each coordinate axis by its own standard deviation, so a field that is isotropic in raw coordinates is represented with a slightly anisotropic kernel whenever the two axes of an arrangement have unequal spreads. To keep this mismatch small and equal across designs, the arrangement generator rejects and redraws a spread design, together with the pairs design derived from it, or a set of cluster centres, whenever the 100 sites of the spread, pairs or clustered design have axis standard deviations differing by more than 10 percent, counting the rejections. The grid is exact. The design table records the ratio of axis standard deviations for every arrangement, and the oracle uses the true isotropic kernel at the sites, so the oracle-versus-binary-control contrast includes this representational difference alongside parameter estimation.
- The production revision used for the full fits will be appended here by the freeze step before the first fit.
- **01 October 2026, amendment 2.** Production code frozen at main revision `9af859795974839ee6038d00b6d8422ae77b047f` and installed into the study library before the first full fit.
- **1 October 2026, amendment 3, after the fits.** The final whole-branch review found implementation deviations from this protocol, recorded here after the fact; none changes a fit, a draw or a scoring function. (a) The longer-run selection check used the mean ESS where the protocol says bulk ESS. This selected one extra fit, `rep03-spread-two_stage` (occupancy species mean ESS 96.3, bulk ESS 109.9), for the longer run, a conservative difference. (b) The oracle's own lattice prediction (outcome 4 for the oracle) was stored for every cell but tabulated only after the fits, from the stored oracle results, as `results/oracle-lattice.csv`. (c) Three validation expectations that the protocol names for `test-score.R` (the distance bins partition the lattice, the no-spatial prediction equals `predictNewSites` with `useSpatial = FALSE`, and one hand-computed group MAE) were added after the fits, with no change to the scoring code. (d) The frozen revision `9af8597` in amendment 2 is a branch commit whose package tree is identical to main at `1526c26`, and amendment 2's "01 October" means 1 October. (e) Small implementation choices: the seeds use the labels `spatial-design-sites`, `-data`, `-detection`, `-survey-ARR`, `-oracle-ARR-SPECIES`, `-oracle-lattice-ARR-SPECIES` and `-fit-ARM-ARR`; a pair partner is placed by redrawing its angle until it falls inside the square rather than by reflection; the seed but not the RNG state is saved; the fits run as `mclapply` forks under one `nohup` parent; `verify.R` audits the results and the hashes of the result files but not the source or library hashes; and the lesson verifier is named `verify-lesson.R`.
