# Targeted spatial recovery checks

## Question and scope

Doug approved extending the spatial validation on 27 September 2026 to feasible two-stage sampling, deliberately rare species, support-point sensitivity and current production code. This study retains the previous informative binary checks as historical references. It does not change package code, default priors or release targets. PRs #13 and #14 remain outside the study.

## Design fixed before fitting

Nine independent communities: three independently seeded datasets at each generating spatial range 0.1067, 0.1711 and 0.2356 (grid positions 4, 6 and 8, on the fitter's standardized coordinate scale). Each community has 100 uniformly distributed sites and eight species, with two species apiece averaging 1%, 5%, 25% and 75% occupancy probability over those sites. Species intercepts are solved to achieve these averages for the generated environment and spatial fields. This is conditional finite-site prevalence, not a claim about prevalence across a larger region. Zero-occupancy and zero-detection species are retained.

Independent unit-SD Gaussian spatial fields are generated directly from the full squared-exponential covariance. One standardized environmental predictor has alternating +0.4 and -0.4 slopes. No traits or residual factors are included; this isolates the spatial and detection questions and does not validate their interactions with traits or latent factors. There are no additional independent ecological observations sharing a location.

Two-stage surveys have two field samples per site, two primers, and six PCR replicates per primer: 200 field samples and 2,400 PCR rows. Collection intercept probabilities are drawn uniformly from 0.2 to 0.5, with alternating +1 and -1 standardized collection slopes. False-positive collection probabilities are uniform from 0.02 to 0.1. PCR true-detection probabilities are uniform from 0.3 to 0.6. Low-contamination q is uniform from 0.01 to 0.05; high-contamination q is uniform from 0.15 to 0.3. The high case is a deliberate difficult scenario in which some q and p values can be close.

Observations are direct Bernoulli detections, so p and q are positive-observation probabilities. They must not receive the lognormal positive-read correction used by the separate historical read-count study. Low and high contamination share biological truth, occupied/collected states and observation uniforms. The binary control receives those same true occupied/unoccupied states. This control isolates the additional difficulty from detection; it does not know the occupancy probabilities, coefficients or spatial field.

Every community is fitted under these three observation arms with 20, 50 and 100 spatial support points, for 81 fits. Twenty is the package default at 100 unique sites; 100 uses all observed locations. The support comparison is paired on data and fitting seeds. As selected support sets need not be nested, it compares the implemented settings, not a nested basis experiment. The generator uses the established simstudy_seed(label, replicate) rule with new spatial-specific labels.

## Numerical protocol

Freeze production revision d3d710e406c5b022df9cb84b7e29e783b4b225c4 and install it in a dedicated library. Record production, compiled-library, generator, runner and scoring hashes. At most four independent worker processes run, each with one sampler thread. Save each full fit before scoring and preserve inputs, RNG states, warnings and logs. Existing archives are read-only.

Seven short pilots check dimensions, truth alignment, reconstruction and runtime; pilot summaries are not recovery evidence. Initial fits use two chains, 3,000 burn-in and 5,000 retained iterations per chain. Select longer checks before comparing scientific outcomes when any scored group, species or element has Rhat above 1.05; a primary occupancy-group or species mean has ESS below 100; a non-range scored element has an unavailable Rhat; or the fitter reports a convergence warning. Longer checks use four chains, 6,000 burn-in and 12,000 retained iterations. Preserve both runs, report every unresolved flag, and assess whether substituting longer fits changes the community-level results. Constant range traces alone are reported and examined separately, not treated as proof of either failure or convergence. Do not continue escalating indefinitely if this finite protocol leaves poor mixing.

## Outcomes and interpretation

Primary outcomes are mean occupancy-probability error, MAE and RMSE at fitted sites, separately below 20%, from 20% through 80%, and above 80% truth. Whole-species prevalence groups are separate outcomes, including absolute and relative bias at 1% and 5%. Report occupied-site counts and detections. Compute each posterior probability by transforming every linear-predictor draw before averaging. Independently rebuild spatial bases, compare all ten with native package matrices, and compare selected posterior draws with the native linear predictor. Compare saved package probability means where available; binary fits do not store that array. A separate audit reconstructs all draws through the native basis projection and rechecks means, intervals and group scores.

Give communities equal weight. Treat the three ranges as fixed design strata: average their three-community means, estimate variance from variation among communities within each range, and use Satterthwaite degrees of freedom for approximate 95% t intervals. Report paired effects of support size and detection and the three-dataset summaries within each range. The nine-community average mixes three prespecified spatial ranges and is not nine replicates at each range. Keep signed errors alongside absolute errors. Report range and spatial-SD recovery, spatial-field attenuation, and the best least-squares representation of the true field at the true range for each selected support set.

The provisional spatial target remains average signed probability error within five percentage points in sufficiently informative cases, assessed separately by probability band. These operational scenarios test whether that target is achieved; they are not assumed informative enough in advance. No rare-species pass threshold is invented: a five-point error is substantial relative to 1% or 5% occupancy. Interval containment and widths are secondary descriptive checks, computed from the saved draws without additional refitting. Nine independent communities cannot establish nominal 95% coverage precisely. New-site prediction, categorical traits, factor interactions, larger sampling designs and a comprehensive coverage grid remain outside this targeted extension.

## Reproduction

Run from a checkout containing these research scripts, after archiving and installing the chosen production revision under STUDY/source-main and STUDY/library. STUDY/source-revision.txt records the exact revision. Use a fresh directory for a new study.

```sh
Rscript dev/simstudy/spatial-targeted-recheck/run.R --repo=. --study=STUDY --mode=prepare
Rscript dev/simstudy/spatial-targeted-recheck/run.R --repo=. --study=STUDY --mode=pilot --workers=4
Rscript dev/simstudy/spatial-targeted-recheck/run.R --repo=. --study=STUDY --mode=initial --workers=4
```

After the initial batch, use `summarise.R --mode=select` with the same repo/study arguments to write the diagnostic-selected keys. Run `run.R --mode=long --keys=KEYS` for that comma-separated list, then `summarise.R --mode=initial` and `--mode=final`. The summary outputs include per-community and aggregate substitution sensitivity. `verify.R --phase=selected --study=STUDY --selection=STUDY/summary-final/selected-fits.csv` audits every selected fit. `plot.R --summary=STUDY/summary-final` draws the figures.

The research regression tests run with `testthat::test_file("test-study.R", stop_on_failure=TRUE)` from this directory. Outputs are kept separately from tracked scripts and compact reports.
