# Spatial fitting corrections and validation

**Status: DRAFT, ALEX TO REVIEW. The spatial point-bias gate remains open.** The comparison baseline is the RNG-safe branch, which also includes collection covariate alignment. No prior, feature switch, coordinate transformation or range grid is changed by this work.

## Demonstrated inconsistencies

The fitter stores spatial coefficients `Bs` with a range-independent normal prior. Its spatial field is `H_l %*% Bs`, where the basis depends on the current length scale. The old range update scores this reconstructed field under a full GP covariance at a proposed range, then keeps `Bs` and changes the basis. That changes the field itself, without scoring the resulting observation likelihood. A valid standalone range test with a fixed supplied GP field therefore does not validate this fitting transition. This finding concerns the representation passed by the full fitter; it does not reopen the historical amplitude or log-determinant experiments.

The previous fitter also kept only five columns of the whitened basis, selected by proximity to physical knots. A whitened column represents a combination of knots. Dropping such columns does not implement a local-knot approximation and makes the implied covariance depend on knot order. Prediction used all columns and a different numerical jitter. An independent deterministic nine-knot probe found a maximum fitted/predicted field difference of 0.499 with identical coefficients and identical jitter; reversing knot order changed fitted covariance entries by up to 0.450.

Finally, each joint coefficient sampler computed the spatial prior mean from observed and latent traits, but filled the spatial slots of its Gaussian prior mean vector with zero. Other spatial trait and variance updates assumed the computed mean. The R implementation, live serial C++ implementation and compatibility wrapper now use the same nonzero mean.

## Correct conditional and basis

Write `eta = Z_l b`, with `Z_l = [1, X, U, H_l]` and coefficient block `b = [B0, B, L, Bs]` for each species. Its prior mean `m` includes the current trait-predicted means of `B` and `Bs`. Its diagonal prior precision `D` has entries `1`, `1/sigma_b^2`, `1` and `1/sigma_bs^2` for those four blocks. Integrating this Gaussian coefficient block gives the range log weight

```text
Q_l,s = Z_l' diag(Omega_s) Z_l + D
h_l,s = Z_l' kappa_s + D m_s
log_weight(l) = sum_s(-logdet(Q_l,s)/2 + h_l,s' Q_l,s^-1 h_l,s/2)
                + log GammaDensity(l; a_l_s, b_l_s).
```

For continuous responses, `kappa = y * Omega` and `Omega = 1/tau^2`. For the binary occupancy states, `kappa = z - 1/2` and `Omega` is the Polya-Gamma precision. These are the Gaussian observation kernel and the [Polson, Scott and Windle augmented logistic kernel](https://users.wpi.edu/~balnan/PolsonScottWindle-JASA2013.pdf). Terms from the Gaussian prior normalizer and `m' D m` are constant across ranges and cancel. The code evaluates all ten grid values, draws from their categorical conditional, then immediately draws the entire integrated coefficient block at the selected range. Hyperparameters and factor scores are updated afterward. This ordering is essential: no update can condition on stale coefficients between range selection and their redraw. The neighboring-grid proposal and its endpoint asymmetry are removed.

An initial coherent update conditioned on `Bs` instead of integrating it out. Full fits exposed a practical mixing bottleneck: in one controlled continuous dataset, an independently integrated Gaussian calculation assigned about 40% probability to a neighboring grid value, while the chains visited it in only 2.5% to 5.3% of retained draws. The final blocked update addresses this evidence without changing the prior or stationary model.

Both fitting and prediction now use `H_l = K_nm L_mm^-T`, with `L_mm L_mm' = K_mm + 1e-5 I`. All selected support-point columns are retained. This gives the usual [subset-of-regressors covariance](https://gaussianprocess.org/gpml/chapters/RW8.pdf), `H_l H_l' = K_nm (K_mm + 1e-5 I)^-1 K_mn`, independently of knot order. The number of support points remains configurable. The approximation still needs enough support points to represent the spatial scale in the data.

The full GP inverse/Cholesky summaries are no longer computed or allocated during fitting. The internal `full_gp = TRUE` option retains them for standalone diagnostics. This avoids an unused `n_unique_sites` squared by ten allocation. Matrix dimensions are preserved for one support point, and prediction preserves the empty factor dimension for a spatial model without residual factors.

## Regression and full-fit evidence

The final focused spatial checks pass 50 assertions. All four fixes merge cleanly, and their combined source suite passes 681 assertions with only the opt-in coverage study skipped. New regressions fail against the preceding implementation and pass with the corrections. They check the integrated Gaussian likelihood against an independently derived observation-space covariance, the augmented binary likelihood against numerical quadrature, the exact grid prior without observations, and the actual range-then-coefficient ordering in a full fit. They also cover nonzero spatial prior means in all three coefficient samplers, fit/prediction basis agreement, analytic SoR covariance, knot-permutation invariance, repeated-coordinate row mapping, one-knot fitting/prediction and omission of unused dense matrices.

Full-basis arithmetic uses BLAS when every site has the canonical complete coefficient index layout; sparse and permuted historical layouts retain the original loops. Twenty-five independent algebra expectations cover both representations, precision matrices and conditional mean shifts. Benchmarks reduce 200 coefficient draws with 100 sites and 99 knots from 7.575 to 0.268 seconds, and 100 draws with 480 sites and 79 knots from 10.473 to 0.434 seconds. Floating-point summation order can change slightly; the conditional distribution is unchanged.

A paired full-fit probe using the current public simulator reproduces the reported upper-boundary behavior. It uses 120 sites, six species, 25 support points, nonzero environmental and spatial effects, raw spatial range 0.05 and continuous noise SD 0.2. Simulation seed 911 and fitting seed 912 are paired across versions; each fit uses two chains with 100 burn-in and 100 retained iterations. The two fitted coordinate axes imply ranges 0.1725 and 0.1689, inside the fitted 0.01 to 0.30 grid. Both old chains retain only the largest range: 200 of 200 draws. A state captured after actual coefficient updates shows that the old full-GP score increases to the upper boundary, while the correct conditional likelihood of those coefficients peaks at a different grid value.

The initial conditional-on-coefficients correction did not clear this weak-support-point probe: 150 of 200 retained draws still selected the upper boundary; field RMSE was 0.810 and correlation 0.687. This intermediate result motivated the additional mixing and information checks. Final blocked-update checks use exactly the standardized coordinates that the fitter receives, fixed support-point counts and distinct interior generating ranges. Their results and remaining limitations are recorded below.


## Remaining point-recovery limitations

The final blocked update also does not clear the original 25-knot probe: it retains the largest range in 156 of 200 draws, with field RMSE 0.809 and correlation 0.688. An oracle least-squares projection of the true field onto each of the ten fitted bases cannot achieve RMSE below 0.737 with these 25 knots (0.765 at the near-true range). This demonstrates inadequate support-point resolution independently of the sampler or prior, and the test remains an explicit failed recovery example.

Three controlled continuous datasets then use 100 sites, 12 independent spatial fields with SD 1, observation SD 0.1, 99 fixed support points, one environmental covariate, nonzero intercepts and no traits or residual factors. Their generating ranges are exactly grid indices 4, 6 and 8 on the standardized coordinates used by the fitter. They use two chains, 400 burn-in and 400 retained iterations per chain, and the existing priors. Each cell is paired between the initial conditional-on-coefficients implementation and the final blocked implementation.

- True range 0.1067: final mean range 0.1447, no upper-boundary draws, field RMSE 0.608, correlation 0.941 and estimated-on-true field slope 0.415.
- True range 0.1711: final mean range 0.2190, no upper-boundary draws, field RMSE 0.422, correlation 0.943 and field slope 0.654.
- True range 0.2356: final mean range 0.2973, 91.6% upper-boundary draws, field RMSE 0.371, correlation 0.943 and field slope 0.725.

The sampler now explores the intended posterior substantially better: retained range transitions across the two chains increase from 24, 4 and 0 to 232, 335 and 118. An independent observation-space Gaussian integration at saved variance draws gives neighboring-range probabilities 0.533/0.550 for the middle cell, versus empirical frequencies 0.508/0.520. Before blocking, the corresponding integrated probabilities were 0.404/0.432 but empirical frequencies only 0.025/0.053. This supports the mixing correction. The remaining attenuated field magnitudes are a separate limitation; correlations and pooled signed bias would conceal them.

The hardcoded continuous-noise prior is inverse-Gamma(5,5) on variance. At 100 sites it puts substantial pressure against generating variance 0.01, even with negligible residual error. Mean posterior noise SDs are 0.787, 0.619 and 0.552 in these three cells. An independent known-mean likelihood/prior diagnostic also shifts its preferred range upward. This evidence implicates prior sensitivity in the residual point bias; it does not establish that weakening this prior alone would solve every spatial case without transferring bias elsewhere. No prior was changed to force a recovery result.


## Informative continuous range check

A final controlled comparison uses 100 unique locations with ten distinct sites at each, 12 independent GP species fields, unit spatial and observation SDs, 99 fixed support points and the same nonzero environmental coefficients/intercepts. Locations, GP innovations, observation-noise draws and knot seeds are shared across three generating ranges on the fitted scale. Each final fit uses two chains, 200 burn-in and 200 retained iterations per chain, the unchanged priors and the cached blocked implementation. Full fits, truth and provenance were successfully saved using the validated portable runner.

- At true range 0.1067, both chains retain that range throughout. Field RMSE is 0.316, correlation 0.949 and field recovery slope 0.902; eta RMSE is 0.296. Mean spatial/noise SDs are 1.0045/1.0079.
- At true range 0.1711, both chains retain that range throughout. Field RMSE is 0.307, correlation 0.951 and field recovery slope 0.909; eta RMSE is 0.284. Mean spatial/noise SDs are 0.9980/1.0073.
- At true range 0.2356, both chains retain that range throughout. Field RMSE is 0.299, correlation 0.953 and field recovery slope 0.915; eta RMSE is 0.266. Mean spatial/noise SDs are 0.9972/1.0066.

An independent Gaussian calculation uses location means and an observation-space covariance to integrate the coefficients at saved variance draws. It assigns more than 99.9% probability to the true grid point in every cell/chain. The constant traces therefore reflect strong posterior concentration, unlike the earlier conditional-sampler bottleneck. Field mean errors are -0.0184, -0.0033 and +0.0050, while eta mean errors are about -0.011. These checks show distinct-range recovery and useful field pattern/level recovery when the data and support-point resolution are informative. They do not erase the low-noise, sparse-support or binary attenuation examples above and below, nor establish nominal coverage or absence of every finite-sample bias. The three cells share random inputs and are not independent simulation replicates.

## Binary range and probability check

A binary dataset uses 80 unique locations, six distinct sites at each location, eight independent GP species fields with SD 1, 79 support points, range 0.1711, modest nonzero intercepts/slopes, one environmental covariate and no traits or residual factors. The GP is generated on the exact coordinate scale subsequently fitted. The final blocked implementation uses the same complete dataset and seed as the intermediate conditional-on-coefficients run, with two chains and 100 burn-in plus 100 retained iterations per chain.

Both final range medians equal 0.1711, with means 0.1747 and 0.1792 and 95% grid intervals 0.1389 to 0.2033. Neither chain visits the upper boundary. Retained range transitions increase from 11/13 in the paired intermediate pilot to 50/57. Field RMSE is 0.660 and correlation 0.747. Occupancy probability RMSE is 0.138, and overall mean error is -0.0057. This pooled mean conceals attenuation: mean errors are +0.141 below true probability 0.2 and -0.126 above 0.8. Fitted species means closely track observed prevalence, within 0.0028, while these truth-stratum errors remain.

A longer intermediate fit with 500 burn-in and 500 retained draws per chain showed similar attenuation, so it cannot be dismissed as short-run fluctuation. The final range trace has R-hat 1.020, but other diagnostics warn about spatial coefficient convergence and spatial-scale ESS. These checks support the new range update and public binary fitting path; they do not establish general convergence or clear the spatial point-bias gate.

A planned 30-sites-per-location binary run was stopped after about ten minutes when a measured weighted-crossproduct benchmark projected roughly an hour for the full check. Inputs and provenance were retained, but it produced no fitted object or posterior results. It supplies no recovery evidence. Its portable generator remains available for a later adequately resourced run.

## Package and integration checks

The final combined branch includes collection alignment, RNG safety, factor correlation preservation and the spatial changes. It merges without source conflicts and passes 681 source assertions. Its installed-package check finishes with zero errors, three warnings and four notes; tests and examples pass. Warnings concern the unchanged clang/R-header warning flag, the undocumented `verbose` argument in `predictNewSites.Rd`, and GNU extensions in Makevars. The notes concern the worktree `.git` file, unavailable time verification, existing LICENSE metadata and existing undefined globals. GitHub reported no CI check runs on the new PRs at verification time. Testing was on arm64 macOS with R 4.5.0; cross-platform validation remains separate.

## Reproduction and release status

The tracked runners `validate_spatial_continuous.R` and `validate_spatial_binary.R` accept explicit precompiled source and fresh output directories. They save generated data and truth, all fitting settings and seeds, source/DLL hashes, warnings, fitted objects, point metrics, range distributions and session information. Their generated inputs were checked for exact equality with the original diagnostic scripts. The source path determines the implementation; changing a run label does not change model code.

From each source checkout, compile once with `Rscript -e 'pkgbuild::compile_dll(debug=FALSE, force=TRUE)'`. Then, from this branch:

```sh
Rscript dev/simstudy/validate_spatial_continuous.R --source=/path/to/compiled-source --out=/path/to/new-output --label=blocked --cell=2 --burn=400 --iter=400
Rscript dev/simstudy/validate_spatial_continuous.R --source=/path/to/compiled-source --out=/path/to/new-repeat-output --label=blocked-repeat10 --cell=2 --repeats=10 --tau=1 --shared-seed=true --family=repeat10 --burn=200 --iter=200
Rscript dev/simstudy/validate_spatial_binary.R --source=/path/to/compiled-source --out=/path/to/new-binary-output --label=blocked6 --grid-index=6 --repeats=6 --burn=100 --iter=100
```

Repeat continuous cells 1 to 3 for the distinct ranges. The original continuous cells use simulation seeds 202609111 to 202609113 and fitting seeds 202619111 to 202619113. The repeated-coordinate cells share seeds 202609111 and 202619111, and therefore provide a controlled range comparison rather than three independent replicates. The binary range-6 cell uses data seed 817006 and fitting seed 918006. Use the aligned/RNG-safe baseline source for an old-code comparison. That baseline cannot fit the repeated-coordinate continuous design: the old range score combines a 1,000-row field with a covariance for 100 unique coordinates and fails with non-conformable arguments. This is a failed baseline fit, not a posterior-recovery comparison.

Local detailed evidence is in the task's `work/spatial-probe`, `work/spatial-review-validation`, `work/spatial-binary-validation`, `work/spatial-dense-validation` and `work/spatial-runner-validation` directories. Only the compact report and reproduction scripts are tracked; generated fits and datasets stay outside the repository. The code baseline for the final blocked checks is `e9fd2d3`; `cba8a39` reuses Gaussian crossproducts without changing the model or priors. Snapshots record which variant each fit used.

The PR remains a draft because fixing the representation and sampler does not by itself clear the requested point-bias gate. Alex should review the algebra and update order, then resolve support-point adequacy and the demonstrated low-noise prior sensitivity before beta. Previously fitted spatial objects should be refitted. Full-basis and blocked updates can cost more for large spatial datasets despite the arithmetic improvements. These focused checks do not establish nominal interval coverage, overall convergence for every parameter, cross-platform reproducibility, or the separate collection/high-false-positive-rate release gates.
