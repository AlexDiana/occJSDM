# Diagnosis of poor spatial-field recovery

Requested by Doug on 28 September 2026 while the matched half-Cauchy longer fits were running. This is a separate conditional diagnosis; the frozen prior experiment and its selection rules remain unchanged.

## Questions and controlled comparisons

Reuse all nine saved binary communities and all eight species, including species with no occupied sites. First quantify the true kernel's neighbour correlations (excluding each site itself), effective covariance rank tr(K)^2/tr(K^2), and the proportion of true-field squared magnitude projected onto the intercept and environmental covariate. Report expected Bernoulli likelihood information, psi*(1-psi), by prevalence. These describe the actual geometry rather than assuming that 100 sites give many strongly correlated neighbours. Small ordinary projections do not rule out posterior confounding.

Use an independent elliptical slice sampler for a zero-mean Gaussian spatial field, fixing the spatial range to its generating value. Four conditional experiments change one assumption at a time relative to experiment B:

- A: supply the true intercept and environmental coefficient; fix spatial amplitude at the median of the current inverse-gamma variance prior, sqrt(1 / qgamma(0.5, shape = 10)).
- B: supply the true intercept and environmental coefficient; fix spatial amplitude at its generating value, 1.
- C: as B, but estimate the intercept under the unchanged Normal(0, 1) prior. Keep the environmental coefficient known. This assesses intercept estimation and its existing prior together, not an alternative intercept prior.
- D: as B, but supply ten independent ecological occupancy realisations per site: the original binary state plus nine reproducibly generated Bernoulli states with the same true occupancy probability. These are independent ecological states sharing a fixed spatial field, not field samples or PCR replicates of one occupancy state. This is an information diagnostic, not a proposed field design.

Use the package's full-support covariance at the true range, constructed independently from the kernel and whitening matrices. Check it against native projection and quantify its difference from the generating full GP covariance. This avoids silently changing the spatial model in the conditional comparison.

## Sampling, scoring and validation

Each species and configuration uses four chains with 1,000 burn-in and 2,000 retained iterations, initialized at zero and independent prior draws. Check every field entry, its across-site average/RMS/projection onto the true centred pattern, and the intercept where sampled using rank-normalized/folded Rhat, bulk ESS, median ESS and 2.5%/97.5% quantile ESS. Exclude fixed amplitude and range from diagnostics. Require Rhat <= 1.05 and every ESS >= 100; missing values fail. If any check fails, repeat that species/configuration once with four chains, 2,000 burn-in and 8,000 retained iterations. Preserve the initial summary and report any unresolved flags; do not extend indefinitely. Cap this separate diagnostic computation at two single-threaded workers while the eight original fitting workers continue.

Use pointwise posterior medians for spatial maps, then centre each median map across sites. Compare raw and centred RMSE, correlation, attenuation slope, interval containment and bounded occupancy-probability means. Include the deterministic zero-field RMSE baseline, whose correlation is undefined. Aggregate equally across communities, report prevalence-specific results and retain all species. The original prior-study comparison remains separate; these conditional experiments do not estimate an additive causal decomposition or a formal upper bound on recovery. Range-learning and environmental-coefficient uncertainty remain unisolated because both are held known.

Validate the independent sampler against deterministic univariate logistic-normal quadrature, a joint field/intercept case with analytically reducible Gaussian structure, and a constant-likelihood correlated Gaussian target. Check likelihoods against R's binomial log density, the matrix dimensions, exact reproducibility and source/input fingerprints. Retain source, seeds, settings, draws, diagnostics and numerical checks in a separate raw archive. Review conclusions against the actual conditional targets and convergence evidence before reporting causes.

The sampler follows Murray, Adams and MacKay (2010), [Elliptical slice sampling](https://proceedings.mlr.press/v9/murray10a.html). This investigation does not alter the production sampler, change a default prior, undertake the deferred full-model occupancy-intercept widening study, or automatically approve the two-stage half-Cauchy extension.
