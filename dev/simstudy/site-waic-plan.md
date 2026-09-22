# Proposed observed-data, site-level WAIC

Status: implementation and numerical checks complete, for Alex to review. This does not close the Lesson 3 factor-count selection exercise.

## Scientific target

Compare predictions for the complete observed community survey at a new, independent site, conditional on its measured covariates and sampling design. Support non-spatial binary, occupancy and two-stage models. Integrate occupancy and collection states exactly and integrate the shared Gaussian site factors jointly, after multiplying conditional species likelihoods. Never use fitted site scores, sampled occupancy states or sampled collection states as observations.

The existing stored scalar combines complete-data terms. Preserve it only through an explicitly named legacy extraction option. Make the corrected calculation the default of `extractWAIC()` and provide a detailed `computeSiteWAIC()` result for uncertainty and diagnostics. Unsupported spatial and continuous targets fail explicitly, without changing fitting or the legacy extraction option.

## Implementation checklist

- [x] Regression tests: distinguish observed from latent data; analytic one-factor integration; shared-species dependence; exact enumeration of two-stage states; missing and unbalanced observations; extreme probabilities; threshold provenance; matched posterior draws; paired comparison uncertainty.
- [x] Post-processing implementation in a separate R file, with one small deterministic C++ integration kernel. Stable log-space arithmetic throughout. Tensor Gaussian quadrature, refined until all site log likelihoods pass a declared tolerance; fail if the requested grid cannot pass. No sampler or prior changes.
- [x] Save the observation threshold in new fits. Require an explicit threshold for older replicated-observation fits and reject a conflicting threshold in new fits.
- [x] Return pointwise contributions, draws-by-site log likelihoods, quadrature diagnostics and high-variance WAIC warnings. Compare only identical scored observations and site grouping; give a paired standard error, without declaring a winner automatically.
- [x] Check saved one- and two-factor teaching fits, full tests, installed-package check and independent review. Record numerical evidence and remaining validation in a companion report.
- [x] Update public help, lesson caveats and TODO review queue.

Delivery: commit and publish as a separate PR for Alex. Leave merging and acceptance of the scientific target to review.

## Formula

For sample m and species j, let a be the PCR likelihood if DNA was collected and b the likelihood if it was not, multiplying over the PCRs and their actual primer identities. A missing response contributes a factor of one. Conditional on site occupancy z, the sample likelihood is theta(z) * a + (1 - theta(z)) * b. Multiply these sample likelihoods within a site to obtain g(0) and g(1). Conditional on the shared site factor vector u, the species likelihood is (1 - psi(u)) * g(0) + psi(u) * g(1). Multiply across species, then integrate u ~ Normal(0, sigma_h^2 I). Apply WAIC to this per-site observed-data log likelihood across matched posterior parameter draws.

For site i, lppd_i = log(mean(exp(log_lik[, i]))), p_i = var(log_lik[, i]), and WAIC_i = -2 * (lppd_i - p_i). Sum over sites. The paired standard error for a difference is sqrt(number_of_sites * var(WAIC_first_i - WAIC_second_i)). This is uncertainty across sites, not MCMC or quadrature uncertainty. A corrected likelihood does not make WAIC automatically reliable; high pointwise variance and convergence problems still require actual held-out-site checks.

## References

- [Merkle, Furr and Rabe-Hesketh: conditional versus marginal likelihoods for latent-variable model comparison](https://arxiv.org/abs/1802.04452).
- [Vehtari, Gelman and Gabry: practical Bayesian model evaluation using LOO and WAIC](https://arxiv.org/abs/1507.04544).
