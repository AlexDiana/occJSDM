# Lesson N: comparing four JSDMs when presence is observed perfectly

Design prepared 21 September 2026. Status: proposed pilot design; no comparison fits have been run.

## What we want to learn

Give occJSDM, gllvm, sjSDM and Hmsc the same simulated community, with a correct presence/absence record for every species at every sampled site. How accurately do they recover the underlying probabilities of occurrence? How well do they predict at new sites?

This removes field collection, PCR detection and false positives from the comparison. A species with a true occurrence probability of 20% still produces a 0 or a 1 when we visit a site. Even a perfectly observed survey does not tell us that probability directly. All four models must estimate it from the community data.

The first run is a small working example for students and a check that the four fitting and prediction workflows are comparable. It cannot establish which package generally performs best. The user has accepted this staged approach: begin with a shared-data pilot, then consider replicated communities and larger sample sizes. The non-spatial pilot does not depend on the spatial changes in PR #8.

## Which sjSDM version to use

Use Doug's **v0.2.1 release, explicitly running the PyTorch CPU backend**, for the main pilot. Freeze commit `d2ca508853a6e39df493f87c21d9c0136cfe652b`. This fork's release number is separate from its R package version: its `DESCRIPTION` says sjSDM 1.0.7. Record both.

There are three practical choices:

| Choice | Implication | Decision for this pilot |
|---|---|---|
| Fork v0.1.0, PyTorch | Contains the original Apple Silicon repairs but predates further fixes to Python dimensions passed from R. | Keep as a historical baseline, not the first choice. |
| Fork v0.2.1, PyTorch | Includes the later compatibility repairs and can explicitly bypass Mojo. | Main comparison. |
| Fork v0.2.1, Mojo | Accelerates the Monte Carlo likelihood; adds a compiler and worker binary whose versions must also be recorded and checked. | Separate optional backend comparison after the main pilot works. |

The v0.1.0 to v0.2.1 source diff shows integer conversion added to random-draw dimensions in both `model_sjSDM.py` and `dist_mvp.py`. The v0.2.1 fitter checks `SJSDM_MOJO_BACKEND`; the value `"0"` selects PyTorch. Set this before reticulate starts, in a fresh R process, and record the effective Python setting. Do not rely on automatic backend selection. The current checkout has the same R and Python package code as the v0.2.1 tag, but the experiment must use an immutable snapshot rather than depend on a moving checkout.

Use the existing working Python environment initially: Python 3.12.13 and PyTorch 2.5.1, reached through the repository's reticulate venv. Record its resolved executable and dependency versions at execution. Keep the R package snapshot and experiment outputs separate from the user's sjSDM checkout, which has an uncommitted tutorial edit.

Mojo **1.1.0 was released on 17 September 2026**, according to the [official releases page](https://mojolang.org/releases/). **Doug decided on 21 September not to upgrade.** Preserve the inspected local environment's Mojo 1.0.0 and MAX 26.5.0. The manifest allows version ranges; it does not itself guarantee those exact installed versions, so record the resolved versions and avoid an environment update. Any optional Mojo backend comparison will use this existing toolchain and the repository's numerical, memory and end-to-end parity checks. No compiler upgrade is planned.

## Pilot data and model scope

Use one newly simulated, non-spatial community with **100 fitting sites, 300 independent test sites and 10 species**. Use two measured environmental predictors and two normally distributed, unmeasured site factors. Species have different intercepts, environmental slopes and factor loadings, including positive and negative associations. Specify and save the complete parameter table and seed before any fit; do not select a community because one method performs particularly well or badly.

Generate probabilities with a logistic link and then sample perfectly observed binary occurrences. Retain both the probabilities and the actual 0/1 occurrences. Use an explicit, short simulator that is independent of all four fitting packages. This dataset is a new experiment, not a rerun of the older 11.8-point sample-size result or the existing teaching dataset.

All four packages receive exactly the same training occurrence matrix and environmental table, with explicit site and species identifiers. The test occurrences and hidden factors are never supplied to a fitter or used to tune it. Learn any centring or scaling from the training predictors only and apply that transformation to test predictors and response-curve grids. Verify each package's additional internal scaling and undo it when reconstructing predictions.

Use species-specific linear environmental responses, an intercept, and a non-spatial residual community component. Do not supply traits, phylogeny, spatial coordinates, nonlinear smooths or neural networks in this baseline. Species traits would give some models extra information and require a separate lesson.

| Package | Pilot configuration |
|---|---|
| occJSDM | Binary JSDM mode, two site factors, no observation stages or spatial field. Disable optional latent-trait structure for this baseline. Record any unavoidable coefficient hierarchy. |
| gllvm | Installed version 2.0.15, binomial logit model with two unconstrained latent variables. Use and record a supported approximation method; verify the selected method in the fitted object. |
| sjSDM | Fork v0.2.1, PyTorch CPU, linear environmental model, binomial logit, explicit `bioticStruct(df = 2)`. Record regularisation and all optimisation settings. |
| Hmsc | Binary probit model, one non-spatial site random level, two factors using fixed minimum/maximum factor counts. No trait or phylogenetic predictors. Freeze the installed version before fitting. |

This matches information and intended residual rank, not every statistical assumption. Priors, shrinkage and fitting methods differ. sjSDM exposes a covariance factor dimension even though its usual description contrasts its approach with latent-variable methods. Verify actual fitted dimensions in every package instead of assuming their defaults are equivalent. If a requested restriction is unsupported, report it and revise the configuration before scoring.

Hmsc uses a probit link for this binary model whereas the first simulator uses logit. Both convert a continuous ecological score to a probability, but their curves differ. State this asymmetry alongside the pilot results. Before any general performance ranking, repeat the experiment with a probit-generating community and matched generating effects on the probability scale. Raw slopes from logit and probit fits are not directly interchangeable.

## Two different probability questions

### 1. How well do we reconstruct the sampled sites?

At a sampled site, all species' observed presences can help infer its hidden conditions. Compare each fitted site's probability with the simulated probability that includes that site's actual hidden factors. This is recovery at sites used for fitting, not an independent prediction test.

For occJSDM and Hmsc, use posterior mean occurrence probabilities, averaging probabilities across draws. For gllvm and sjSDM, check whether native outputs average over uncertainty in the site's hidden factors or simply substitute a point estimate or zero. Where needed, compute the conditional probability mean using a documented study helper, integrating the fitted two-dimensional factor distribution against that site's observed community. Validate this integration against higher-accuracy numerical integration on small examples. For Bayesian fits, never multiply the data likelihood into joint posterior draws a second time.

Show a four-package recovery comparison only after this target has been verified for all four. If a comparable extraction is not achieved, explicitly mark that result unavailable; do not replace it with an environmental-only prediction under the same label. Keep native package outputs available separately when the study calculation differs from the ordinary prediction method.

### 2. How well do we predict a new, unsurveyed site?

At a new site we know the environment, but none of its species occurrences or hidden site factors. Predict by averaging over possible hidden conditions. Compare this with the true probability averaged over those same possible conditions. This is the primary common prediction comparison.

For example, a location's environment might imply an average 40% occurrence probability, while that particular simulated site's unmeasured conditions raise its probability to 70%. A model given only the environment cannot be expected to discover that particular 70%. The lesson must show both quantities and explain which one is being scored.

Setting a hidden factor to zero and then converting to a probability is generally different from averaging probabilities over hidden factors. The [gllvm prediction documentation](https://jenniniku.github.io/gllvm/reference/predict.gllvm.html) explicitly says `level = 0` sets latent variables to zero. The inspected sjSDM Python prediction path also applies the inverse link to the environmental score without integrating its fitted covariance. Therefore ordinary prediction calls cannot automatically be treated as the common marginal prediction.

Use analytical integration where available and checked numerical integration otherwise. For a probit predictor with normal hidden contribution of variance `v`, the marginal probability is `pnorm(mu / sqrt(1 + v))`; verify the package's definition of `v`. For a logit predictor, integrate `plogis(mu + sqrt(v) * u)` over standard-normal `u`. For Bayesian fits, do this within each parameter draw before averaging. Confirm prediction accuracy by tightening the numerical integration until changes are below 0.01 percentage points. This is a numerical tolerance, not an ecological release target.

Also score the held-out 0/1 occurrences with the Brier score and log score. These assess predictions of observable outcomes and are separate from error against known underlying probabilities. Test-site identifiers must be new to the fitted random-effect structure; test occurrences must never be supplied as conditional observations.

## What students will see

Write `vignettes/occJSDM-lesson-N.Rmd` with visible, readable tidyverse code and the existing teaching CSS. Show simulation, fitting and extraction code with explanatory variable names and blank lines between steps. Long fits will be displayed but not run while knitting; rendering will use a compact, verified results bundle.

1. **The same community given to all four models.** Show a small presence/absence matrix beside the known probabilities. Explain why a 0 is not a known probability of zero.
2. **True versus estimated probabilities.** Use identical 0% to 100% axes and a diagonal equality line, with separate panels for sampled-site recovery and new-site prediction. Label every panel's target and information available to the model.
3. **Direction and size of errors.** Report mean signed error and mean absolute error in percentage points, overall, by species, and in true-probability bands below 20%, 20% to 80%, and above 80%. Include cell counts. Explain with an explicitly hypothetical example why positive and negative errors can cancel in the signed mean while absolute error remains large. Any plotted numerical result must come from saved fits.
4. **Environmental responses.** Show fitted and true marginal occurrence curves across each environmental gradient, keeping the other measured predictor at its training mean. Compare curves on the probability scale rather than raw link-scale coefficients.
5. **What a small pilot can tell us.** Report fitting diagnostics, numerical prediction checks and runtime with the results. Distinguish one community's errors from evidence of systematic bias across repeated communities.

The first lesson will focus on probability recovery and prediction. Residual-correlation and variation-partitioning comparisons require agreed definitions and denominators across packages and belong in subsequent sublessons. Factor axes can rotate and change sign without changing the community model, so raw loading matrices are not a valid shared truth score.

## Fitting, evidence and checks

Run each package in its own R process. Use an isolated R library for packages that need installing; preserve the user's global library. The inventory currently has gllvm 2.0.15 and occJSDM installed, but not Hmsc or sjSDM as installed R packages. Use current occJSDM main revision `8654ff19cc67d6ac9e8f8bc477e5f2228f61dd98` as the intended source baseline, installing a matching snapshot rather than assuming the global installed copy matches it.

Archive inputs, truth, parameter tables, package/source hashes, seeds, complete fit settings, warnings, diagnostics, full fits and elapsed times. Save fitted and true probabilities with site/species identifiers in a compact bundle for the lesson. Rendering must neither refit models nor silently regenerate truth. Do not reuse the older occJSDM-only numbers as the new baseline.

Set fitting budgets and optimisation restart rules before seeing truth errors. Use multiple chains for the Bayesian models and inspect convergence of reported probabilities as well as identifiable parameter summaries. For optimised fits, inspect objective histories, convergence information and repeated starts. Retain failed or unstable attempts in the run log; do not select an optimiser start by which best matches the simulated truth. Equal iteration counts are not equal computational effort or equal convergence.

Before publishing results, verify: identical inputs across adapters; correct identifier order; training-only preprocessing; no test-data leakage; correct probability scales and shapes; comparable prediction targets; reproducible extraction after restarting R; and independent recomputation of every displayed error summary. Report unresolved fitting problems beside the relevant results. Neither converged fitting nor these checks proves that a model is unbiased.

## Subsequent stages

After the pilot and its target checks work, expand to replicated communities and nested sets of 100, 300 and 1,000 training sites. Keep generating ecological relationships fixed within each community, retain the original fitting sites for a paired recovery comparison, and keep an independent test set. Include both logit and probit generating scenarios before making broad comparative claims. Fix the replicate count and computing budget before launching that expansion.

An optional sjSDM backend check can then compare PyTorch and Mojo 1.0.0 on identical inputs, with more than one fitting seed and the actual backend recorded. Differences between optimisation runs must not automatically be attributed to the backend.

Spatial models, dispersal examples, observation-error comparisons, species-trait models and a comprehensive benchmark remain later work. No accuracy target or beta-release decision is being set by this pilot.
