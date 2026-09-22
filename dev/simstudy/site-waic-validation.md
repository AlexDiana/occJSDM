# Observed-data site WAIC: proposed correction for Alex

Status: **ALEX TO REVIEW**, 22 September 2026. Branch `codex/site-waic`, based on main `3a97267`. No changes to sampling, priors, fitted occupancy estimates or stored teaching fits.

## What this fixes, in plain words

A model comparison should ask how well each candidate model explains the same observations. The previous scalar also scored the unobserved occupancy and collection states that each model had invented during fitting. Those hidden states are not extra observations, and different models can invent different states.

The new calculation asks: **how probable is this site's complete observed survey, allowing for all the hidden processes that could have produced it?** For a PCR dataset, that means its detections and nondetections, with the correct field samples and primers. A whole site's community is one scoring unit. This is the unit relevant to predicting an independently sampled site with the same covariates and sampling design.

Occupancy and collection states are summed out exactly. The unmeasured site factors are averaged over numerically. Species share those factors, so the calculation keeps their dependence: it combines species probabilities before averaging over the shared factors. Averaging each species separately and then multiplying would lose that information.

The WAIC formula applied afterwards is conventional. Its uncertainty and reliability checks are just as necessary after correcting its inputs.

## Review map

1. `R/site-waic.R`: exact observation/collection likelihood, posterior draw mapping, integration refinement, WAIC summaries and comparison guards. `src/site_waic.cpp`: small deterministic log-space integration kernel. The derivation is in `site-waic-plan.md`.
2. `R/output.R`: `extractWAIC()` now defaults to the observed-data score. This changes both its meaning and its cost: it calculates a score instead of instantly retrieving a scalar. `extractWAIC(fit, type = "legacy")` explicitly retrieves the old scalar. `results_output$WAIC` remains legacy for compatibility, and the existing teaching archive still stores those original values.
3. `R/runOccJSDM.R`: save the threshold used to turn read counts into detections. Older saved occupancy/two-stage fits must supply their actual fitting threshold, because that cannot be inferred reliably from the reads. New fits reject an inconsistent threshold at scoring time.
4. `tests/testthat/test-site-waic.R`: independent enumeration and integration references, numerical and mapping regressions. Public help and Lesson 3 explain how to use the corrected score and why it does not yet settle the teaching example's factor count.

## Scope and limits

Supports non-spatial binary, occupancy and two-stage models, including environmental/collection covariates, traits already represented in stored species coefficients, multiple primers, unequal replication and missing PCRs. Hidden-state summaries are sufficient; full latent-state draws are not needed. Completely unobserved sites are omitted and named.

Spatial and continuous models are rejected by the new score. Their fitting and explicit legacy-score retrieval remain available. Spatial validation needs a target appropriate to correlated sites; this implementation must not be applied to them by pretending sites are independent.

Tensor quadrature grows quickly with factor count. Default settings can check up to three factors; some three-factor draws need a finer grid than the default node budget allows. It stops if the grid budget or refinement test fails. It does not silently return an unchecked score. Successive-grid agreement is a numerical check, not a rigorous bound, so a stricter independent grid sequence is part of the verification below.

Model comparisons require identical scored responses, response order, threshold, sample/primer mapping and site identities. Covariates and factor counts may differ. In particular, binary fits' synthetic `siteNames=1:n` must not hide a permutation of the actual sites in `data_info$Site`. Independent review found that case; a regression test failed before the identity correction and passes afterwards.

## Tests and independent numerical checks

- **Source suite:** 692 passing expectations, zero failures/errors/test warnings; the existing opt-in long coverage study is skipped. The new file contributes 51 expectations.
- **Mathematical references:** explicit enumeration of every occupancy and collection state in a small unbalanced, partly missing, multi-primer dataset; adaptive one-dimensional integration of shared-factor species likelihoods; known symmetric probabilities; analytic occupancy-only mixtures. Joint integration retains species dependence and agrees after an orthogonal factor rotation and equivalent SD/loading rescaling.
- **Numerical and API cases:** log likelihoods of -2,000 without underflow; impossible observations and nonfinite arrays rejected; posterior iterations/chains paired correctly; unchanged RNG state; old-fit threshold requirements; all-missing sites; incompatible model comparisons; three real short fitting paths with summarized latent states.
- **Reference software:** score, effective-parameter penalty and standard error agree with `loo::waic()` on both teaching fits' log-likelihood matrices. This cross-check validates the final WAIC arithmetic; the independent enumeration and integration checks validate its likelihood inputs.
- **Installed package:** `R CMD check`, with manual and vignettes disabled, completes with zero errors, three existing warnings and three notes. Installed tests and examples pass. Warnings concern the R header's unsupported clang warning option, the existing undocumented `predictNewSites(verbose)` argument and existing GNU Makefile extensions. Notes concern the linked worktree's `.git`, accepted LICENSE metadata and existing undefined globals. No new diagnostic names this implementation.
- **Documentation:** roxygen help generated, pkgdown reference index passes, and Lesson 3 renders to both HTML and Markdown. The old archive exporter/verifier explicitly read the legacy stored scalar and remain usable with the original package version. No saved fits or plots were regenerated.
- **Independent code review:** no remaining high/medium findings after the binary site-identity correction. Alex's scientific/code review is still required.

## Saved teaching-fit check

The existing one-factor and two-factor fits use the same observations: 100 sites, 10 species, two field samples per site, two primers and six PCR replicates per primer. The generating simulation has two factors. We did not refit either candidate.

This is a numerical validation using **1,000 of the 24,000 retained draws per fit**: 250 evenly spaced draws from each of four chains. It is not a final model-selection analysis.

| Quantity | One factor | Two factors |
|---|---:|---:|
| Observed-data site WAIC | 15,829.80 | 15,827.60 |
| Across-site standard error of the score | 292.92 | 292.53 |
| WAIC effective-parameter penalty | 86.03 | 84.54 |
| Sites with pointwise penalty above 0.4 | 98/100 | 98/100 |
| Approximate Monte Carlo SE of WAIC | 1.04 | 1.04 |
| Largest site log-likelihood change under stricter integration | 3.59e-7 | 1.08e-6 |

**Interpretation:** lower WAIC nominally favours two factors by 2.20 points. The standard error of that *paired difference* is 1.24 points. It is much smaller than either score's individual standard error because both models face the same easy and difficult sites. Monte Carlo uncertainty is separate: the approximate delta-method MCSE is about 1.04 per score, or about 1.47 for their difference if the two fitting runs are independent. Numerical integration changes either total score by less than 2e-8, so quadrature is not what limits this comparison.

Most importantly, **98 of 100 sites trigger the usual WAIC reliability diagnostic in each model**. A pointwise penalty above 0.4 is a warning about the approximation, not proof that this site's ecology is wrong or that the software failed. It tells us to check actual held-out-site predictions before trusting a ranking. The nominal preference happens to match the generating count of two, but this is not evidence that the selection procedure reliably recovers it.

A preliminary 200-draw check also passed numerical refinement, with a difference of 2.45 and paired SE of 2.32. Its changing estimates reinforce why a small posterior subset is appropriate for numerical validation, not a final factor-count decision.

## Reproduction

Run from this branch's repository root with R, `pkgload` and the unchanged teaching-fit archives. `loo` and `posterior`, when installed, provide the independent WAIC arithmetic check and approximate Monte Carlo SE. The script does not install packages.

```sh
Rscript dev/simstudy/validate_site_waic.R \
  /path/to/parent-of-teaching-fit-archives \
  /path/to/output-directory \
  250
```

The parent directory contains `prediction-lesson-20260922/one-factor-fit.rds` (MD5 `418ccf94af92351350602b3db54fc626`) and `vignette-lesson-20260919/default-fit.rds` (MD5 `67a364872efeb742de94e705194a7c33`). The script saves pointwise likelihoods, comparisons, input hashes and numerical diagnostics outside the package. It checks default quadrature orders 15/25/41/65 at tolerance 1e-4 against an independently started 31/51/81/121 sequence at tolerance 1e-6, on the same posterior draws. Both teaching fits required order 25 for 999 selected draws and order 41 for one draw under the default settings.

## What remains after this PR

Alex should review the target, likelihood and API change. For the Lesson 3 factor-count example, score held-out *observed surveys* at independent sites, check MCMC precision using more retained draws, and assess whether WAIC differences agree with that validation. The existing lesson comparison against true occupancy probabilities remains useful but scores a different target. Do not mark factor-count selection complete or promise recovery of the generating count. Spatial model comparison remains separate work.
