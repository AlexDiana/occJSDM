# Observed-data site WAIC: proposed correction for Alex

Status update, 7 October 2026: Alex approved PR #14 at `3eb9cd3` on 6 October. Conflict resolution against main `66ffa78` was verified; PR #14 merged on 7 October as `a56f554`. The September checks below remain historical evidence from the two-sample survey; former Lesson 3 model comparison is now in Lesson 4, whose current survey has three samples per site. Factor-count reliability validation remains open.

## What this fixes, in plain words

A model comparison should ask how well each candidate model explains the same observations. The previous scalar also scored the unobserved occupancy and collection states that each model had invented during fitting. Those hidden states are not extra observations, and different models can invent different states.

The new calculation asks: **how probable is this site's complete observed survey, allowing for all the hidden processes that could have produced it?** For a PCR dataset, that means its detections and nondetections, with the correct field samples and primers. A whole site's community is one scoring unit. This is the unit relevant to predicting an independently sampled site with the same covariates and sampling design.

Occupancy and collection states are summed out exactly. The unmeasured site factors are averaged over numerically. Species share those factors, so the calculation keeps their dependence: it combines species probabilities before averaging over the shared factors. Averaging each species separately and then multiplying would lose that information.

The WAIC formula applied afterwards is conventional. Its uncertainty and reliability checks are just as necessary after correcting its inputs.

## Review map

1. `R/site-waic.R`: exact observation/collection likelihood, posterior draw mapping, integration refinement, WAIC summaries and comparison guards. `src/site_waic.cpp`: small deterministic log-space integration kernel. The derivation is in `site-waic-plan.md`.
2. `R/output.R`: `extractWAIC()` now defaults to the observed-data score. This changes both its meaning and its cost: it calculates a score instead of instantly retrieving a scalar. `extractWAIC(fit, type = "legacy")` explicitly retrieves the old scalar. `results_output$WAIC` remains legacy for compatibility, and the existing teaching archive still stores those original values.
3. `R/runOccJSDM.R`: save the threshold used to turn read counts into detections. Older saved occupancy/two-stage fits must supply their actual fitting threshold, because that cannot be inferred reliably from the reads. New fits reject an inconsistent threshold at scoring time.
4. `tests/testthat/test-site-waic.R`: independent enumeration and integration references, numerical and mapping regressions. Public help and current Lesson 4 explain how to use the corrected score and why it does not yet settle the teaching example's factor count.

## Scope and limits

Supports non-spatial binary, occupancy and two-stage models, including environmental/collection covariates, traits already represented in stored species coefficients, multiple primers, unequal replication and missing PCRs. Hidden-state summaries are sufficient; full latent-state draws are not needed. Completely unobserved sites are omitted and named.

Spatial and continuous models are rejected by the new score. Their fitting and explicit legacy-score retrieval remain available. Spatial validation needs a target appropriate to correlated sites; this implementation must not be applied to them by pretending sites are independent.

Tensor quadrature grows quickly with factor count. Default settings can check up to three factors; some three-factor draws need a finer grid than the default node budget allows. It stops if the grid budget or refinement test fails. It does not silently return an unchecked score. Successive-grid agreement is a numerical check, not a rigorous bound, so a stricter independent grid sequence is part of the verification below.

Model comparisons require identical scored responses, response order, threshold, sample/primer mapping and site identities. Covariates and factor counts may differ. In particular, binary fits' synthetic `siteNames=1:n` must not hide a permutation of the actual sites in `data_info$Site`. Independent review found that case; a regression test failed before the identity correction and passes afterwards.

## Historical tests and independent numerical checks (22 September 2026)

- **Source suite:** 692 passing expectations, zero failures/errors/test warnings; the existing opt-in long coverage study is skipped. The new file contributes 51 expectations.
- **Mathematical references:** explicit enumeration of every occupancy and collection state in a small unbalanced, partly missing, multi-primer dataset; adaptive one-dimensional integration of shared-factor species likelihoods; known symmetric probabilities; analytic occupancy-only mixtures. Joint integration retains species dependence and agrees after an orthogonal factor rotation and equivalent SD/loading rescaling.
- **Numerical and API cases:** log likelihoods of -2,000 without underflow; impossible observations and nonfinite arrays rejected; posterior iterations/chains paired correctly; unchanged RNG state; old-fit threshold requirements; all-missing sites; incompatible model comparisons; three real short fitting paths with summarized latent states.
- **Reference software:** score, effective-parameter penalty and standard error agree with `loo::waic()` on both teaching fits' log-likelihood matrices. This cross-check validates the final WAIC arithmetic; the independent enumeration and integration checks validate its likelihood inputs.
- **Installed package:** `R CMD check`, with manual and vignettes disabled, completes with zero errors, three existing warnings and three notes. Installed tests and examples pass. Warnings concern the R header's unsupported clang warning option, the existing undocumented `predictNewSites(verbose)` argument and existing GNU Makefile extensions. Notes concern the linked worktree's `.git`, accepted LICENSE metadata and existing undefined globals. No new diagnostic names this implementation.
- **Documentation:** roxygen help generated, pkgdown reference index passes, and Lesson 3 renders to both HTML and Markdown. The old archive exporter/verifier explicitly read the legacy stored scalar and remain usable with the original package version. No saved fits or plots were regenerated.
- **Independent code review:** no remaining high/medium findings after the binary site-identity correction. Alex's scientific/code review was still required at this date; he approved the PR on 6 October.

## Historical saved teaching-fit check (22 September 2026)

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

Alex approved the target, likelihood and API change on 6 October 2026. For the current Lesson 4 factor-count example, score held-out *observed surveys* at independent sites, check MCMC precision using more retained draws, and assess whether WAIC differences agree with that validation. The existing lesson comparison against true occupancy probabilities remains useful but scores a different target. Do not mark factor-count selection complete or promise recovery of the generating count. Spatial model comparison remains separate work.


## Integration validation (7 October 2026)

Merged main `66ffa783991f8b19f9eb85e01ba499ae1db5508c` into the PR branch based on Alex-approved `3eb9cd33c2a7bee1980eb39ba8260e7fe86a3302`. All five conflicts are resolved. `R/site-waic.R` and `src/site_waic.cpp` are byte-for-byte unchanged from that approved head. The integration retains main's sampler, priors, threshold conversion and output fixes. Native exports were regenerated with Rcpp and the compiled registrations checked: both `sample_BBsL` routines have 17 arguments and `site_loglik_cpp` has 6.

- **Source tests:** 1,050 passing expectations, zero failures or test warnings, one explicitly opt-in coverage study skipped. The focused site-WAIC file passes 67 expectations, including 16 new assertions checking threshold-2 raw reads against pre-binarized threshold-1 observations in occupancy and two-stage fits. Both paths use the same seed and nondefault intercept-prior SD of 2; posterior results and scored likelihoods agree, raw reads and threshold metadata are retained, and explicit legacy extraction is unchanged.
- **Installed-package check:** clean integrated and clean recorded-main source archives were each built and checked with `--no-manual --no-build-vignettes` and vignette building disabled at the build step. Both report zero errors, five warnings and two notes. The shared warnings are the R header's unsupported clang warning option, undocumented `predictNewSites(verbose)`, GNU Makefile extensions, and two missing-vignette-output warnings caused by disabling vignette builds. Shared notes concern LICENSE metadata and existing R-code globals/analysis. The integrated installed suite passes 1,014 expectations, with zero failures/test warnings and eight standard skips (six CRAN-gated tests and two source-only checks); main passes 947 under the same settings. No new diagnostic category appears.
- **Help and lessons:** roxygen2 8.0.0 regenerates help without unrelated changes. Lessons 3 and 4 render in fresh sessions to Markdown and HTML; the changed sections were visually inspected. Numerical tables are unchanged apart from the explicit legacy-score label. Two regenerated Lesson 3 figures differed only in rendering details (legend order and a one-level colour difference); main's tracked figures were retained. All teaching bundles and publication settings match main. Lesson-link and archive-source-check regressions pass, as do all 10 lesson-site JavaScript tests.
- **Scope and provenance:** no archive was rescored, no fits were rerun for teaching, and no stored hashes or source allowlists were changed. Exporter changes preserve the legacy scalar but change source identity; the build README explains why historical archives require matching historical source and library. September's two-sample numerical validation is not evidence for the current three-sample fits.
- **Independent final review:** no material integration findings. The reviewer checked native interfaces, main-side preservation, threshold provenance, lesson placement, historical evidence and archive gates. Review did not certify current-fit scientific reliability, which remains open.
- **Diff checks:** no conflict entries remain. Unstaged changes and the staged PR diff against main pass whitespace checks. The full staged merge relative to the old PR head reports whitespace already present in main's archive/generated files; those unrelated files were preserved.

The conflict-resolution commit updates the existing PR only. It does not merge into main, publish lessons or complete factor-count selection. Alex's 6 October approval applies to the original head, not a fresh review of this integration.


## Current three-sample teaching survey after merge (7 October 2026)

The non-spatial teaching survey and all six fits were rebuilt from merged `a56f5548e68e09d244513d09fcc04ad52f0caa19` in a fresh archive, `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-post-pr14-20261007`, with its own installed library. The old archive remains intact. Every rebuilt posterior-results object is identical to its predecessor, with unchanged data, seeds, priors and MCMC schedules. New source hashes, retained thresholds and fit-file hashes supply current provenance. All eleven existing bundle exporters and their independent verifiers passed. The refreshed quickstart detection figure also reflects the accessible palette already present in merged code. This work changes no package estimator or fitting code and does not refresh the separate shipped `sampleresults` or historical spatial/four-package studies.

The current survey has 100 sites, 10 species, three field samples per site, two primers and six PCRs per primer/sample. The default two-factor fit and one-factor comparison use identical observations and 24,000 retained draws each. The following computational checks use evenly spaced subsets from all four chains. They are not a completed model-selection analysis.

| Quantity, 4,000 draws per fit | One factor | Two factors |
| --- | ---: | ---: |
| Observed-data site WAIC | 23,121.60 | 23,120.29 |
| Across-site SE of score | 375.22 | 375.31 |
| Effective-parameter penalty | 88.84 | 88.69 |
| Sites with pointwise penalty above 0.4 | 99/100 | 99/100 |
| Approximate Monte Carlo SE of score | 0.514 | 0.515 |
| Largest site log-likelihood change with stricter quadrature | 2.09e-6 | 1.97e-6 |

| Draws per fit | Difference, one minus two | Paired across-site SE | Approximate difference MCSE |
| --- | ---: | ---: | ---: |
| 1,000 | 0.4713 | 1.3201 | 1.3961 |
| 4,000 | 1.3114 | 0.7556 | 0.7274 |

The MCSE uses the same delta-method influence calculation as the historical check, preserving iteration order within chains. Difference MCSE assumes independent fitting runs. Both sample sizes agree with `loo::waic()` for score, effective-parameter penalty and across-site SE. Each used default quadrature and a separately started stricter sequence; the largest change in either total score over both checks is 1.91e-8. A larger posterior subset reduces estimated Monte Carlo uncertainty, but the score difference still changes appreciably and 99 of 100 sites remain flagged in each fit.

**Decision: this example does not justify choosing a factor count using WAIC.** Two factors have a small nominal advantage, but numerical correctness does not establish reliability of the WAIC approximation. No held-out observed-survey refits were run here. Their predictive assessment remains necessary before a positive model-selection claim. The independent-site occupancy exercise in Lesson 4 remains a different target. These results neither replace September's historical evidence nor complete the broader release gate.

Reproduce the scoring, export and independent compact-bundle verification from the repository root, with the new archive's library first in `R_LIBS`:

```sh
Rscript dev/simstudy/vignette-lesson/validate-waic.R ARCHIVE NEW_OUTPUT_1000 250
Rscript dev/simstudy/vignette-lesson/validate-waic.R ARCHIVE NEW_OUTPUT_4000 1000
Rscript dev/simstudy/vignette-lesson/export-waic.R ARCHIVE
Rscript dev/simstudy/vignette-lesson/verify-waic.R ARCHIVE
```

For the exporter, use `ARCHIVE/waic-1000` and `ARCHIVE/waic-4000` as the two fresh output directories. The validator rejects archive/library or fit/source mismatches and refuses to overwrite an output directory. The September validator is unchanged. The compact `vignettes/teaching-data/site-waic-lesson.rds` records both draw checks, pointwise results, source and input hashes, survey settings and the unresolved reliability conclusion. Machine-readable tables and fit comparisons are in `vignette-lesson/post-pr14-validation/`; full likelihood matrices remain in the archive. The new verifier independently recalculates the WAIC summaries and paired SE from those matrices and checks the compact bundle against its source files and actual fits.

The post-merge delivery checks also pass: lesson helper tests; all eleven existing bundle verifiers and the new WAIC verifier; Markdown and HTML renders of lessons 0-4 and the quickstart; lesson navigation and archive-source-check regressions; and all ten lesson-site JavaScript tests. An independent final reviewer checked all six fits, actual fit hashes, likelihood-matrix arithmetic, across-site SEs against `loo::waic()`, chain grouping, provenance and the statistical interpretation, with no actionable findings. The reviewer did not rerun fitting or quadrature and did not assess held-out predictive performance. Package estimator code and publication configuration are unchanged.
