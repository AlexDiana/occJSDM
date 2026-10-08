# Targeted interval calibration follow up

Doug selected a targeted follow-up on 8 October 2026. This work corrects the historical q comparison, examines collection-coefficient intervals in the already selected nonspatial fits, and checks continuous-model B0 intervals with spatial effects absent. It does not change production code or priors and does not replace the paper's larger calibration study.

## Historical q

Preserve the 21 August 2026 four-cell archive and its intervals. Recompute q truth as nominal contamination-event probability times the probability that rounded lognormal reads are at least one: the normal upper-tail probability at log(1.5), with mean 1.5 and SD 1. Retain nominal and corrected coverage separately. Use all 30 communities in each low/high-q K3/K30 cell. K30 remains diagnostic, outside the practical limit of six PCR replicates per primer. This is corrected historical scoring, not a current-code MCMC rerun.

## Collection coefficients and current q

Reuse all 60 selected ordinary two-stage fits from the current-main recheck: ten communities under low/high contamination at K6 in each baseline, four-field-sample and 300-site design. Exclude the known-U conditional control. Preserve the original selected longer fits and their convergence flags. Transform true collection intercepts and slopes using each fit's actual covariate centers and scales, and verify the raw and fitted linear predictors agree. Report slopes separately for negative, zero and positive truth. Compute 95% equal-tailed intervals, width, signed error, coverage, directional misses, rank-normalized Rhat and bulk/tail ESS. Ten communities provide a targeted diagnostic, not precise evidence of universal calibration. These fits use production revision 2a75bf1, with later maintenance outside their provenance.

## Continuous B0

Freeze main c9ad954 in a private source/library archive, using an isolated R 4.5 launcher and its existing dependencies. Generate 100 independent communities using the historical continuous seed family and defaults: 100 sites, ten species, two environmental covariates, two measured and two unmeasured traits, noise SD one, and no spatial field. Each community supplies a paired d2 fit and a d0 control: preserve its intercepts, slopes, traits and Gaussian residuals, removing only the true hidden-factor contribution for d0. Fit both with default priors, two chains, 1,000 burn-in and 2,000 retained iterations per chain, with one sampler thread and at most four independent workers. Store data, truth, input hashes, RNG state, warnings, elapsed time, posterior summaries and fits. This changes both the generating and fitted factor structure; it is a mechanism check, not a same-data comparison.

Operational lesson added after execution: preserve the executed research scripts and keep them unchanged until every process exits. Fresh archives use atomic startup saves and per-shard status files. This does not change the prespecified generating scenarios, sampling schedules or selection rules.

Before the main run, verify both pilot models are inferred as continuous, have zero spatial support and no detection arrays, and reconstruct each generating mean on the fitted environmental scale. Any fitting warning, nonfinite Rhat, Rhat above 1.05 or bulk/tail ESS below 400 for B0 or tau triggers a longer fit using four chains, 3,000 burn-in and 6,000 retained iterations. Keep initial and selected results. Selection depends on diagnostics, never coverage. Report remaining flags and a sensitivity summary with flagged fits retained and excluded; do not remove them from the primary summary.

## Uncertainty and verification

Communities, not species or primers, are independent replicates. Calculate mean coverage and its Monte Carlo SE across per-community coverage fractions, with a 95% t interval clipped to [0,1]. Report point error and interval width together; negligible signed bias alone does not establish correct interval width. Independently reconstruct the corrected q decisions and the current collection truth from raw input values. Check all continuous coefficient/noise truth dimensions, hash identity and complete job counts. Preserve the old spatial continuous arm (200 communities, B0 coverage 0.879) as historical evidence; the clean configurations cannot establish a universal continuous-model result or explain that arm without further controlled work.
