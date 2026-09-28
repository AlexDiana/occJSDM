# Implementation and numerical validation

Production fitting code and the protocol were frozen in `4509629` before any scientific half-Cauchy fit. The report/analysis tooling was added in `cccca5d`. A later documentation clarification says that the intervention leaves all *other* occupancy/detection priors unchanged; the spatial-amplitude prior itself belongs to the occupancy model. This wording correction does not change the frozen installed fitting code.

## Package checks

The new sampler was tested against direct numerical integration of its conditional density at residual scales 0.0001, 0.3 and 2, including nonzero observed and latent trait means. Each posterior mean agrees within the prespecified 1.5% Monte Carlo tolerance. Further tests check trait-mean subtraction, unit rescaling, explicit prior activation, malformed values, nonspatial misuse and exact equivalence of implicit versus explicit default fits within one installation. All 24 new expectations pass.

The complete source test suite passed. `devtools::check(document = FALSE, cran = FALSE, args = '--no-manual')` completed with zero errors, three existing warnings and four notes. Installed tests passed 858 expectations, with zero failures/test warnings and the three expected source-only/opt-in skips. Examples and rebuilt vignettes passed. The warnings concern the existing compiler flags, undocumented `predictNewSites()` `verbose` argument and GNU extensions in `src/Makevars`. Logs are retained in the raw study archive under `package-tests.log`, `package-check.log` and `check/`.

## Reference compatibility

A saved binary input and RNG state were fitted with the frozen original and experimental installations in separate R sessions using the default prior. The maximum absolute difference across posterior outputs was 2.220446049250313e-15; warnings and final RNG state were exactly identical. This establishes numerical compatibility at the declared 1e-12 tolerance, not cross-installation bit identity. The saved inverse-gamma controls remain applicable.

## Research and independent review

All 25 field/initialization research expectations and 17 analysis expectations pass. They include randomized field reconstruction against a direct iteration-wise oracle, invariance to species offsets after centring, diagnostic failure cases, cloning the real fitting function while preserving every other expression, symmetric longer selection, gate completeness, and rejection of changed prior scales/source fingerprints. Independent review confirmed the auxiliary sampler algebra, matrix ordering, centring, update order and trait residual handling. It also confirmed the real initializer clone, matched budgets, complete gate inputs and atomic preparation of shared continuation settings.

The short half-Cauchy pipeline pilot took 14.9 seconds for 200 total MCMC iterations. Its independent native reconstruction reproduced field metrics within 2.05e-15, traces within 5.97e-15, probabilities/intervals within 8.05e-15 and group metrics within 2.78e-17. Diagnostics recomputed from the stored traces agree exactly. Folded-rank Rhat differs by up to 0.00345 between the two machine-precision reconstructions because of median-distance tie handling; no threshold crosses 1.05. These short pilot estimates are not scientific recovery evidence. The final audit must report the same checks and all threshold crossings for every selected scientific fit.

## Compute-only amendment

After the first four initial fits, Doug authorized using more than four CPUs. The continuation scheduler now caps the combined original initial pool and continuation workers at eight independent single-threaded fits on the 10-core, 32-GB machine. It queues longer pairs from diagnostic decisions as soon as each initial result becomes available. No sampling budget, prior, seed, diagnostic threshold or recovery criterion changed; the original fitting source and script fingerprints remain frozen.

## Outstanding scientific validation

The scientific initial fits, prescribed longer checks, initialization sensitivity, paired results and final independent audits must be complete before claiming a scientific result. The two-stage extension is conditional on the amended binary gate and satisfactory final numerical audit. No default change or merge follows automatically from testing the prior.

An operational logging fix restarted four newly dispatched continuation checks after about four minutes, before any of those four had saved a fit. R's `system2(stdout = file, stderr = TRUE)` captured output rather than returning the numeric exit status; the corrected adapter sends both streams to the same file. Four execution regressions verify stdout/stderr preservation, paths/arguments with spaces, and successful/failing exit codes. The original initial pool continued uninterrupted, all completed fits were preserved, and the restarted jobs used their original input RNG states and full prescribed budgets.

The already-required longer inverse-gamma control (`range6-rep03-binary-k100`) also passes the full native audit. Two historical script paths pointed into the worktree removed during the earlier user-requested cleanup; their current repository copies have exactly the original recorded hashes. The audit resolves only missing research-script paths this way and still requires matching content fingerprints.

## Posterior-moment finding and amended scoring

Before inspecting half-Cauchy recovery outcomes, a direct mathematical argument established that the full-rank binary likelihood retains a positive limit as the common spatial scale tends to infinity. Its half-Cauchy posterior tail is therefore proportional to inverse scale squared, leaving the amplitude posterior mean infinite and spatial-field means non-integrable. An independent reviewer confirmed the argument and the actual pilot's full-rank bases. The [dated amendment](ESTIMAND-AMENDMENT.md), committed at `162a577`, defines pointwise spatial medians for both priors, centring after medians, retained bounded occupancy means, and rank/quantile diagnostics. Legacy averages and the original gate remain archived but cannot support a finite posterior-mean comparison.

All 21 new robust-scoring expectations pass, including independently indexed medians/intervals, a skewed-draw distinction from means, explicit centring order, quantile ESS, missing-diagnostic rejection and pointwise diagnostic gate failures. A second independent randomized reconstruction with three species, three chains, 311 draws and multiple range blocks matched medians/endpoints/mean/RMS traces exactly and projection traces within 8.9e-16. The amended audit separately requires exact reproduction of all stored diagnostic summaries and numerical agreement with native spatial reconstruction; unresolved threshold crossings block the final extension decision.

The amended full-draw audit passed on the first original binary control: native pointwise field reconstruction differs by at most 1.41e-13, field medians/intervals/metrics by 8.97e-14, and probability summaries by 1.83e-14. Both trace and individual-field diagnostic reproductions agree exactly; no diagnostic threshold crosses. The final audit also recomputes the extension gate from the actual audited results, includes the separate initialization fit, invalidates stale decisions before reruns, and verifies the separately frozen scoring-dependency manifest. The complete research suite passes 67 expectations (25 original metrics, 17 original analysis, four execution and 21 amended metrics).
