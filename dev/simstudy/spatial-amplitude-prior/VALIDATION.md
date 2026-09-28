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

The scientific initial fits, prescribed longer checks, initialization sensitivity, paired results and final independent audits must be complete before claiming a scientific result. The two-stage extension is conditional on the prespecified binary gate and satisfactory final numerical audit. No default change or merge follows automatically from testing the prior.
