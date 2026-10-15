# Non-spatial candidates for B's undercoverage, 15 October 2026

The spatial 2x2 ([its report](../spatial-field-2x2-20261015/REPORT.md)) found that B's 95% intervals undercover most when there is no spatial component at all: 0.906 (SE 0.008) in `no_spatial`, and 0.685 for true effects over 3 in absolute value. The range test on branch `claude/range-inside-grid` also left B unchanged. This run tests the non-spatial candidates, each as an arm paired with `no_spatial` on the base cell's settings, with no spatial field simulated or fitted.

## Design

Six cells at 100 replicates each, all at seed label `base`:

-   `no_spatial`: the control, rerun so every arm comes from one revision. Its replicates should match the 2x2's `no_spatial` bit for bit.
-   `no_spatial_long`, Monte Carlo error: `no_spatial`'s datasets fitted with 5,000 burn-in and 5,000 sampling iterations per chain, against 1,000 and 1,000.
-   `no_factors`, latent factors competing with the covariates: no latent factors in the truth or the fit (d = 0).
-   `no_factors_fit2`: `no_factors`' datasets fitted with 2 latent factors. Against `no_factors` it isolates the cost of fitting factors that are not there, on identical data.
-   `no_spatial_notraits`, B's hierarchical prior: `no_spatial`'s datasets fitted without the traits.
-   `no_spatial_gt4`, a looser hierarchical prior on B: `no_spatial`'s datasets fitted with 4 latent traits instead of 2.

The brief asked for a looser prior on `sigma_b` if `listPriors` allows it. It does not: `runOccJSDM()` hard-codes an InvGamma(10, 1) prior on the variance `sigma_b^2`, and package code is not to be changed here. In the fit, B's prior is normal around `Tr G + A C` with standard deviation `sigma_b`, with N(0, 1) priors on the elements of G, A and C. So each latent trait adds 1 to every B element's prior variance, and `no_spatial_gt4` raises it from about 4.1 to about 6.1. That is the nearest looser prior available without changing package code. It loosens the prior's spread across species, not the `sigma_b` layer.

All runs use the study runner with 5 replicate processes and one thread per fit, 2 chains. The harness gains `fit_gt`, the number of latent traits fitted, and a per-cell `MCMCparams` that replaces the runner's settings for that cell only. The runner prints and records any per-cell settings and adds them to the checkpoint key only when a cell has them, so every earlier grid keeps its checkpoint name. Production cells leave both unset and fit as before.

## Pairing

Every arm draws `base`'s truths. The simulator draws B0, G, A, C and B before the factor scores and loadings, so B, B0 and G are identical in all six cells, `no_factors` included. `no_spatial_long`, `no_spatial_notraits` and `no_spatial_gt4` fit exactly `no_spatial`'s datasets, and `no_factors_fit2` exactly `no_factors`'. The two d = 0 arms have different datasets from the others, because the factors leave the linear predictor and the random draws after them shift. The harness tests in `tests/testthat/test-simstudy-b-candidates.R` check this at replicate 1.

## Run

Launched detached from this directory with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`. The log, provenance and results are ignored while the run is in flight so the runner records a clean tree.
