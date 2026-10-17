# The prior on sigma_b, 17 October 2026

The two B reports ([non-spatial candidates](../b-nonspatial-candidates-20261015/REPORT.md), [4 latent traits with long chains](../b-gt4-long-20261016/REPORT.md)) traced B's 95% undercoverage in `no_spatial` (0.906 against 0.95) mainly to shrinkage under B's hierarchical prior, with unconverged chains as a smaller part, and reached 0.956 by fitting 4 latent traits instead of the simulated 2 at 5,000 + 5,000 iterations per chain. Both named one lever they could not test: the prior on `sigma_b`, the spread of each species' environmental responses around their trait prediction. `runOccJSDM()` hard-coded an InvGamma(10, 1) prior on `sigma_b^2`, whose shape counts like 20 pseudo-observations, as many as the 20 B elements in this design (10 species, 2 covariates), and `sigma_b`'s posterior sat at the prior's 0.33 in every arm, against a simulated truth of 0.35. This run adds a `listPriors` option for that prior and tests whether loosening it restores B's coverage, on the same datasets as before.

## The option

`listPriors$sigma_b_prior` is `"inverse_gamma"` (the default, with shape `a_sigmab`, default 10, and rate `b_sigmab`, default 1, now user-set) or `"half_cauchy"` on `sigma_b` itself with `sigma_b_scale` (default 1). The half-Cauchy draw is the auxiliary-variable Gibbs step already used for the spatial amplitude option. Default fits are bit-identical to the previous code, which `tests/testthat/test-sigma-b-prior.R` checks, along with the option reaching the one place `sigma_b` is drawn and leaving the spatial SD prior alone. The applied prior is saved in `infos$sigma_b_prior`.

The three priors compared here, on `sigma_b` with the truth at 0.35:

-   **InvGamma(10, 1)** on `sigma_b^2`, the default: prior mean of `sigma_b^2` 1/9, so `sigma_b` about 0.33 (median 0.32), 95% prior interval about 0.24 to 0.46. Shape 10 counts like 20 pseudo-observations, so with 20 B elements the data and the prior have equal weight on `sigma_b`.
-   **Half-Cauchy(1)** on `sigma_b`: median 1, no upper bound, 50% of the prior between 0.41 and 2.4. Weakly informative in the usual sense; the chains start `sigma_b` at 1 as before.
-   **InvGamma(2, 1/9)** on `sigma_b^2`: the same prior mean of `sigma_b^2` as the default, 1/9, with shape 2, which counts like 4 pseudo-observations and leaves the prior variance infinite. The prior median of `sigma_b` is 0.26 and its 95% prior interval about 0.14 to 0.68. It holds the location and changes only the weight, so it separates the prior's weight from its form.

## Design

Five cells at 100 replicates each, all at seed label `base` and all fitting `no_spatial`'s datasets with no spatial field simulated or fitted:

-   `no_spatial_sbhc`: half-Cauchy(1) on `sigma_b`, 1,000 + 1,000 iterations per chain. Pairs with `no_spatial`.
-   `no_spatial_sbig2`: InvGamma(2, 1/9) on `sigma_b^2`, 1,000 + 1,000. Pairs with `no_spatial` and with `no_spatial_sbhc`.
-   `no_spatial_sbhc_long`: half-Cauchy(1), 5,000 + 5,000, the package's default chain length. Pairs with `no_spatial_long`. This is the primary arm.
-   `no_spatial_sbig2_long`: InvGamma(2, 1/9), 5,000 + 5,000. Pairs with `no_spatial_long` and with `no_spatial_sbhc_long`.
-   `no_spatial_sbhc_gt4_long`: half-Cauchy(1) with 4 latent traits fitted, 5,000 + 5,000. Pairs with `no_spatial_gt4_long` and with `no_spatial_sbhc_long`, to see whether a free `sigma_b` makes the latent-trait probe redundant or adds to it.

The comparison arms are read from the earlier results files, not rerun: `no_spatial` and `no_spatial_long` from `../b-nonspatial-candidates-20261015/results/simstudy-20261016-004116.rds` (revision `c40f385`) and `no_spatial_gt4_long` from `../b-gt4-long-20261016/results/simstudy-20261016-074749.rds` (revision `815ba6f`, package code identical to `c40f385`). This branch's package code is `c40f385` plus the option, whose default path is bit-identical, so a rerun of `no_spatial` would reproduce the earlier one and is not spent.

All arms run through the study runner with 5 replicate processes and one thread per fit, 2 chains. The long cells carry their own MCMC settings, which replace the runner's.

## Pairing

Every arm draws `base`'s truths and fits exactly `no_spatial`'s datasets, because none of them changes anything the simulator sees. B, B0, G and `sigma_b` truths are therefore identical element for element across all eight arms compared. `tests/testthat/test-simstudy-sigma-b-prior.R` checks at replicate 1 that each cell hands the fitter `no_spatial`'s data and the priors and settings its name says; `mixing.R` re-simulates all 100 replicates, and `analyse.R` checks the truths in the results files before anything else.

## What is recorded

-   B: coverage overall and by absolute true effect (under 1, 1 to 2, 2 to 3, over 3), mean interval width, RMSE, slope of the posterior means on the truth, error-to-width ratio.
-   G and B0: coverage and width, as the cost side.
-   `sigma_b`: where its posterior settles (mean and interval across replicates), the coverage of the truth 0.35 by its interval, and its spread across replicates.
-   Paired contrasts, replicate by replicate, with standard errors, between each arm and its comparison arm, plus the interaction of the prior with chain length, (sbhc_long minus sbhc) minus (long minus no_spatial).
-   Mixing of B, B0 and `sigma_b` (bulk and tail effective sample size, R-hat) in 5 refitted replicates per new arm, as in the two B reports.

## Success criteria, written before the run

Standard errors are across replicates. The primary arm is `no_spatial_sbhc_long`, the half-Cauchy at the package's default chain length, against `no_spatial_long`.

-   **Restores:** the arm's 95% interval for B coverage (mean plus or minus 1.96 SE) contains 0.95, the paired gain over its comparison arm exceeds 2 SE, G coverage does not exceed 0.97, and B0 coverage is not lower than the comparison arm's by more than 2 SE.
-   **Partly restores:** the paired gain exceeds 2 SE but the interval excludes 0.95, or the gain is there and G over-covers past 0.97 or B0 loses more than 2 SE.
-   **No effect:** the paired gain is within 2 SE of zero.
-   **Large effects:** true effects over 3 are "mostly closed" at 0.90 or above (they cover at 0.685 in `no_spatial`, 0.762 in `no_spatial_long`, 0.895 in `no_spatial_gt4_long`).
-   **Form against weight:** if `no_spatial_sbig2` is within 2 SE of `no_spatial_sbhc` (and the same at the long length), the prior's weight is what matters; if the half-Cauchy exceeds the weak inverse-gamma by more than 2 SE, its form matters too.
-   **Chain length still needed:** `no_spatial_sbhc_long` exceeds `no_spatial_sbhc` by more than 2 SE.
-   **The probe is replaced:** `no_spatial_sbhc_gt4_long` is within 2 SE of `no_spatial_sbhc_long`. If instead the two add, the latent-trait count acts on something the `sigma_b` prior does not reach.
-   **Where `sigma_b` goes:** it "moves" if the mean of its posterior means across replicates differs from the default arm's 0.33 by more than 2 SE of the paired difference. A posterior that settles near the truth of 0.35 with intervals covering it at about 0.95 says the data identify `sigma_b` once the prior lets them; one that drifts well above says B's spread across species is wider than the trait regression and the default prior allow.

## Run

Launched detached from this directory, on mains power with Low Power Mode off and no other fits running, with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`. The log, provenance and results are ignored while the run is in flight so the runner records a clean tree. Expected cost from the earlier runs: about 30 seconds per short fit and about 90 seconds per long fit at 5 workers, so roughly two and a half hours for the 500 fits.

## Results

The run finished on 17 October 2026 at 18:11 with all 500 fits and no failures. A half-Cauchy(1) on `sigma_b` partly restores B: +0.015 (SE 0.006) at 1,000 + 1,000 and +0.016 (SE 0.004) at 5,000 + 5,000 on identical datasets, to 0.938 (0.927 to 0.949) at the package default chain length, with no cost to B0 or G, but short of 0.95 and short of the latent-trait probe, to which it adds nothing. The weak InvGamma(2, 1/9) changes nothing: within the inverse-gamma family `sigma_b`'s posterior is its prior. Under the half-Cauchy `sigma_b` is unidentified (intervals from 0.03 to 1.6). Chain length is still needed and adds. The findings are in [REPORT.md](REPORT.md).
