# A community mean for the occupancy intercepts, 16 October 2026

The rare-dominated run ([its report](../rare-dominated-20261015/REPORT.md)) found that in eDNA fits of a community whose species are mostly rare, every rare species is estimated at about a third of sites whatever its true occupancy. The true mean occupancy intercept of -2.72 was estimated at -0.82, and the community's mean occupancy at 0.37 against a true 0.20. The cause is B0's fixed Normal(0, 1) prior, centred on 50% occupancy, which decides the weakly identified split between how often a species is present and how often a sample catches it. That report asked whether B0 should instead shrink towards an estimated community mean, as in Hmsc and spOccupancy. This run tests that change.

## The model change

The experimental option `listPriors = list(b0_prior = "community")`, added on this branch and off by default, gives each species' B0 a `Normal(mu_b0, sigma_b0^2)` prior, with `sigma_b0` still 1 by default, and estimates `mu_b0` under a `Normal(0, 2^2)` hyperprior (`mu_b0_sd`, default 2). Each iteration `mu_b0` is drawn from its conjugate normal conditional given the species' B0 draws, after the coefficient block. The prior mean reaches all five places that encode B0's prior: both C++ coefficient samplers, the R reference sampler, and both spatial range log-weight functions, which integrate B0 out before the range is chosen. With the option off every sampler receives a zero mean, which is what it used before. The collection intercepts keep their fixed Normal(0, 1) prior, and `sigma_b0` is not estimated.

## Design

Six arms at 100 replicates each, all fitted with `b0_prior = "community"` through the harness's fit-only setting `fit_b0_prior`, each paired with a default-prior arm that has already run on identical datasets:

-   `base_rare_cm`, `occupancy_rare_cm` and `binary_rare_cm`: the rare-dominated arms (true B0 drawn around -2.75), against `base_rare`, `occupancy_rare` and `binary_rare` from the rare-dominated run at `b0d5811`.
-   `base_cm`, `occupancy_cm` and `binary_cm`: their parents (true B0 drawn around 0), against replicates 1 to 100 of `base`, `occupancy` and `binary` from the R = 300 paper study at `de4389c`.

Both comparison runs used package code identical to this branch's base, `34eac47`. Every arm takes its counterpart's seed label and B0 shift, so its datasets are bit-identical to that counterpart's (`tests/testthat/test-simstudy-b0-community.R` checks this at replicate 1). The fits use the study runner's settings: 2 chains of 1,000 draws after 1,000 burn-in, 5 replicate processes, one thread per fit.

The runner stores interval rows, not draws. To check mixing, `refit_draws.R` refits replicates 1 to 10 of all twelve arms, the six community arms and their six default counterparts, at the same settings and keeps the draws of `mu_b0`, B0, the collection intercepts, `theta0` and `p`. The default refits also show whether this branch with the option off reproduces the earlier runs' rows exactly, which is the study-level check that the option leaves default fits unchanged.

## Success criteria

Written on 16 October 2026 before any fit ran. "Now" is the default prior's result on the same datasets, from the rare-dominated report. Errors are estimate minus truth; community-level quantities are averaged over a replicate's 10 species first.

Does the community mean remove the block shift in the rare-dominated eDNA arms (`base_rare_cm` against `base_rare`, `occupancy_rare_cm` against `occupancy_rare`)?

-   C1. Community mean occupancy error within ±0.05 (now +0.167 and +0.171).
-   C2. Community mean B0 error within ±0.5 logit (now +1.89 and +1.87).
-   C3. Species occupancy bias under +0.10 in each of the true-occupancy bands up to 0.05, 0.05 to 0.1 and 0.1 to 0.2 (now +0.30, +0.26 and +0.19 in `base_rare`).
-   C4. B0 coverage at least 0.85 (now 0.422 and 0.393).
-   C5. Collection-intercept bias within ±0.15 (now -0.47 and -0.37), so the fit no longer lowers the collection rate to keep B0 near 0.
-   C6. `mu_b0`'s 95% interval contains the generating mean, -2.75, in at least 85% of replicates.

The shift counts as removed if C1, C2 and C5 hold in both arms, and as reduced if each of those errors is at least halved.

In the binary rare arm (`binary_rare_cm` against `binary_rare`):

-   C7. B0 coverage at least 0.85 (now 0.449) and community mean B0 error within ±0.3 (now +0.78).

No harm in the parents (`base_cm`, `occupancy_cm` and `binary_cm` against `base`, `occupancy` and `binary`):

-   C8. B0 coverage no more than 0.02 below the default's (now 0.939, 0.935 and 0.860).
-   C9. Community mean occupancy error no more than 0.02 larger in absolute value (now -0.023 and -0.011).
-   C10. No true-occupancy band's absolute species occupancy bias more than 0.03 larger.
-   C11. Collection intercept, `p` and `theta0`: bias changed by at most 0.05, and coverage no more than 0.02 lower; B coverage no more than 0.02 lower.

Mixing, from the refits of replicates 1 to 10:

-   C12. `mu_b0`'s rank-normalised R-hat under 1.05 in at least 9 of the 10 replicates of each community arm, and its median bulk ESS at least 100 of the 2,000 draws.
-   C13. B0 mixes no worse than under the default prior on the same datasets: in each community arm, B0's median bulk ESS at least 75% of its default counterpart's, and the share of B0 elements with R-hat over 1.05 no more than 0.05 higher.

The option is a success in this study if the shift is removed (C1, C2 and C5 in both rare eDNA arms) with no harm to the parents (C8 to C11) and acceptable mixing (C12 and C13). C3, C4, C6 and C7 describe how far the improvement reaches into individual species.

## Run

Launched detached from this directory, the main run and the draws refit side by side, with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch-draws.sh > launch-draws.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log` and `tail -f draws.log`. Logs, provenance, results and draws are ignored while the runs are in flight so each records a clean tree.

A third run was added once the first refits showed `mu_b0` mixing slowly: `refit_long.R` refits replicates 1 to 5 of the four eDNA community arms and the two default rare arms at 5,000 + 5,000. It was launched from a copy outside the tree, so that the main run would still record a clean tree, with its outputs under the ignored `draws/long/`; `long-provenance.txt` records the copy's checksum, which matches the committed `refit_long.R`:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch-long.sh <path to refit_long.R> > launch-long.log 2>&1 < /dev/null &
```

Three runs of 5 workers each were too many for the laptop, so the long run was paused from 10:07 until the main run finished at 11:06, leaving 5 workers in total. Keep to 5 workers in total, not per run.

## Results

The main run finished on 16 October 2026 at 11:06 with all 600 fits and no failures; the refits finished at 10:05 (120 fits) and 11:13 (30 fits). The community mean removes most of the occupancy and B0 shift in the rare-dominated eDNA arms: community occupancy error +0.17 to -0.05 (two-stage) and -0.01 (one-stage), B0 coverage about 0.4 to about 0.93. It does not meet the success definition, because the collection intercepts' error reverses from -0.47 to +0.28 instead of shrinking, the collection intercepts' own fixed prior now deciding the split. The parents lose a little, and `mu_b0` mixes slowly in the eDNA models at the study's chain length, though longer chains confirm the averaged results. The findings are in [REPORT.md](REPORT.md).
