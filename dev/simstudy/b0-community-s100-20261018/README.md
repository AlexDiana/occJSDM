# B0's community mean at 100 species, 18 October 2026

Simulation studies from 15 to 17 October 2026 found that an estimated community mean for the species' occupancy intercepts, `b0_prior = "community"`, removes a community-wide occupancy bias in rare-dominated communities at no measurable cost to ordinary ones ([rare-M3 report](../rare-m3-20261016/REPORT.md), [ordinary-M3 report](../ordinary-m3-20261016/REPORT.md), [field-sample ladder](../community-mean-ladder-20261016/REPORT.md)). The code is PR #36 on AlexDiana/occJSDM, and its decision 2 is whether the option becomes the default. Every one of those studies simulated 10 species. This is a mixing check at a realistic 100 species before that decision, not a full simulation study.

At 10 species and 3 field samples per site, `mu_b0` mixed slowly at the study chain length of 1,000 + 1,000. In the two-stage arms its rank-normalised R-hat was under 1.05 in 5 of 10 replicates of the ordinary community and 6 of 10 of the rare-dominated one, with median bulk ESS 38 and 36 of 2,000 draws. At the package default of 5,000 + 5,000, on replicates 1 to 5, it converged by R-hat (at most 1.026 ordinary, 1.041 rare), with bulk ESS 115 to 366 (ordinary) and 71 to 199 (rare) of 10,000 draws (`long-chains.csv` in the two M = 3 run directories).

Why 100 species might differ: `mu_b0` is drawn by Gibbs given the species' B0, from a normal whose sd is close to `sigma_b0 / sqrt(S)` (`sigma_b0` is 1), so with 100 species each step moves `mu_b0` about a third as far as with 10. If its posterior narrowed by the same factor, mixing would be unchanged: in a plain normal hierarchy the autocorrelation of this centred sampler depends on how firmly each species' data fix its own intercept, not on the number of species. But here each species' B0 trades off against its collection intercept (for a species found at few sites, high occupancy with a low collection rate and low occupancy with a high one fit the data about equally well), so the direction in which `mu_b0`, all the B0 and all the collection intercepts shift together may not narrow as fast as the Gibbs step does. If it does not, `mu_b0` mixes more slowly at 100 species than at 10.

## Design

-   Arms: `base_cm_s100`, an ordinary community (true B0 drawn from Normal(0, 1)), and `base_rare_cm_s100`, a rare-dominated one (Normal(-2.75, 1)). They are defined in `tests/testthat/helper-simstudy.R` as `base_cm_m3` and `base_rare_cm_m3` with 100 species in place of 10, and checked by `tests/testthat/test-simstudy-s100.R`. Two-stage model, 100 sites, 3 field samples per site, 2 primers with 3 PCR replicates each, 2 latent factors, 2 latent traits, a spatial field, 2 occupancy covariates and 1 collection covariate. Fitted with `b0_prior = "community"`, so `mu_b0` has a `Normal(0, 2^2)` hyperprior, and every other prior at its default.
-   Pairing: both arms take seed label `s100`, so at each replicate the rare arm's true B0 is the ordinary arm's minus 2.75, species by species, and everything else drawn before the occupancy states is shared. At replicate 1 the mean true occupancy is 0.48 in the ordinary community and 0.17 in the rare one.
-   Chains: replicates 1 to 5 of each arm, 2 chains of 5,000 burn-in and 5,000 kept draws (the package default, no thinning), so 10,000 draws per fit.
-   Fitting: `dev/simstudy/theta-community-mean-20261016/refit.R`, which fits as the study runner does (workers seeded with `clusterSetRNGStream(cl, 20260727L)`, one RcppParallel thread, `devtools::load_all()`) and keeps the draws of `mu_b0`, B0, the collection intercepts, `theta0` and `p`. 5 worker processes with `VECLIB_MAXIMUM_THREADS=1`, on R 4.6.1 with Accelerate's BLAS, and no other study running.
-   The 10-species comparison: the long-chain steps of the two M = 3 runs kept their draws only in worktrees that no longer exist, and their committed summaries lack the correlations and R-hat shares reported below. Step `s10` refits replicates 1 to 5 of `base_cm_m3` and `base_rare_cm_m3` at 5,000 + 5,000 with the same script. Since those runs the package has changed only in documentation, so the refits should reproduce those runs' `long-chains.csv` values for these two arms exactly. That is checked before they are used.

## Timing

Step `timing` fits replicates 1 to 5 of both 100-species arms at 200 + 200 on 5 workers before the main run. A fit at 5,000 + 5,000 runs 25 times as many iterations per chain, so its projected time is 25 times the short fit's, and the 10 fits take two rounds on 5 workers. If the projected wall time is over about 3 hours, the main run uses 3 replicates per community instead of 5, and this README records the timing and the reason before it starts.

## Rule

Written on 18 October 2026 before any fit ran.

`mu_b0` mixes well enough at 100 species if, in every replicate of both arms at 5,000 + 5,000, its rank-normalised R-hat (`posterior::rhat`) is under 1.05 and its bulk ESS (`posterior::ess_bulk`) is at least 100 of the 10,000 draws.

The ESS bar is applied to every replicate, the stricter reading. Under it the 10-species rare-dominated arm at 5,000 + 5,000 would itself have failed (bulk ESS 78 and 71 in replicates 2 and 3). So the report also gives the result with the ESS bar applied to the median replicate, and applies both readings to the `s10` refits, so that a failure at 100 species can be told apart from one already present at 10.

## Reported, with no pass mark

For each 100-species arm, with the `s10` refits alongside:

-   `mu_b0`: R-hat, bulk ESS and tail ESS per replicate, and trace plots of both chains.
-   B0 and the collection intercepts: per replicate, the median R-hat, the share of species with R-hat over 1.05, and the median and minimum bulk ESS.
-   The slow direction: within each fit, the correlation of `mu_b0` with the species' mean B0 and of the mean B0 with the mean collection intercept, computed within each chain and averaged over the two, and the bulk ESS and R-hat of those two means.
-   `mu_b0`'s posterior mean against the realised mean of the true B0, and its posterior sd, expected to be about a third of the 10-species value (1 / sqrt(10) is 0.32).
-   The community occupancy error, as a sanity check: the mean estimated minus the mean true occupancy over sites and species.
-   Elapsed time per fit.

## Run

`launch.sh` runs one step, named by its first argument (`timing`, `s10` or `main`, where `main` takes the number of replicates per arm as its second argument). Each step refuses to start if its output folder exists, writes a provenance file (revision, dirty-file count, R version, `extSoftVersion()[["BLAS"]]`, the BLAS thread setting R sees, power state, settings) and keeps the laptop awake with `caffeinate -i` while it runs. Logs, provenance and draws are ignored while the runs are in flight, so each step records a clean tree. Launched detached from this directory:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh timing > launch-timing.log 2>&1 < /dev/null &
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'sh launch.sh s10; sh launch.sh main 5' > launch-main.log 2>&1 < /dev/null &
```

`analyse.R` then writes the summaries and trace plots, and the findings go in REPORT.md.
