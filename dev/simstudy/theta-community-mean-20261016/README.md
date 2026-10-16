# Community means for both intercepts, 16 October 2026

The community-mean run ([its report](../b0-community-mean-20261016/REPORT.md)) estimated the centre of B0's prior, `mu_b0`, instead of fixing it at 0. In rare-dominated eDNA communities that removed most of the occupancy shift (community occupancy error +0.17 to -0.05 in the two-stage arm, -0.01 in the one-stage arm), but the collection intercepts' error reversed, from -0.47 to +0.28 and from -0.37 to +0.22. With B0's centre free, the occupancy-collection split at community scale was decided by the collection intercepts' own fixed Normal(0, 1) prior, which sits above 98 to 99% of their true values (mean -0.64). This run estimates the collection intercepts' centre as well, so that neither side has a fixed prior centre, and asks whether the data at 2 field samples per site then identify the split at all.

## The model change

The experimental option `listPriors = list(theta_intercept_prior = "community")`, added on this branch and off by default, gives each species' collection intercept, `beta_theta[1, s]`, a `Normal(mu_theta, 1)` prior in place of `Normal(0, 1)`, and estimates `mu_theta` under a `Normal(0, 2^2)` hyperprior (`mu_theta_sd`, default 2, the same as `mu_b0_sd`). Each iteration, after the collection coefficients are drawn, `mu_theta` is drawn from its conjugate normal conditional given the species' collection intercepts, and the next iteration's collection draw uses it as the intercepts' prior mean. The collection intercepts' prior enters the sampler in one place, the prior-mean vector `b_betatheta` passed to the collection sampler; no other step integrates the collection coefficients out, so there is no second site to keep in agreement. The name is `mu_theta`, not `mu_theta0`, because `theta0` is already the field-contamination rate.

One detail differs from the default sampler. `sample_betatheta_cpp_parallel()` leaves a species' collection coefficients unchanged in an iteration where the species is absent from every site, although with no collection data their conditional is their prior. Under the community option a stale intercept would feed `mu_theta`'s draw, so such species are drawn from the prior the sampler used. Default fits do not take this path.

The draws are returned in `results_output$mu_theta_output`, a niter by nchain array present only when the option is on; `returnConvergenceDiagnostics()` reports them, and `infos$collection_intercept_prior` records the prior. The option can be combined with `b0_prior = "community"`, which is how this run uses it. With the option off nothing changes: the sampler receives the same zero mean, no random number is drawn for `mu_theta`, and no output or `infos` element is added. [check_default_identity.R](check_default_identity.R) checks this on twelve seeded fits against the branch base, `cf9c133`, B0-only community fits included: all twelve whole fit objects are identical apart from `infos$execution`, which records timings and process ids, and a control fit that uses the option differs, as it must ([default-identity.csv](default-identity.csv)).

## Design

Six arms at 100 replicates each, all fitted with `b0_prior = "community"` and `theta_intercept_prior = "community"` through the harness's fit-only settings `fit_b0_prior` and `fit_theta_intercept_prior`:

-   `base_rare_cm_both` and `occupancy_rare_cm_both`: the rare-dominated eDNA arms (true B0 drawn around -2.75).
-   `base_cm_both` and `occupancy_cm_both`: their parents (true B0 drawn around 0).
-   `base_rare_cm_both_wide` and `occupancy_rare_cm_both_wide`: the rare arms again with both hyperpriors widened from sd 2 to sd 5 (`mu_b0_sd = 5`, `mu_theta_sd = 5`), to test whether the data or the hyperpriors decide the split.

Every arm takes its counterpart's seed label and B0 shift, so its datasets are bit-identical to those of two arms already run (`tests/testthat/test-simstudy-theta-community.R` checks this at replicate 1):

-   The default-prior arms: replicates 1 to 100 of `base` and `occupancy` from the R = 300 paper study at `de4389c`, and `base_rare` and `occupancy_rare` from the rare-dominated run at `b0d5811`.
-   The B0-only community arms `base_cm`, `occupancy_cm`, `base_rare_cm` and `occupancy_rare_cm` from the community-mean run at `62b3138`.

Binary cells are not run: they have no collection stage. The fits use the study runner's settings: 2 chains of 1,000 draws after 1,000 burn-in, 5 replicate processes, one thread per fit.

The runner stores interval rows, not draws. To check mixing, `refit.R` refits replicates 1 to 10 of the six arms at the same settings and keeps the draws of `mu_b0`, `mu_theta`, B0, the collection intercepts, `theta0` and `p`. Their default and B0-only counterparts' draws were kept by the community-mean run and summarised in its [mixing.csv](../b0-community-mean-20261016/mixing.csv), which this run compares against. The same step refits replicates 1 to 3 of the eight comparison arms on this branch; matching their stored rows exactly is the study-level check that the change leaves default and B0-only fits as they were. If `mu_b0` or `mu_theta` mixes poorly at 1,000 + 1,000 (C12 fails), replicates 1 to 5 of the four main arms are refitted at the package default of 5,000 + 5,000.

The truth for `mu_theta` is the realised mean of the species' true collection intercepts, on the fit's centred scale: the simulator draws each species' collection rate uniformly between 0.2 and 0.5, not intercepts around a mean.

## Success criteria

Written on 16 October 2026 before any fit ran. "Default" and "B0 only" are the results on the same datasets from the community-mean report. Errors are estimate minus truth; community-level quantities are averaged over a replicate's 10 species first. The collection-intercept error is the mean of the posterior means on the logit scale, as before.

Do both community means remove the block shift in the rare-dominated eDNA arms (`base_rare_cm_both`, `occupancy_rare_cm_both`)?

-   C1. Community mean occupancy error within ±0.05 (default +0.167 and +0.171; B0 only -0.050 and -0.013).
-   C2. Community mean B0 error within ±0.5 logit (default +1.89 and +1.87; B0 only -0.46 and -0.01).
-   C3. Species occupancy bias under +0.10 in each of the true-occupancy bands up to 0.05, 0.05 to 0.1 and 0.1 to 0.2 (default +0.30, +0.26 and +0.19 in `base_rare`; B0 only +0.03, +0.00 and -0.03 in `base_rare_cm`).
-   C4. B0 coverage at least 0.85 (default 0.422 and 0.393; B0 only 0.917 and 0.940).
-   C5. Collection-intercept bias within ±0.15 (default -0.47 and -0.37; B0 only +0.28 and +0.22).
-   C6. `mu_b0`'s 95% interval contains the generating mean, -2.75, in at least 85% of replicates (B0 only 0.92 and 0.97).
-   C14. `mu_theta`'s 95% interval contains the realised mean collection intercept in at least 85% of replicates, here and in the parents.

The shift counts as removed if C1, C2 and C5 hold in both arms. It counts as reduced relative to B0 only if C1 and C2 hold and the absolute collection-intercept error is at most half of B0 only's (0.14 and 0.11).

No harm in the parents (`base_cm_both` and `occupancy_cm_both` against `base` and `occupancy`):

-   C8. B0 coverage no more than 0.02 below the default's (0.939 and 0.935).
-   C9. Community mean occupancy error no more than 0.02 larger in absolute value (default -0.023 and -0.011).
-   C10. No true-occupancy band's absolute species occupancy bias more than 0.03 larger.
-   C11. Collection intercept, `p` and `theta0`: bias changed by at most 0.05, and coverage no more than 0.02 lower; B coverage no more than 0.02 lower.

Mixing, from the refits of replicates 1 to 10 of the four main arms:

-   C12. `mu_b0` and `mu_theta` each have rank-normalised R-hat under 1.05 in at least 9 of the 10 replicates of each arm, and median bulk ESS at least 100 of the 2,000 draws.
-   C13. B0 and the collection intercepts mix no worse than under the default prior on the same datasets: in each arm, median bulk ESS at least 75% of the default counterpart's, and the share of elements with R-hat over 1.05 no more than 0.05 higher.

Do the data at 2 field samples per site identify the split? If they did not, the community-wide position along the occupancy-collection trade-off would be decided by the two hyperpriors. Both are centred at 0, so widening both by the same factor would leave the posterior centred near the same place but spread it along the trade-off roughly in proportion to the hyperprior sd. If the data decide the split, widening changes little. From the `_wide` arms against `base_rare_cm_both` and `occupancy_rare_cm_both`, paired by replicate:

-   C15. Widening both hyperpriors from sd 2 to sd 5 changes the community mean occupancy error by less than 0.03 and the collection-intercept error by less than 0.15, on average over the 100 replicates, and widens `mu_theta`'s 95% interval by less than 25% on average (an unidentified split would widen it by up to 150%).
-   C16, descriptive. `mu_theta`'s posterior sd against its hyperprior sd of 2, against the floor of about 0.31 it would have if every species' intercept were known (10 species, each with prior sd 1), and the within-fit posterior correlation of `mu_b0` and `mu_theta`, from the refits.

The option pair is a success in this study if the shift is removed (C1, C2 and C5 in both rare eDNA arms) with no harm to the parents (C8 to C11) and acceptable mixing (C12 and C13). The data identify the split if C15 holds in both rare arms. C3, C4, C6 and C14 describe how far the improvement reaches into individual species and the two community means.

## Run

Each step is launched detached from this directory with `launch.sh`, one at a time so the laptop never runs more than 5 fitting workers:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'sh launch.sh draws; sh launch.sh main' > launch.log 2>&1 < /dev/null &
```

and, if C12 fails, after the main run:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh long > launch-long.log 2>&1 < /dev/null &
```

Monitor with `tail -f draws.log`, `tail -f main.log` and `tail -f long.log`. Logs, provenance, results and draws are ignored while the runs are in flight, so each records a clean tree.
