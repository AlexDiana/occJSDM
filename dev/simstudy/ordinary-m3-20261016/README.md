# B0's community mean in ordinary communities at 3 field samples per site, 16 October 2026

At 3 field samples per site, the rare-M3 run ([its report](../rare-m3-20261016/REPORT.md)) found that `b0_prior = "community"` alone removes the rare-dominated community's occupancy shift by every criterion set for it. Whether the option could become the default also depends on what it costs communities that are not rare-dominated. At 2 field samples per site the community-mean run ([its report](../b0-community-mean-20261016/REPORT.md)) found small losses in the parents: in `base_cm` the community occupancy error went from -0.023 to -0.047 (C9 failed) and the collection-intercept bias from -0.002 to +0.052 (C11 failed narrowly), and in `occupancy_cm` B0 coverage went from 0.935 to 0.909 (C8 failed). This run refits the two eDNA parents at 3 field samples per site under the default prior and under B0's community mean.

## Design

Four arms at 100 replicates each:

-   `base_m3` and `occupancy_m3`: the two-stage and one-stage parents (true B0 drawn from the simulator's own Normal(0, 1)) at 3 field samples per site, default prior.
-   `base_cm_m3` and `occupancy_cm_m3`: the same datasets with `b0_prior = "community"`.

Within each model the two priors fit identical datasets. Each arm takes its parent's seed label (`base`, `occupancy`), so, as `tests/testthat/test-simstudy-ordinary-m3.R` checks at replicate 1, the M = 3 datasets share everything drawn before the field samples with the M = 2 parents' datasets: B0, the occupancy predictor, the true occupancy states, `p`, `q`, `theta0` and the collection intercepts. Only the field samples, the collection covariate and the collection slopes differ. The M = 2 counterparts are replicates 1 to 100 of `base` and `occupancy` from the R = 300 paper study and `base_cm` and `occupancy_cm` from the community-mean run. The binary parent is not refitted: it has no field samples to vary.

The fits use the study runner's settings: 2 chains of 1,000 draws after 1,000 burn-in, 5 replicate processes, one thread per fit. To check mixing, the both-means run's `refit.R` refits replicates 1 to 10 of the four arms at the same settings, keeping the draws of `mu_b0`, B0, the collection intercepts, `theta0` and `p`; matching the stored rows exactly is the check that the refits are the study's fits. If C12 fails, replicates 1 to 5 of all four arms are refitted at 5,000 + 5,000. The default arms are included so that the community-minus-default comparison can be repeated at the longer chain length.

## Success criteria

Written on 16 October 2026 before any fit ran. "Default" is the default prior's result on the same M = 3 datasets. The M = 2 values in brackets are the community-mean run's, default first. Errors are estimate minus truth; community-level quantities are averaged over a replicate's 10 species first; the collection-intercept error is the posterior mean minus the truth on the logit scale.

No harm to ordinary communities (`base_cm_m3` against `base_m3`, `occupancy_cm_m3` against `occupancy_m3`):

-   C8. B0 coverage no more than 0.02 below the default's (M = 2: 0.939 -\> 0.927 and 0.935 -\> 0.909).
-   C9. Community mean occupancy error no more than 0.02 larger in absolute value (M = 2: -0.023 -\> -0.047 and -0.011 -\> -0.015).
-   C10. No true-occupancy band's absolute species occupancy bias more than 0.03 larger, over the bands 0 to 0.05, 0.05 to 0.1, 0.1 to 0.2, 0.2 to 0.3, 0.3 to 0.5, 0.5 to 0.7, 0.7 to 0.9 and 0.9 to 1 that hold at least 10 species-fits (M = 2, the band that worsened most: 0.7 to 0.9 in `base`, -0.169 -\> -0.191; 0.5 to 0.7 in `occupancy`, -0.050 -\> -0.050).
-   C11. Collection intercept, `p` and `theta0`: bias changed by at most 0.05, and coverage no more than 0.02 lower; B coverage no more than 0.02 lower (M = 2: collection-intercept bias -0.002 -\> +0.052 and +0.042 -\> +0.056; B coverage 0.936 -\> 0.930 and 0.925 -\> 0.920). The one-stage model has no `p`.

Mixing, from the refits of replicates 1 to 10:

-   C12. In each community-mean arm, `mu_b0` has rank-normalised R-hat under 1.05 in at least 9 of 10 replicates and median bulk ESS at least 100 of the 2,000 draws (M = 2: 5 and 4 of 10 replicates, median ESS 27 and 29).
-   C13. B0 and the collection intercepts mix no worse than under the default prior on the same M = 3 datasets: median bulk ESS at least 75% of the default arm's, and the share of elements with R-hat over 1.05 no more than 0.05 higher (M = 2, B0: median ESS 61 -\> 40 and 62 -\> 52, share 0.30 -\> 0.50 and 0.31 -\> 0.45).

Descriptive, with no pass mark:

-   D1. The third field sample: under each prior, the change from M = 2 to M = 3 in community occupancy error, B0 coverage and collection-intercept bias, paired by replicate, and whether the option's cost (community minus default) is smaller at M = 3 than at M = 2.
-   D2. `mu_b0`'s coverage of the generating mean, 0, and its estimate against the realised mean of the true B0.
-   D3. If the long chains run: how far each replicate's errors move between 1,000 + 1,000 and 5,000 + 5,000, and the community-minus-default differences of C8 and C9 at 5,000 + 5,000 on replicates 1 to 5.

The criteria are judged on the point estimates, as in the M = 2 run; paired standard errors are reported alongside, and a failure within two of them of its threshold is reported as narrow but still counts as a failure. `b0_prior = "community"` counts as harmless to ordinary communities at M = 3 if C8 to C11 hold in both arms. Together with the rare-M3 result, that would make it a candidate for the default. If C12 or C13 fail at 1,000 + 1,000 but the long chains leave the averaged results in place (D3), it would be a candidate that needs longer chains than the default prior to mix.

## Run

Launched detached from this directory, the draws refit and then the main run, so the laptop never runs more than 5 fitting workers:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'sh launch.sh draws; sh launch.sh main' > launch.log 2>&1 < /dev/null &
```

and, if C12 fails, after the main run:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh long > launch-long.log 2>&1 < /dev/null &
```

Logs, provenance, results and draws are ignored while the runs are in flight, so each records a clean tree.
