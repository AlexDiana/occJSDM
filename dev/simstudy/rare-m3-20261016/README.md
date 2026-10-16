# Rare-dominated communities at 3 field samples per site, 16 October 2026

The both-means run ([its report](../theta-community-mean-20261016/REPORT.md)) fitted rare-dominated eDNA communities with 2 field samples per site. With both community means estimated, the community's occupancy came out close to right, but the collection intercepts undershot (-0.52 two-stage, -0.24 one-stage). The error sat in the rarest species, whose own occupancy-collection split each species' few detections could not settle. Three field samples per site is a more realistic design, and a third sample is exactly the information that separates a species' occupancy from its collection rate. This run refits the two rare-dominated arms at M = 3 under all three priors: the default, B0's community mean alone, and both community means.

## Design

Six arms at 100 replicates each:

-   `base_rare_m3` and `occupancy_rare_m3`: the two-stage and one-stage rare-dominated arms (true B0 drawn around -2.75) at 3 field samples per site, default prior.
-   `base_rare_cm_m3` and `occupancy_rare_cm_m3`: the same datasets with `b0_prior = "community"`.
-   `base_rare_cm_both_m3` and `occupancy_rare_cm_both_m3`: the same datasets with `b0_prior` and `theta_intercept_prior` both `"community"`.

Within each model the three priors fit identical datasets. Each arm also takes the seed label of its M = 2 counterpart (`base_rare`, `occupancy_rare`, their `_cm` and `_cm_both` arms). Measured at replicates 1 to 3 before this run and tested in `tests/testthat/test-simstudy-rare-m3.R`, that makes the M = 3 datasets share everything drawn before the field samples with the M = 2 datasets: B0, the occupancy predictor, the true occupancy states, `p`, `q`, `theta0` and the collection intercepts. Only the field samples, the collection covariate and the collection slopes differ. So species occupancy and B0 can be compared across M replicate by replicate, while the collection intercepts' truth on the fit's centred scale differs slightly between M = 2 and M = 3. The M = 2 arms come from the rare-dominated run (default), the community-mean run (B0 only) and the both-means run.

The fits use the study runner's settings: 2 chains of 1,000 draws after 1,000 burn-in, 5 replicate processes, one thread per fit. To check mixing, the both-means run's `refit.R` refits replicates 1 to 10 of the six arms at the same settings, keeping draws. If C12 fails, replicates 1 to 5 of the four community-mean arms are refitted at 5,000 + 5,000.

## Success criteria

Written on 16 October 2026 before any fit ran. "M = 2" is the same prior's result on the paired M = 2 datasets. Errors are estimate minus truth; community-level quantities are averaged over a replicate's 10 species first; the collection-intercept error is the mean of the posterior means on the logit scale.

Do both community means remove the block shift at M = 3 (`base_rare_cm_both_m3`, `occupancy_rare_cm_both_m3`)?

-   C1. Community mean occupancy error within ±0.05 (M = 2: +0.033 and +0.040).
-   C2. Community mean B0 error within ±0.5 logit (M = 2: +0.40 and +0.47).
-   C3. Species occupancy bias under +0.10 in each of the true-occupancy bands up to 0.05, 0.05 to 0.1 and 0.1 to 0.2 (M = 2: +0.144, +0.106 and +0.053; +0.136, +0.104 and +0.057).
-   C4. B0 coverage at least 0.85 (M = 2: 0.936 and 0.919).
-   C5. Collection-intercept bias within ±0.15 (M = 2: -0.52 and -0.24).
-   C6. `mu_b0`'s 95% interval contains -2.75 in at least 85% of replicates (M = 2: 0.91 and 0.84).
-   C14. `mu_theta`'s 95% interval contains the realised mean collection intercept in at least 85% of replicates (M = 2: 0.89 and 0.95).

The shift counts as removed by both means if C1, C2 and C5 hold in both arms.

Does B0's mean alone suffice at M = 3 (`base_rare_cm_m3`, `occupancy_rare_cm_m3`)? The same C1, C2, C3, C4, C5 and C6 (M = 2: occupancy error -0.050 and -0.013, B0 error -0.46 and -0.01, collection-intercept bias +0.28 and +0.22). It counts as sufficient if C1, C2 and C5 hold in both arms.

What does the third field sample do (each M = 3 arm against its M = 2 counterpart, paired by replicate)?

-   C17. Under both means, the absolute collection-intercept error at M = 3 is at most half its M = 2 value (0.26 and 0.12).
-   C18, descriptive. Under each prior, the change in community occupancy error, in the low-band occupancy bias and in B0 coverage. With the default prior at M = 2 the community occupancy error was +0.167 and +0.171.

Mixing, from the refits of replicates 1 to 10:

-   C12. In each community-mean arm, `mu_b0` (and `mu_theta` in the both-means arms) has rank-normalised R-hat under 1.05 in at least 9 of 10 replicates and median bulk ESS at least 100 of the 2,000 draws.
-   C13. B0 and the collection intercepts mix no worse than under the default prior on the same M = 3 datasets: median bulk ESS at least 75% of the default arm's, and the share of elements with R-hat over 1.05 no more than 0.05 higher.

Both means are a success at M = 3 if C1, C2 and C5 hold in both arms and C12 and C13 hold. The parents are not rerun at M = 3, so this run cannot show that the community means leave common species unharmed.

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

## Results

The draws step finished on 16 October 2026 at 14:10 (60 fits), the main run at 15:46 with all 600 fits and no failures, and the long-chain step at 15:59 (20 fits). At 3 field samples per site, B0's community mean alone meets every criterion set for it in both rare arms: community occupancy error -0.022 and -0.007, B0 error -0.11 and +0.07, collection-intercept bias +0.14 and +0.13 against a limit of 0.15. Both community means still fail C5, the collection intercepts undershooting by 0.38 and 0.22, and the rarest species are again estimated too common. The default prior's bias barely moves with a third field sample. Mixing improves but still fails at the study's chain length. The findings are in [REPORT.md](REPORT.md).
