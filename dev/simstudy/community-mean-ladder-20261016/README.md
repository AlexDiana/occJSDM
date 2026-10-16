# The field-sample ladder with B0's community mean, 16 October 2026

Under the default prior, field samples per site undo the occupancy squeeze in an ordinary community: the slope of each species' estimated occupancy on its true occupancy is 0.44, 0.66, 0.86 and 0.95 at 2, 5, 10 and 20 field samples per site, while extra PCR replicates do not help ([the field-sample ladder](../field-sample-ladder/REPORT.md)). B0's community mean (`b0_prior = "community"`) removes the community-wide shift in rare-dominated communities and at 3 field samples per site costs ordinary communities nothing measurable ([rare-m3](../rare-m3-20261016/REPORT.md), [ordinary-m3](../ordinary-m3-20261016/REPORT.md)), but a squeeze remains, with a slope of 0.60 in the ordinary and 0.62 in the rare-dominated community at M = 3. Under the default prior the rare-dominated community's squeeze barely responds to field samples (slope 0.27 at M = 2, 0.26 at M = 3), and the ladder has never been run on such a community.

This run asks: once the community mean is on, how many field samples per site remove the remaining squeeze, in an ordinary and in a rare-dominated community, and how does that compare with the default prior?

## Design

Nine new arms, all two-stage (the `base` cell family: 100 sites, 10 species, 3 PCR replicates per field sample, a spatial field, traits), at 100 replicates with 2 chains of 1,000 draws after 1,000 burn-in, to match the existing arms:

-   Ordinary community (true B0 drawn from Normal(0, 1), as in the paper study), B0's community mean: `M5_cm`, `M10_cm` and `M20_cm`, at 5, 10 and 20 field samples per site.
-   Rare-dominated community (harness setting `b0_mean = -2.75`), default prior: `base_rare_m5`, `base_rare_m10` and `base_rare_m20`.
-   Rare-dominated community, B0's community mean: `base_rare_cm_m5`, `base_rare_cm_m10` and `base_rare_cm_m20`.

Existing arms, reused rather than refitted, replicates 1 to 100 of each:

-   Ordinary, default prior, M = 2, 5, 10 and 20: `M2`, `M5`, `M10` and `M20` from the field-sample ladder (`../field-sample-ladder`, merged into this branch from `claude/field-sample-ladder`). Its `K30` arm is not used.
-   Ordinary, M = 3, both priors: `base_m3` and `base_cm_m3` from `../ordinary-m3-20261016`.
-   Ordinary, M = 2, community mean: `base_cm` from `../b0-community-mean-20261016`, with its default-prior counterpart `base` from the R = 300 paper study (`../paper-study-20261015`).
-   Rare-dominated, M = 2: `base_rare` (default prior) from `../rare-dominated-20261015` and `base_rare_cm` from `../b0-community-mean-20261016`.
-   Rare-dominated, M = 3, both priors: `base_rare_m3` and `base_rare_cm_m3` from `../rare-m3-20261016`.

### Pairing

Every new arm takes its parent's seed label. The ordinary community-mean arms take the ladder's, `mladder`, so `M5_cm` fits exactly `M5`'s datasets, and likewise at 10 and 20. The rare arms take `base_rare`'s, `base`, so at every M they share B0, eta and the true occupancy states with `base_rare` and `base_rare_m3`, and each rare community-mean arm fits exactly its default arm's datasets. So within each community and M the two priors are compared on identical datasets, and along the rare ladder (M = 2, 3, 5, 10, 20) and the ordinary ladder (M = 2, 5, 10, 20 under the default prior, 5, 10, 20 under the community mean) the true occupancies are the same at every rung, so changes with M are paired by species.

The ordinary community's M = 2 community-mean rung (`base_cm`) and both M = 3 rungs come from runs at seed label `base`, a different set of 100 truths from the ladder's. The default prior has both at M = 2 (`base` and `M2`), which shows how far the change of truths alone moves each quantity; the ordinary community-mean ladder is read with that in mind.

`tests/testthat/test-simstudy-cm-ladder.R` checks the pairing at replicate 1 without fitting. `analyse.R` checks it on every replicate from the stored truths, as the earlier runs did, and writes `pairing-check.csv`: within each M the two priors' truths must agree in every block (B0, B, G, the collection coefficients, theta0, `p`, `q`) and in each species' true occupancy; across M within a seed label, B0 and the true occupancy must agree.

## What is reported

For each community and prior, by M (2, 3, 5, 10, 20), over replicates 1 to 100. Occupancy is the share of a replicate's 100 sites where a species is present: the truth from the simulated states, the estimate from the posterior probability of presence.

-   The slope of each species' estimated occupancy on its true occupancy, a least-squares fit over the 1,000 species-fits.
-   Occupancy bias (estimate minus truth) by band of true occupancy: 0 to 0.05, 0.05 to 0.1, 0.1 to 0.2, 0.2 to 0.3, 0.3 to 0.5, 0.5 to 0.7, 0.7 to 0.9 and 0.9 to 1, with standard errors that treat each replicate as one cluster.
-   The mean absolute error of each species' occupancy.
-   The community average occupancy error: per replicate, the mean estimated minus the mean true occupancy over the 10 species, then averaged over replicates.
-   B0 coverage of the 95% intervals, with a standard error across replicates, and the mean interval width. Under the default prior, B0 coverage fell to 0.873 at M = 20 in the ordinary community.
-   The collection-intercept bias: posterior mean minus truth on the logit scale, averaged over species and replicates.
-   The `p` bias: posterior mean minus truth.
-   In the community-mean arms, the rank-normalised R-hat and bulk ESS of `mu_b0`, from refits of replicates 1 to 10 that keep the draws (see *Mixing*). The M = 2 and M = 3 values come from the earlier runs' mixing tables.

Alongside, at each M, the community-mean arm minus the default arm on identical datasets, paired by replicate, for the community occupancy error, the mean absolute error, B0 coverage and the collection-intercept bias.

## When is the squeeze removed?

Written on 16 October 2026 before any fit ran. The ladder is discrete, so the answer for each community and prior is the smallest rung (2, 3, 5, 10 or 20 field samples per site) that meets each rule:

-   **Mostly removed**: occupancy slope at least 0.80, and in every band holding at least 10 species-fits the mean occupancy bias is within 0.10 of zero.
-   **Removed**: occupancy slope at least 0.90, and every such band within 0.05 of zero.

On the ordinary community under the default prior, these rules restate the ladder report's reading: mostly removed at 10 field samples per site (slope 0.86), removed at 20 (0.95). They are judged on the point estimates; band standard errors are reported alongside, and a failure within two of them of its threshold is reported as narrow but still counts as a failure. If no rung meets a rule, the answer is "more than 20". The REPORT also gives the community occupancy error, the mean absolute error and B0 coverage at the rung each rule names, since a slope near 1 does not by itself mean the intervals are right.

## Mixing

`../theta-community-mean-20261016/refit.R` refits replicates 1 to 10 of the six new community-mean arms at the study's settings, keeping the draws of `mu_b0`, B0, the collection intercepts, theta0 and `p`. Every refit must reproduce its arm's stored rows exactly, which shows the refits are the study's fits. `mu_b0` counts as converged in an arm if its R-hat is under 1.05 in at least 9 of the 10 replicates and its median bulk ESS is at least 100 of the 2,000 draws (criterion C12 of the M = 3 runs, which failed there at 1,000 + 1,000).

`mu_b0_check.R` applies that rule after the main run. If any arm fails, replicates 1 to 5 of each failing arm are refitted at 5,000 + 5,000, and the REPORT says how far each replicate's community occupancy error, mean absolute error, occupancy slope and B0 coverage move between the two chain lengths. The default-prior arms are not refitted: they have no `mu_b0`.

## Cost

An estimate, not a measurement: M = 3 fits took about 40 s and the ladder's 500 fits 137 minutes on 5 workers, which suggests about 1, 1.7 and 3.2 minutes per fit at M = 5, 10 and 20. That puts the 900 main fits at about 6 hours, the 60 mixing refits at about 25 minutes and the long chains, if needed, at up to an hour, on mains power. On battery a closed lid pauses the run without changing it.

## Run

`launch.sh` runs one step: `draws` (the mixing refits), `main` (the nine arms through the study runner) or `long` (the long chains for the arms listed in `draws/long-arms.txt`). Each step refuses to start while another simulation study or refit is running, if its output exists, or if the package code, tests or harness have uncommitted changes, and records its provenance. Launched detached from this directory, one step at a time, so the laptop never runs more than 5 fitting workers:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'sh launch.sh draws && sh launch.sh main && { /usr/local/bin/Rscript mu_b0_check.R . || sh launch.sh long; }' > launch.log 2>&1 < /dev/null &
```

The draws step ran from 22:53 to 23:10 on 16 October and the main step started at 23:10, both on 5 workers. Because the run would go overnight, Doug asked at about 23:30 for 6 workers. The main step was stopped just after a checkpoint and resumed from it on 6 workers, with the long step to follow on 6:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'RESUME=1 CORES=6 sh launch.sh main && { /usr/local/bin/Rscript mu_b0_check.R . || CORES=6 sh launch.sh long; }' >> launch.log 2>&1 < /dev/null &
```

The runner seeds each fit from its scenario and replicate alone, so neither the resume nor the worker count changes any fit; the refits, made on 5 workers, must still reproduce the stored rows exactly, which checks this.

Logs, provenance, results and draws are ignored while the runs are in flight, so each step records a clean tree. Then, from the repository root:

```sh
Rscript dev/simstudy/community-mean-ladder-20261016/analyse.R
Rscript dev/simstudy/community-mean-ladder-20261016/mixing.R
```
