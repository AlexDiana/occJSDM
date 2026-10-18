# Which latent term lets the community's occupancy intercepts drift, 18 October 2026

The 100-species mixing check ([b0-community-s100-20261018](../b0-community-s100-20261018/REPORT.md), branch claude/b0-community-s100) found that with `b0_prior = "community"` the community mean `mu_b0` failed to converge in 3 of 10 two-stage fits at the package default of 5,000 + 5,000 iterations: ordinary replicate 2 (R-hat 1.060, bulk ESS 32) and rare replicates 3 and 4 (R-hat 1.128 and 1.189, bulk ESS 13 and 8). This is a diagnostic of those three fits, run before any sampler move is written, to say what a move should actually move. The sampler is not changed.

## What the kept draws already show

In rare replicates 3 and 4 the two chains sit apart for all 5,000 kept draws. The difference between the chains is a near-uniform shift of +0.15 to +0.18 in every species' B0 (sd across species about 0.12, uncorrelated with the true occupancy or the true B0), with no compensating shift in the collection intercepts (-0.01), logit theta0 (-0.03) or logit p (0.00). In ordinary replicate 2 the chains differ by only 0.03 on average but by -0.40 in the species with true occupancy 0.05 to 0.1, a slow drift rather than two modes. The centred `mu_b0` Gibbs step is not the bottleneck: each species' B0 is pinned by its own data with precision about 20 (100 sites of Polya-Gamma weights of about 0.2) against a prior precision of 1, so the sd of the `mu_b0` draw given the B0, about 0.1, is already close to the posterior sd of `mu_b0`, 0.15 to 0.19.

So the partner of the shift, if there is one, is a term of the occupancy linear predictor that the run did not keep. The predictor is eta = B0 + X B + Ks Bs + U L (`computePsiCoef()` in R/jsdmfun.R): the species' intercepts, the covariate term, the spatial term (the spatial basis Ks at the drawn range, times the species' spatial coefficients Bs) and the latent factor term (site scores U times loadings L). The other possibility is that nothing compensates: the chains then differ in the number of occupied sites, and the direction is identified only weakly, through the latent states z.

## Design

-   Fits: `base_rare_cm_s100` replicates 3 and 4 and `base_cm_s100` replicate 2, refitted with `refit.R`, a variant of `dev/simstudy/theta-community-mean-20261016/refit.R`. Same datasets and same chains as the main run: workers seeded with `clusterSetRNGStream(cl, 20260727L)`, the same seeds through `draw_truth()`, nburn = 5000, and niter = 500 with nthin = 10, so each chain runs the main run's 10,000 iterations and stores every 10th draw after burn-in. `summarisedLatentPresences = FALSE`, so z and psi are stored per draw. Kept: `mu_b0`, B0, B, Bs with the range index `idx_ls`, U, L, the three scale parameters, the collection intercepts, theta0, p, z and psi, plus what rebuilds the terms (X, the spatial basis at every range of the grid, the support-point index) and the truth.
-   Storage check, before the run (`check-storage.R`, step `check`): the same dataset fitted twice in one process, once as the main run fitted it (nthin = 1, summarisedLatentPresences = TRUE) and once as `refit.R` fits it (nthin = 10, a tenth of the kept draws, summarisedLatentPresences = FALSE), for the same number of iterations. The draws at every 10th iteration must be identical, and the stored psi must equal logistic(eta) rebuilt from the stored terms (`terms.R`). Two settings: the 10-species rare arm, replicate 1, at 50 + 100, and the 100-species rare arm, replicate 3, at 20 + 40.
-   Reproduction check, after the run (`analyse.R`): the stored `mu_b0`, B0, collection intercept, theta0 and p draws against the main run's at every 10th iteration (draws/main in the worktree .claude/worktrees/nostalgic-mendel-8b7449, regenerable with that run's launch.sh). The differences must be 0.
-   3 worker processes with `VECLIB_MAXIMUM_THREADS=1` and one RcppParallel thread each, on R 4.6.1 with Accelerate's BLAS, with no other study running; about 35 minutes.

## What is computed

For each fit, from the stored draws (`analyse.R`):

-   The chain-2 minus chain-1 difference in the mean over sites, species and draws of each additive term of eta (B0, X B, Ks Bs, U L), of eta itself, of psi and of the total number of occupied sites; also of `mu_b0`, the mean collection intercept, mean logit theta0, mean logit p and the column means of U.
-   The same difference species by species: the shift in each species' B0 against the shift in the site-mean of each term for that species, with the correlation and the regression slope across species. A term that compensates has a slope near -1 and leaves the species' eta shifts near 0.
-   Within each chain, the correlation over draws of the mean B0, and of `mu_b0`, with the mean of each term, with mean psi and with the number of occupied sites, averaged over the two chains; and the rank-normalised R-hat and bulk ESS of each mean term.
-   Batch means (25 stored draws, 250 iterations) of each mean term per chain, plotted.

## Reading the result

Written before the run.

-   If a term's chain difference is about minus the B0 difference, species by species (slope near -1, eta shifts near 0), that term is the partner, and a sampler move should move `mu_b0` and all the B0 jointly against that term's mean, or centre that term.
-   If no term compensates, the shift goes straight into eta and psi and the chains differ in the number of occupied sites: the direction is identified only through z, and a move should propose a shift of `mu_b0` and all the B0 with z integrated out, or the conclusion is that only longer chains or a tighter `mu_b0` hyperprior help.
-   The recommendation goes in REPORT.md. The move itself is not implemented here.

## Run

`launch.sh` runs one step, named by its first argument (`check` or `main`). Each step refuses to start if its output exists, writes a provenance file (revision, dirty-file count, R version, `extSoftVersion()[["BLAS"]]`, the BLAS thread setting R sees, power state, settings) and keeps the laptop awake with `caffeinate -i` while it runs. Logs, provenance and draws are ignored, so each step records a clean tree. The check runs in the foreground; the main step is launched detached:

```sh
sh launch.sh check
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh main > launch-main.log 2>&1 < /dev/null &
```

`analyse.R` then writes the tables and plots, and the findings go in REPORT.md.
