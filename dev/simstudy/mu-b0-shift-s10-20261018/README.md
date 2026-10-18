# The collapsed shift step at 10 species: design, 18 October 2026

The collapsed shift of `mu_b0` and all the occupancy intercepts ([mu-b0-shift-20261018](../mu-b0-shift-20261018/REPORT.md), commit `a43682f`) was tested only at 100 species. Before it goes into PR #36, this run checks it in the paper's own design, 10 species, where `mu_b0` was already close to the convergence bar without it. It asks two things: does the step leave the posterior where it was, and does it leave `mu_b0`'s mixing at least as good.

## Design

-   Arms: `base_cm_m3` (ordinary community, true B0 drawn from Normal(0, 1)) and `base_rare_cm_m3` (rare-dominated, Normal(-2.75, 1)), both two-stage, 100 sites, 10 species, 3 field samples per site, `b0_prior = "community"`, with the shift step at its defaults.
-   Replicates 1 to 5 of each, 2 chains of 5,000 burn-in and 5,000 kept draws, draws kept, refitted with `../mu-b0-shift-20261018/refit-s100.R` unchanged.
-   The comparison is the step-off refits of the same replicates at the same settings, the `s10` step of the 100-species check (`../b0-community-s100-20261018`, run 18 October 2026 11:08 to 11:15 at `fd2689f`). Both use the same datasets and chain seeds. With the step on, the chains draw extra random numbers, so they are different chains of the same posterior.
-   One RcppParallel thread and one BLAS thread per fit (`VECLIB_MAXIMUM_THREADS=1`), 5 workers, R 4.6.1 with Accelerate BLAS, launched detached by `launch.sh` from a clean tree.

## Rules, set before any fit

-   **R1, mixing by the 100-species check's rule.** `mu_b0` has R-hat under 1.05 in every fit, and bulk ESS of at least 100 in 10,000 draws in every fit (strict reading) or in the arm's median (median reading). Reported for both arms and both readings, with and without the step.
-   **R2, mixing not worse.** In each arm, the median bulk ESS of `mu_b0` with the step is at least that without it.
-   **R3, posterior unchanged.** In every fit, the posterior means of `mu_b0` and of the species' mean B0 with and without the step differ by less than 3 combined Monte Carlo standard errors, `sqrt(mcse_on^2 + mcse_off^2)`.
-   **Descriptive, not a rule.** The community occupancy error and B0 coverage with and without the step (50 species-fits per arm are too few for a rule on coverage), and the step's acceptance rates and adapted proposal sd.

The step passes this check if R2 and R3 hold in both arms. R1 is reported for comparison with the 100-species check; the step-off fits already fail its strict reading in the rare arm.

## Files

-   `launch.sh` runs the step; `analyse.R` writes `fits.csv`, `before-after.csv` and `rule.csv` from `draws/s10` and the step-off draws.
-   Logs, provenance and draws are listed in `.gitignore`, so the launch records a clean tree.
