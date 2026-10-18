# The collapsed shift step at 10 species: report, 18 October 2026

**Verdict: the step passes.** In the paper's design, 10 species at 3 field samples per site, the collapsed shift of `mu_b0` and all the occupancy intercepts leaves the posterior where it was and makes `mu_b0` mix better in both communities. Both rules set in [README.md](README.md) hold: R2 (mixing not worse) and R3 (posterior unchanged). R1, the 100-species check's rule, now passes in both communities on both readings, where the rare-dominated community failed its strict reading without the step.

## What ran

Replicates 1 to 5 of `base_cm_m3` and `base_rare_cm_m3` (two-stage, 100 sites, 10 species, 3 field samples per site, `b0_prior = "community"`) at 2 chains of 5,000 + 5,000 with the step at its defaults, refitted with `../mu-b0-shift-20261018/refit-s100.R` unchanged. They ran on 18 October 2026 from 16:30 to 16:39 on 5 workers from a clean tree at `55da26a`, one RcppParallel and one BLAS thread per fit (`VECLIB_MAXIMUM_THREADS=1`), R 4.6.1 with Accelerate's `libBLAS.dylib`, mains power with Low Power Mode off (`s10-provenance.txt`). No fit failed. The comparison is the step-off refits of the same replicates at the same settings in the 100-species check's `s10` step (`../b0-community-s100-20261018`), on the same datasets and chain seeds.

## `mu_b0`, with the step off and on

Rank-normalised R-hat, then bulk ESS of 10,000 draws, step off then step on, for replicates 1 to 5 (`before-after.csv`):

-   Ordinary community: R-hat 1.01, 1.02, 1.00, 1.03, 1.02 against 1.00, 1.00, 1.00, 1.00, 1.01; bulk ESS 115, 183, 366, 167, 141 against 258, 530, 427, 331, 447. Median 167 against 427.
-   Rare-dominated community: R-hat 1.04, 1.03, 1.01, 1.03, 1.02 against 1.00, 1.00, 1.01, 1.00, 1.00; bulk ESS 110, 78, 71, 142, 199 against 235, 203, 128, 177, 411. Median 110 against 203.
-   The species' mean B0, the slow direction: bulk ESS 50 to 137 without the step and 98 to 242 with it.

By the rules (`rule.csv`):

-   R1, R-hat under 1.05 and bulk ESS of at least 100 in every fit (strict) or in the median (median): without the step the ordinary arm passes both readings and the rare arm only the median; with the step both arms pass both.
-   R2: the median bulk ESS rises from 167 to 427 and from 110 to 203. Passes in both arms.
-   R3: in every fit the step-on and step-off posterior means differ by at most 1.47 combined Monte Carlo standard errors for `mu_b0` and 1.43 for the mean B0, against a limit of 3, and by at most 0.16 posterior sds. Passes in both arms.

## Descriptive

-   Community occupancy error, mean over the five fits: +0.014 without the step and +0.016 with it (ordinary), -0.012 and -0.012 (rare). B0 coverage over the 50 species-fits: 1.00 and 1.00, 0.94 and 0.96.
-   The step accepts 40 to 49% of its proposals (target 0.44), with adapted proposal sds of 0.27 to 0.46, against `mu_b0` posterior sds of 0.42 to 0.59.
-   Cost: 3.84 to 3.96 minutes per fit with the step against 3.64 to 3.85 without it, 2 to 6% more.

## What it means

The step was built for the stuck chains at 100 species, and this run shows it also helps where the paper's simulations live: at 10 species it roughly doubles `mu_b0`'s effective sample size and lifts the rare-dominated community clear of the convergence bar, with no change to the estimates. With the 100-species rerun ([mu-b0-shift-20261018](../mu-b0-shift-20261018/REPORT.md)), that supports adding it to PR #36 together with the community mean.

## Files

-   `README.md` (design and rules, written before any fit), `launch.sh`, `analyse.R`.
-   `fits.csv`, `before-after.csv`, `rule.csv`, written by `analyse.R`. The draws (`draws/`) are not committed; `launch.sh` regenerates them.
