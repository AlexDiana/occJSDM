# The field-sample ladder on current code: when does the prior stop deciding occupancy?

The R = 300 paper study found each species' occupancy squeezed towards the middle in the eDNA cells, with rare species estimated too common and common species too rare. The squeeze comes from the occupancy intercept `B0`, which the data barely separate from the rate at which samples catch a species' DNA, so its prior decides the split. This run repeats the 29 July field-sample ladder on current code to find how much replication undoes the squeeze, and of which kind. The design and the launch are in [README.md](README.md).

## Summary

-   **Field samples undo the squeeze; PCR replicates do not.** The slope of each species' estimated occupancy on its true occupancy is 0.44 with 2 field samples per site, 0.66 with 5, 0.86 with 10 and 0.95 with 20, but only 0.52 with 2 field samples of 30 PCR replicates each. Species with true occupancy under 0.2, averaging 0.16, are estimated at 0.35, 0.30, 0.23 and 0.16 as field samples rise, and at 0.33 with 30 PCR replicates.
-   **In this design the squeeze is mostly gone by 10 field samples per site and essentially gone by 20.** The mean absolute error of a species' occupancy falls from 0.085 to 0.065, 0.040 and 0.020. With 30 PCR replicates it is 0.077.
-   **The trade-off behind the squeeze weakens with field samples only.** The correlation between a species' `B0` error and its collection-intercept error is -0.36, -0.35, -0.28 and -0.16 as field samples rise, and -0.42 with 30 PCR replicates.
-   **PCR replicates fix the lab detection rate instead.** The overestimate of `p` falls from +0.038 at 2 samples to +0.004 at 20, and to +0.002 with 30 PCR replicates.
-   **As the data take over, `B0`'s intervals become too narrow.** `B0` coverage is 0.942, 0.916, 0.889 and 0.873 as field samples rise, and 0.937 with 30 PCR replicates. The intervals narrow from a mean width of 3.11 to 1.40 faster than the bias at the extremes shrinks: at 20 samples, true `B0` below -1.5 is still overestimated by 0.50 on the logit scale and covered 69% of the time. This is the pattern the paper study found in the `binary` cell, where presence is observed directly.
-   **B's coverage does not move.** It is 0.89 to 0.91 in every arm, although its estimates track the truth better with more samples, with slopes rising from 0.68 to 0.83. That fits the non-spatial candidates report, which put B's shortfall on its hierarchical prior and on short chains.

## Run

All 500 fits ran on 16 October 2026 from 01:01 to 03:18 at revision `49b51c1`, with a clean tree at start and at write and no failures: 5 arms of 100 replicates, 5 replicate processes, one thread per fit, on mains power with Low Power Mode off, under R 4.6.1 with Accelerate BLAS. It was launched only after the rare-dominated community run and the non-spatial candidates run had finished, as Doug asked.

## By arm

Coverage standard errors are computed across replicates. Occupancy is each species' estimated occupancy against its true occupancy, from the results file's occupancy rows.

-   `M2`, 2 field samples per site: occupancy slope 0.44, mean absolute error 0.085, `B0` coverage 0.942 (SE 0.007), `B0` slope 0.35, trade-off correlation -0.36, `p` error +0.038.
-   `M5`: occupancy slope 0.66, error 0.065, `B0` coverage 0.916 (0.009), `B0` slope 0.54, trade-off -0.35, `p` error +0.015.
-   `M10`: occupancy slope 0.86, error 0.040, `B0` coverage 0.889 (0.010), `B0` slope 0.68, trade-off -0.28, `p` error +0.007.
-   `M20`: occupancy slope 0.95, error 0.020, `B0` coverage 0.873 (0.012), `B0` slope 0.75, trade-off -0.16, `p` error +0.004.
-   `K30`, 2 field samples with 30 PCR replicates: occupancy slope 0.52, error 0.077, `B0` coverage 0.937 (0.008), `B0` slope 0.41, trade-off -0.42, `p` error +0.002.

`theta0` coverage moves from 0.974 at 2 samples to 0.948 at 20, and is 0.985 with 30 PCR replicates. `q` covers 0.94 to 0.96 throughout.

## Species occupancy by true occupancy

Estimated occupancy, averaged within each band of true occupancy, for `M2`, `M5`, `M10`, `M20` and `K30` in that order:

-   Under 0.2, true mean 0.16 (19 of 1,000 species-fits per arm): 0.35, 0.30, 0.23, 0.16, 0.33.
-   0.2 to 0.3, true mean 0.26 (61): 0.38, 0.35, 0.31, 0.28, 0.38.
-   0.3 to 0.5, true mean 0.42 (420): 0.44, 0.44, 0.43, 0.42, 0.45.
-   0.5 to 0.7, true mean 0.59 (427): 0.52, 0.56, 0.58, 0.58, 0.54.
-   0.7 to 0.8, true mean 0.75 (58): 0.59, 0.68, 0.73, 0.73, 0.62.
-   Over 0.8, true mean 0.84 (15): 0.66, 0.70, 0.81, 0.80, 0.72.

## B0 by true value

Mean error on the logit scale, with coverage, for the most extreme species, in the same arm order:

-   True `B0` below -1.5 (65 per arm): +1.19 (0.82), +0.91 (0.71), +0.67 (0.77), +0.50 (0.69), +1.16 (0.77).
-   True `B0` above 1.5 (68 per arm): -1.39 (0.59), -0.96 (0.71), -0.60 (0.79), -0.53 (0.74), -1.19 (0.65).

The error at the extremes shrinks with field samples, but the intervals shrink faster, so coverage there stays poor. The full table is [b0-by-band.csv](b0-by-band.csv).

## Against the 29 July ladder

Replicates 1 to 50 share their truths with the July ladder; the largest difference is 9e-16, which is rounding. On the same data, current code gives `B0` slopes a little closer to 1 than July's in every arm, but at 5, 10 and 20 field samples its `B0` intervals are 18, 25 and 28% narrower, and coverage is 3.8 (SE 1.0), 7.2 (SE 1.2) and 5.8 (SE 1.2) points lower. At 2 samples and with 30 PCR replicates the coverage changes are within two standard errors. The calibration report traced the same narrowing in `B` and in `binary`'s `B0` to the spatial commit `f1d2e12`, which this run did not isolate. The paired table is [b0-paired-july.csv](b0-paired-july.csv).

## What this does not show

-   **Few rare species.** The truths come from the model's own prior, so only 19 of 1,000 species-fits per arm have true occupancy under 0.2. The rare-dominated community run, on branch `claude/rare-dominated-occupancy`, showed that the bias at a given true occupancy is the same in a rare community, but a ladder on a rare community would show directly how many field samples rare species need.
-   **One design.** The answer is for 100 sites, 10 species, 3 PCR replicates per sample, a spatial field and traits, and true sample-capture rates of 0.2 to 0.5 at average conditions. Lower capture rates would probably need more field samples.
-   **Convergence.** It was not checked. The non-spatial candidates run found `B` unconverged at this chain length, so the `B` and `B0` interval results may be partly sampling.

## What this changes

-   **Survey design.** With 2 field samples per site, the estimates for rare and common species mostly reflect the prior rather than the data. If species-level occupancy matters, about 10 field samples per site removed most of the squeeze in this design. Extra PCR replicates sharpen the lab detection rate but do not separate occupancy from capture.
-   **The paper.** Occupancy bias should be reported by true occupancy and by field replication, and the advice above belongs with it.
-   **Calibration.** At high replication the `B0` intervals undercover at the extremes, so more data exposes a calibration problem rather than removing it. That bears on the community-mean intercept question for Alex in `TODO.md`.

## Files and reproduction

-   [analyse.R](analyse.R) writes [ladder-by-arm.csv](ladder-by-arm.csv), [occupancy-by-band.csv](occupancy-by-band.csv), [b0-by-band.csv](b0-by-band.csv), [pairing-across-arms.csv](pairing-across-arms.csv) and [b0-paired-july.csv](b0-paired-july.csv). The July results file it pairs with lives in the main checkout's `dev/simstudy/results/` and is not committed.
-   `results/simstudy-20261016-031826.rds` holds the run; `run.log` and `provenance.txt` record it.

Run from the repository root:

```sh
Rscript dev/simstudy/field-sample-ladder/analyse.R
```
