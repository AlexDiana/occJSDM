# Count model stage 1: what the sampler costs and how well the size is recovered

Date: 25 November 2026. Written by Claude from one run of `dev/count-model/stage1-checks.R` on the branch `claude/count-model-stage1` at `bf18dd3`. The design is `dev/superpowers/specs/2026-11-12-count-model-design.md`; this report answers the two measurements its stage 1 asks for (question 9 and the recovery of the size for rare species). Nothing here is a fit: the pieces were run alone, with the log mean known.

## Recommendation

- **No saddlepoint sampler is needed.** On the realistic matrix the route for shapes from 13 to 170 took 22% of the count sampler's draw time, under the quarter that the stage 1 plan set as the line. It is close to the line, and it would cross it in a survey with larger counts.
- **The cost to look at first is the fractional part of the shape, not the shapes from 13 to 170.** Cells with a shape up to 13 are 92% of the matrix and took 78% of the draw time, and most of what each of them costs is the 20-term series for the fractional part. See *Time per draw by shape*.
- **The size of a species with fewer than 10 non-zero counts is not recovered.** Its posterior is close to the prior. See *Size recovery by rarity*.

## How it was run

- Command: `Rscript dev/count-model/stage1-checks.R`, once, in the foreground, in one process with no parallel workers. Seeds: `set.seed(20261125)` and `setOccJSDMSeed(20261125L)`.
- Machine: Apple M4 laptop, 10 cores. R version 4.6.1 (2026-06-24).
- Power state: on battery ("Now drawing from 'Battery Power'"), `lowpowermode 0`.
- One other R job was running on the laptop at the time, a single-threaded script from another session (`slab_prior.R`) using one core. The timings below may be a little slower than on an idle machine. The shares of time, which the recommendation rests on, compare routes timed within a few seconds of each other.

## Time per draw by shape

Seconds for 200,000 draws of `samplePGcount_cpp()`, and microseconds per draw.

- Shape 0.5: 0.655 s at predictor 0 (3.3 microseconds per draw); 0.772 s at predictor -3 (3.9).
- Shape 1.5: 0.697 s (3.5); 0.871 s (4.4).
- Shape 5.5: 0.903 s (4.5); 1.057 s (5.3).
- Shape 13.5: 1.335 s (6.7); 1.606 s (8.0).
- Shape 50.5: 3.139 s (15.7); 4.322 s (21.6).
- Shape 120.5: 6.529 s (32.6); 9.584 s (47.9).
- Shape 169.5: 10.445 s (52.2); 12.239 s (61.2).
- Shape 171, the normal route: 0.017 s (0.09); 0.019 s (0.10).

What the numbers show:

- **Each whole unit of shape costs about 0.26 microseconds at predictor 0.** That is the slope from shape 1.5 to 5.5, and it is one draw of the existing unit sampler.
- **The fractional part costs about 3.2 microseconds, as much as 12 unit draws.** A draw at shape 0.5 has no whole-number part and takes 3.3 microseconds. Every cell whose shape is not a whole number pays this, which is every cell unless the size is a whole number.
- **Cost rises in proportion to the shape up to 170 and then drops by a factor of about 600** when the normal route takes over.

## Where the time goes on a realistic matrix

The matrix: 250 sites by 300 species, one sample per site at unit effort. Species' centres are normal with mean -1 and standard deviation 2; each site's log mean is normal around its species' centre with standard deviation 1.5; sizes are log-normal with log mean 0 and log standard deviation 1. The shape of each cell is its count plus the species' size, and the predictor is the log mean minus the log size.

- 63% of the counts are zero. The largest count is 14,617.
- The shape has median 1.72, 90th percentile 9.63, 99th percentile 121 and 99.9th percentile 1,023.
- The predictor runs from -10.4 to 10.8, so no cell is beyond 12.

Share of cells and share of draw time by route, each route timed by drawing its cells alone 20 times:

- **Extreme predictor (beyond 12 in absolute value):** 0 cells, no time.
- **Shape up to 13:** 69,360 cells (92.5%); 5.709 s; 77.9% of the time.
- **Shape from 13 to 170:** 5,100 cells (6.8%); 1.617 s; 22.1% of the time.
- **Shape above 170:** 540 cells (0.7%); 0.002 s; 0.03% of the time.

One draw of the whole matrix took 0.40 s (the mean of 20), which is 5.3 microseconds per cell. At that rate the Polya-Gamma draws of a count fit of this size would take about 66 minutes per 10,000 iterations on this laptop, before anything else in the sampler runs.

Applying the plan's rule: the route for shapes from 13 to 170 takes under a quarter of the draw time, so no saddlepoint sampler is needed. Two cautions on that reading:

- 22% is close to 25%. The share depends on how many counts fall between 13 and 170, so a survey with larger counts could cross the line. If it does, the first thing to try is the 20-term series for the whole shape, which is already built. It costs about 3.2 microseconds whatever the shape, against 6.7 at a shape of 13.5 and 52 at a shape of 169.5 for the exact route, and it is an approximation where the exact route is exact. A saddlepoint sampler would come only after that.
- The larger share, 78%, is in cells with a shape up to 13, and there the fractional part dominates. This is not what question 9 asked about, and nothing was changed in response. It is recorded here because it is where a speed-up would pay most: for example fewer terms in the series, which trades a little variance for time. That is a decision for Doug and Alex.

## Size recovery by rarity

The data: 60 species by 250 sites, one sample per site at unit effort. Species' centres are evenly spaced from -4 to 3, each site's log mean is normal around its species' centre with standard deviation 1 (the plan does not give this spread; 1 is what the unit test uses), and the true sizes alternate between 0.5 and 5. The log mean is known to the update. 5,000 iterations of `count_size_step()` from a starting size of 1, the first 2,000 adapting the step, with the prior standard deviation at 1. Summaries are from the 3,000 kept draws. Acceptance rates over the kept draws were between 0.38 and 0.50 for every species, against a target of 0.44.

For reference, the prior alone (a half-normal with standard deviation 1 on the inverse square root of the size) has a median size of about 2.2 and a 95% interval from about 0.2 to about 1,000.

By number of non-zero counts, all species:

- **Under 10 (5 species):** median error of the posterior median on the log scale +0.70; 95% intervals cover the truth for 5 of 5; the median ratio of the interval's upper to lower limit is about 2,900.
- **10 to 49 (17 species):** median error +0.29; cover 17 of 17; median ratio about 780.
- **50 to 149 (17 species):** median error +0.04; cover 17 of 17; median ratio 2.6.
- **150 and over (21 species):** median error -0.01; cover 20 of 21 (95%); median ratio 1.6.

The same, split by the true size, because the two behave differently:

- **True size 0.5.** Under 10 non-zero counts (3 species): median error +0.88, ratio about 2,900. 10 to 49 (10 species): +0.33, ratio 5.8. 50 to 149 (9 species): -0.004, ratio 1.9. 150 and over (8 species): +0.06, ratio 1.5, cover 7 of 8.
- **True size 5.** Under 10 (2 species): median error -0.92, ratio about 4,900. 10 to 49 (7 species): +0.13, ratio about 3,300. 50 to 149 (8 species): +0.13, ratio 106. 150 and over (13 species): -0.08, ratio 2.0, cover 13 of 13.

What the numbers show:

- **With fewer than 10 non-zero counts the size is not recovered.** The posterior median was between 1.0 and 2.0 for all five species, whether the true size was 0.5 or 5, and the 95% intervals ran from a lower limit of 0.13 to 0.25 up to an upper limit of 300 to 1,500. That is the prior, barely changed. The intervals cover the truth only because they are so wide. The update did not fail for any species, including the one with a single non-zero count.
- **A large size needs far more data than a small one.** At a true size of 0.5 the interval is already narrow (a ratio under 6) with 10 to 49 non-zero counts. At a true size of 5 it is still over three orders of magnitude wide there, and still about 100-fold wide with 50 to 149 non-zero counts. With a size of 5 the counts of a rare species are close to Poisson, and the data cannot tell a size of 5 from a much larger one, so the upper limit is set by the prior.
- **With 150 or more non-zero counts the size is recovered.** The median error is within 8% at both true sizes, and the interval is 1.6-fold wide in the median.
- **For the first use this matters.** The beetle survey has several hundred species, many seen a handful of times. For those the size will be set by its prior, which pulls towards the Poisson as intended. Whether that is acceptable, or whether the sizes should share a community prior, is question 5 in the design, still open for Alex.

## What this run does not show

- It is one simulated dataset for each measurement, with one seed. The group summaries rest on 5 to 21 species each.
- The log mean was known. In a fit it is estimated together with the size, and latent factors also create overdispersion, so recovery in a fit will be worse than here. Stage 2 measures that.
- Timings are from one laptop on battery with one other job running. They say nothing about a Linux server.
