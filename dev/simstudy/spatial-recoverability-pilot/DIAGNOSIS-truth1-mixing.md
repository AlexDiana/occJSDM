# Diagnosis: the truth 1 oracle fails its convergence screen

Written 10 October 2026, after the pilot's fit and score steps finished and before any report. Status: diagnosis complete, no fix applied, nothing committed. This note is written for an independent review. It separates what was measured from what is inferred, and lists what has not been checked.

## Summary

The pilot ran to completion under the frozen design: 240 jobs fitted, no job errors, `verify.R` passing all 17 checks. The map reading rule returned `stop`, but the cause is not a finding about recoverability. All 180 truth 1 jobs failed the MCMC convergence screen even after the prescribed extension, so all nine truth 1 cells are `unresolved` under the spec, which makes them neither pass nor fail. The `stop` branch fired only because no cell passed.

The failure is a mixing problem in the truth 1 joint elliptical-slice sampler, and it worsens with site count. A blockwise variant of the same sampler, with one ellipse angle per block, reached about ten times the median effective sample size at equal iterations on two test communities, and agreed with the joint sampler on the posterior where both could be compared. The evidence is strong enough to support a sampler change, and thin enough that it should be reviewed before one is made.

## What was run

The study directory is `dev/simstudy/results/spatial-recoverability-20261010/`. Truth 2 launched at 15:13 BST and completed at 16:22. Truth 1 launched at 16:23 and completed at 20:19. Scoring completed at 20:31. The timing probe had projected 36.90 CPU-hours at the upper bound, 4.61 wall-clock hours across eight workers, which passed the fixed 12-hour rule.

The budget and screen are the frozen ones from `config.rds`: four chains, 3,000 burn-in and 4,000 retained iterations, with one extension to 8,000 burn-in and 16,000 retained iterations if the screen fails. The screen requires rank-normalised Rhat at most 1.01 and bulk and tail ESS at least 400 for three quantities per species (field average, field-truth projection, probability average). For a truth 1 joint job that is 36 quantities, and a species is `unresolved` if any of its three fails. The gate marks a cell `unresolved` if any of its ten communities has any unresolved species among those it uses.

## What was measured

All numbers in this section come from `scores/screen.csv`, `scores/jobs.csv` and `scores/species.csv` in the study directory, and from the stored job files. They are reproduced by the scripts in `diagnosis/`.

- All 180 truth 1 jobs used the extension and all 180 ended `unresolved`. Of 2,160 species-level screen results for truth 1, 160 passed.
- Truth 2 is mostly resolved. Of 720 species-level results, 629 passed. The 91 failures are in the n=300 cell (10) and the n=900 cell (81). The n=100 cell passed every screen result.
- Truth 1 mixing degrades with site count and does not depend on range. Median bulk ESS across the 36 quantities, by cell, for range 0.03, 0.08 and 0.15: at n=100 it is 389, 451 and 470; at n=300 it is 147, 153 and 157; at n=900 it is 49, 45 and 66. Median Rhat rises from about 1.008 at n=100 to about 1.025 at n=300 and 1.06 to 1.08 at n=900. The Rhat condition was the less binding one for n=100 and the ESS condition was binding everywhere above it.
- Truth 2, which uses the sweep's unmodified per-species sampler, has median bulk ESS of 1,413, 1,256 and 1,047 at n=100, 300 and 900.
- Mean proposals per slice step for truth 1 rise from 6.7 at n=100 to 7.8 at n=300 and 8.9 at n=900. For truth 2 they are 4.1, 5.2 and 6.3. More proposals per step means the accepted angle is contracting.
- Per block, in the stored thinned draws of three example jobs (range 0.08, replicate 1, uniform, at n=100, 300 and 900), the weakly informed independent blocks are the slowest. At n=900 their ESS for the block mean, block SD and projection on the block's posterior mean is mostly 12 to 65, against about 70 to 150 for the two shared fields. The per-block ESS values at n=900 come from 400 thinned draws per chain, so they are noisy.
- The extended chains are not stationary at n=900 in the example job. The difference between the first-half and second-half means of the 36 traces reaches 0.036 against a median trace SD of 0.052. The integrated autocorrelation time of the traces is about 160 iterations at n=100, 430 at n=300 and 1,300 at n=900 (medians), so 16,000 retained iterations per chain are not enough at the largest size.

## Hypothesis

The current sampler draws one ellipse angle for the whole 14-block state (two shared fields and twelve independent parts). The shared fields are informed by all twelve species' observations and are pinned tightly relative to their prior, which forces a small accepted angle. The independent parts are weakly informed and need larger moves, but they are updated with the same small angle, so they diffuse slowly. The effect grows with site count because the likelihood becomes more informative about every block.

This hypothesis fits the per-block ESS pattern and the rise in proposals per step. The accepted angles were not recorded per block, so the mechanism is inferred, not observed.

## Test

A diagnostic blockwise elliptical slice sampler was written in `diagnosis/block-ellipse.cpp`. It cycles over the 14 blocks and gives each its own angle, updating one block at a time conditional on the others, with the same prior factor and likelihood. It does not touch the study code.

The comparison script `diagnosis/compare.R` runs both samplers on the same community and starting values, at the frozen initial budget (3,000 burn-in, 4,000 retained, four chains), and screens the 36 traces with the study's own thresholds.

- Truth 1, range 0.08, n=300, replicate 1, uniform: the joint sampler has median ESS 30.1 (minimum 16.1) and median Rhat 1.107 (maximum 1.181), in 20 seconds. The blockwise sampler has median ESS 370.5 (minimum 141.8) and median Rhat 1.009 (maximum 1.036), in 30 seconds. ESS per second is 1.47 and 12.29.
- Truth 1, range 0.08, n=900, replicate 1, uniform: the joint sampler has median ESS 12.7 (minimum 6.9) and median Rhat 1.250 (maximum 1.587), in 120 seconds. The blockwise sampler has median ESS 402.8 (minimum 114.2) and median Rhat 1.008 (maximum 1.032), in 161 seconds. ESS per second is 0.11 and 2.5. The raw log is `diagnosis/compare-n900.log`.
- Fraction of the 36 traces passing the screen at the initial budget: 0 for the joint sampler at both sizes, and 0.17 at n=300 and 0.44 at n=900 for the blockwise sampler.

Agreement of targets was checked in `diagnosis/agree.R` on truth 1, range 0.08, n=100, replicate 1, uniform, with 40,000 retained iterations per chain for each sampler. The site-field posterior means correlate at 0.9915 with a maximum difference of 0.163 against a median posterior SD of 0.978. The posterior SD ratio, blockwise over joint, has median 1.002. The 36 trace means agree to within 0.0101, with correlation 0.99997.

## Limits of the evidence

- The comparison uses one community per site count and one range. It has not been run across the grid.
- At the initial budget the blockwise sampler still fails the all-36 screen at both sizes, with minimum ESS 142 and 114. Passing would need the extension. A pass at the extended budget is expected, since extended draws are four times as many, but has not been run.
- The agreement check is at n=100 only and compares two samplers that are both Monte Carlo estimates, so a 0.163 maximum difference is consistent with Monte Carlo error but does not prove equality. It does not test n=300 or n=900.
- The blockwise sampler has no unit tests. In particular, it has not been checked against the one-field, no-independent-part reduction to the sweep's sampler that the spec requires of the joint sampler, or against a deterministic quadrature case.
- The mechanism (one angle across blocks of very different information) is a hypothesis that fits the data. An alternative contribution, for example from the 90 and 10 split between shared and independent variance, was not isolated.
- The 91 failed truth 2 screen results, mostly at n=900, are not diagnosed here. They come from the per-species sampler, which is the sweep's unmodified code, and may have a different cause.

## Options

The sampler swap changes frozen code and the frozen budget path, so it is a decision for Doug and not part of the pilot as handed over.

- Accept the result as it stands: truth 1 is `unresolved` at this budget. The report would say so and would not read `stop` as a finding about recoverability.
- Swap the truth 1 sampler to the blockwise one, with tests that include the spec's one-field reduction and a quadrature check, a revision record in the spec, and the screen thresholds unchanged. Move the 180 existing truth 1 job files aside as a record and re-run truth 1 only. A reasonable estimate is four to five hours, as in the first run. Then score, verify and report.
- Keep the current sampler and raise the budget. The integrated autocorrelation time of about 1,300 iterations at n=900 suggests this would need far more than the prescribed extension, and it would also change the frozen budget.

## Review questions

- Is the one-angle explanation the right reading of the per-block ESS pattern, or does the evidence point to another cause?
- Is the blockwise sampler's target the same as the joint sampler's? Check the update in `block-ellipse.cpp`, in particular the incremental update of the per-species effects and the likelihood.
- Does the evidence justify changing the sampler for the full truth 1 run, or is more testing needed first, and if so which?
- Are the truth 2 failures at n=900 likely to share this cause?

## Files

- `diagnosis/screen-by-cell.R`: screen results by cell, quantity and arrangement, and mean proposals.
- `diagnosis/per-block-mixing.R`: trace and per-block mixing for three example jobs.
- `diagnosis/block-ellipse.cpp`: the diagnostic blockwise sampler.
- `diagnosis/compare.R`: joint against blockwise on one community. Arguments are the job name, burn-in and retained iterations.
- `diagnosis/agree.R`: the posterior agreement check at n=100.
- `diagnosis/compare-n900.log`: the raw n=900 comparison output.

These scripts read the study directory and the stored job files under `dev/simstudy/results/spatial-recoverability-20261010/`, which are not committed. Run them through `run-r`, from the worktree root.
