# Simulation of the sigma_b move, 22 November 2026

This records steps 1 and 2 of the evidence in `dev/superpowers/specs/2026-11-22-sigma-b-move-design.md`, run as Task 7 of `dev/superpowers/plans/2026-11-22-sigma-b-move.md`. Step 3, the replay by the AlphaEco project, was run there on 22 and 23 November 2026 and is summarised at the end, in the section "Step 3".

## What was run

Every fit is on data from `simulateOccJSDMData()`: 256 sites, 48 species, 5 measured covariates, 64 embedding columns, 2 traits, 4 latent traits, model "binary", true sigma_b 0.1, data set k simulated with seed 1000 + k and fitted with seed 2000 + k. Each fit has 4 chains of 5,000 burn-in and 20,000 kept sweeps, with `n_factors = 2` and the embedding columns under the grouped prior. The script is `repro.R`, and `summarise.R` reads its output. The result files are in `dev/simstudy/results/sigma-b-move-20261122/`, which git ignores.

- Without the move, label `base`: library built from a2a64a1, started 09:19 GMT and last file written 12:04 GMT, 3 workers, on AC power with Low Power Mode off at the start.
- With the move, label `move`: library built from de7e92d (this branch's head before this record), started 12:05 GMT and last file written 14:42 GMT, 3 workers, on AC power with Low Power Mode off at the start. No other `Rscript` job was running at launch.

A kept sweep is "low" when sigma_b is under 0.02. R-hat is `coda::gelman.diag` (no autoburnin, univariate) over the measured rows of the coefficient output.

## The two readings

1. Baseline: **reproduced**. At least one chain is low for 500 or more kept sweeps in a row: data set 4, chains 1, 2 and 3, with longest low runs of 700, 636 and 1,250 sweeps. The other four data sets have no low sweeps.
2. With the move: **stays are short**. The longest low stay in any chain is 87 sweeps, under a fifth of the baseline's longest (1,250 / 5 = 250), and no stay exceeds 500 sweeps.

The share of kept sweeps low was expected not to fall. It fell: in data set 4 it is 0.068 with the move against 0.080 without (chains 0.126, 0.002, 0.053, 0.092 against 0.114, 0.048, 0.080, 0.077). In the other four data sets it is 0 in both runs.

## Measured-coefficient R-hat, largest value and number over 1.1

At 5,000 kept sweeps, then at all 20,000.

- Data set 1: base 1.155 (8 over), 1.240 (7 over); move 1.021 (0 over), 1.004 (0 over).
- Data set 2: base 1.022 (0 over), 1.010 (0 over); move 1.017 (0 over), 1.005 (0 over).
- Data set 3: base 1.072 (0 over), 1.016 (0 over); move 1.044 (0 over), 1.015 (0 over).
- Data set 4: base 1.566 (51 over), 1.106 (1 over); move 1.272 (22 over), 1.313 (35 over).
- Data set 5: base 1.225 (9 over), 1.072 (0 over); move 1.261 (20 over), 1.071 (0 over).

In data set 4 the R-hat at 20,000 kept sweeps is worse with the move than without it: 35 coefficients over 1.1 against 1. Its chains disagree on the level of sigma_b: the median is 0.635 in chain 2 and 0.085 to 0.120 in the other three (0.106 to 0.169 in all four without the move).

## Acceptance of the sigma_b move's proposals

Per chain, chains 1 to 4. The adapted proposal's step was adapted during burn-in from 0.5 and then frozen; the large proposal's step is 1.

- Kept sweeps, adapted proposal: d1 0.57/0.38/0.45/0.46; d2 0.47/0.43/0.45/0.41; d3 0.46/0.43/0.45/0.44; d4 0.81/0.44/0.78/0.87; d5 0.45/0.28/0.41/0.39.
- Kept sweeps, large proposal: d1 0.02/0.02/0.02/0.02; d2 0.02/0.02/0.02/0.02; d3 0.01/0.01/0.01/0.01; d4 0.22/0.03/0.15/0.21; d5 0.01/0.02/0.01/0.01.
- Burn-in, adapted proposal: 0.43 to 0.45 in every chain of every data set.
- Burn-in, large proposal: 0.010 to 0.030 in every chain of every data set.
- Adapted step at the end of burn-in: d1 0.029/0.054/0.042/0.040; d2 0.041/0.046/0.044/0.048; d3 0.019/0.021/0.021/0.020; d4 0.047/0.045/0.044/0.044; d5 0.025/0.056/0.029/0.031. None is near the adaptation bound of 5.

Nothing is claimed beyond the two readings.

## Step 3: the replay by the AlphaEco project

Run by the AlphaEco project on its own data design, 22 November 2026 16:59 to 23 November 2026 01:26, on Doug's laptop, on AC power with Low Power Mode off, at most 5 workers. Its report is held by the design session and is not in this repository, since it was made on that project's data. What follows is a summary: counts only, no dataset content. The design session recounted every count below from the report's per-fit sections.

- Data sets: the 18 simulated data sets whose grouped fits failed the measured-coefficient check in that project's main run, and 6 controls whose fits passed it. Same simulated species, sites and splits as the main run.
- Libraries: this branch at eb362e1 ("move"), and f8a1a21, an earlier commit of PR 45 whose sampler is the same as at a2a64a1 ("no move"). The move ran in every fit of the first and in none of the second.
- Seeds: a new fit seed for every data set, the same in both arms. The project checked that none of them, and none of the chain seeds drawn from them, had been used in any earlier run.
- Fits: 4 chains of 5,000 burn-in and 5,000 kept sweeps; four data sets with the move run on to 20,000 kept. 48 fits, 0 errors.

The reading fixed in the spec, on the 18 data sets at 5,000 kept sweeps with the move: **the move is enough**. 1 of the 18 fails the measured-coefficient check (at most 2 allowed), and no chain of any fit has a run of more than 500 sweeps with sigma_b under 0.02 (the longest is 79).

Beside it, with no label:

- Without the move, on the same seeds: 5 of the 18 fail, and 5 fits have a chain with a run of more than 500 sweeps under 0.02 (the longest 1,459).
- Paired by data set: 18 pass in both arms, 5 fail only without the move, 1 fails only with it.
- Measured coefficients over 1.1, summed over all 24 data sets: 1 with the move, 18 without. Fits with an R-hat of log sigma_b over 1.1: 5 with the move, 11 without.
- Controls: 0 of 6 fail in either arm.
- The four long fits with the move pass at 20,000 kept sweeps, with no measured coefficient over 1.1 (largest R-hat 1.005 to 1.011).
- The one failure with the move has nothing to do with sigma_b: no chain goes low, and one of its four chains sits in the shared-structure state while the others sit in the leftover state (R-hat of log sigma_emb 2.7). The move does not address that. The same data set and seed pass without the move.
- The chains at different levels of sigma_b seen in data set 4 above did not appear: within a fit, the largest chain median of sigma_b is at most 2.16 times the smallest with the move. Between data sets the level runs from 0.05 to 0.31.
- Acceptance over kept sweeps, per chain: adapted proposal 0.23 to 0.64, large proposal 0.04 to 0.29.
- Share of kept sweeps with sigma_b under 0.02, mean over fits: 0.040 with the move, 0.031 without. As in steps 1 and 2, the move shortens the stays and does not empty the low state.
- Time per short fit: 29 to 52 minutes with the move and 34 to 43 without, with up to 5 fits running at once. No cost is visible at this level of noise.

How strong this is: 5 failures against 1 out of 18 could arise by chance about one time in five, so the count alone is suggestive. The disappearance of the long stays (no fit with the move, 5 without) is the direct evidence for the move.

## Summary without the move, in full

```
base d1 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6249
base d1 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6234
base d1 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6269
base d1 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.8965
base d1 fit: R-hat measured max 1.155 (8 over 1.1) at 5,000 kept; 1.240 (7 over 1.1) at 20000 kept; 92.2 minutes; low share 0.000
base d2 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6324
base d2 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6387
base d2 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6326
base d2 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6349
base d2 fit: R-hat measured max 1.022 (0 over 1.1) at 5,000 kept; 1.010 (0 over 1.1) at 20000 kept; 92.3 minutes; low share 0.000
base d3 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4550
base d3 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4579
base d3 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4585
base d3 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4599
base d3 fit: R-hat measured max 1.072 (0 over 1.1) at 5,000 kept; 1.016 (0 over 1.1) at 20000 kept; 92.1 minutes; low share 0.000
base d4 chain 1: low share 0.114, longest low run 700, runs >= 500: 2, median sigma_b 0.1056
base d4 chain 2: low share 0.048, longest low run 636, runs >= 500: 1, median sigma_b 0.1056
base d4 chain 3: low share 0.080, longest low run 1250, runs >= 500: 1, median sigma_b 0.1687
base d4 chain 4: low share 0.077, longest low run 336, runs >= 500: 0, median sigma_b 0.1263
base d4 fit: R-hat measured max 1.566 (51 over 1.1) at 5,000 kept; 1.106 (1 over 1.1) at 20000 kept; 73.4 minutes; low share 0.080
base d5 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0759
base d5 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0758
base d5 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0766
base d5 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0804
base d5 fit: R-hat measured max 1.225 (9 over 1.1) at 5,000 kept; 1.072 (0 over 1.1) at 20000 kept; 73.2 minutes; low share 0.000
base: longest low run over all chains 1250
```

## Summary with the move, in full

```
move d1 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6251
move d1 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6247
move d1 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6221
move d1 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6250
move d1 fit: R-hat measured max 1.021 (0 over 1.1) at 5,000 kept; 1.004 (0 over 1.1) at 20000 kept; 81.5 minutes; low share 0.000, sigma_b acceptance adapted 0.57/0.38/0.45/0.46, large 0.02/0.02/0.02/0.02, adapted step 0.029/0.054/0.042/0.040
move d2 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6288
move d2 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6324
move d2 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6336
move d2 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 0.6359
move d2 fit: R-hat measured max 1.017 (0 over 1.1) at 5,000 kept; 1.005 (0 over 1.1) at 20000 kept; 81.5 minutes; low share 0.000, sigma_b acceptance adapted 0.47/0.43/0.45/0.41, large 0.02/0.02/0.02/0.02, adapted step 0.041/0.046/0.044/0.048
move d3 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4600
move d3 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4622
move d3 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4529
move d3 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.4590
move d3 fit: R-hat measured max 1.044 (0 over 1.1) at 5,000 kept; 1.015 (0 over 1.1) at 20000 kept; 81.4 minutes; low share 0.000, sigma_b acceptance adapted 0.46/0.43/0.45/0.44, large 0.01/0.01/0.01/0.01, adapted step 0.019/0.021/0.021/0.020
move d4 chain 1: low share 0.126, longest low run 87, runs >= 500: 0, median sigma_b 0.0931
move d4 chain 2: low share 0.002, longest low run 12, runs >= 500: 0, median sigma_b 0.6349
move d4 chain 3: low share 0.053, longest low run 36, runs >= 500: 0, median sigma_b 0.1203
move d4 chain 4: low share 0.092, longest low run 52, runs >= 500: 0, median sigma_b 0.0848
move d4 fit: R-hat measured max 1.272 (22 over 1.1) at 5,000 kept; 1.313 (35 over 1.1) at 20000 kept; 75.8 minutes; low share 0.068, sigma_b acceptance adapted 0.81/0.44/0.78/0.87, large 0.22/0.03/0.15/0.21, adapted step 0.047/0.045/0.044/0.044
move d5 chain 1: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0876
move d5 chain 2: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0584
move d5 chain 3: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0854
move d5 chain 4: low share 0.000, longest low run 0, runs >= 500: 0, median sigma_b 1.0826
move d5 fit: R-hat measured max 1.261 (20 over 1.1) at 5,000 kept; 1.071 (0 over 1.1) at 20000 kept; 75.7 minutes; low share 0.000, sigma_b acceptance adapted 0.45/0.28/0.41/0.39, large 0.01/0.02/0.01/0.01, adapted step 0.025/0.056/0.029/0.031
move: longest low run over all chains 87
```

## Fewer latent traits, 26 November 2026

On PR 46 (25 November 2026) Alex Diana explained the low state as the latent-trait term over-fitting the measured block: with 4 latent traits and 5 measured covariates the model cannot tell the leftover from the latent-trait term, and the leftover's scale goes to zero. He suggested fewer latent traits or more covariates. Doug asked for a test of the first.

What was run:

- `lattrait.R`, with the reading rules fixed in its header before any run, and `lattrait_summarise.R`. The script was not edited after the runs began.
- The five data sets and seeds of `repro.R`, a library without the move built from a2a64a1, 4 chains of 5,000 + 20,000, 5 workers. Arm `gt2` fits 2 latent traits (09:09 to 11:03 GMT, 113 to 114 minutes per fit). Arm `gt4` fits 4, as a control (11:05 to 13:21 GMT, 134 to 136 minutes per fit).
- Power: the laptop went onto battery during the first arm and back onto AC power during the second. Low Power Mode was off throughout. The arms ran one after the other, so their times are not a fair comparison.
- The control arm's `sigma_b` draws are identical to those in the `base` files of `repro.R` in all five data sets.
- The measured block is 6 design columns, not 5: one of the 5 measured covariates is categorical with three levels and enters as two indicator columns.

The readings:

1. Low stays: **long stays gone**. With 2 latent traits no chain of any data set has a kept sweep under 0.02. With 4, data set 4 has four stays of 500 sweeps or more (700, 636 and 1,250 in chains 1 to 3) and 8.0% of its kept sweeps are low; the other four data sets have none.
2. Measured coefficients: 2 fits fail with 2 latent traits (data sets 1 and 2) and 3 with 4 (data sets 1, 4 and 5). At all 20,000 kept sweeps the coefficients over 1.1 are 11 and 60 in data sets 1 and 2 with 2 latent traits, and 7 and 1 in data sets 1 and 4 with 4.
3. Recovery, error with 2 latent traits over error with 4. Embedding contribution: 1.032, 0.988, 1.025, 1.045, 1.126, which reads **no harm**. Measured coefficients: 1.057, 1.023, 1.001, 1.243, 1.037, which reads **no harm**. Coverage of the measured coefficients' 95% intervals with 2: 0.944, 0.944, 0.955, 0.878, 0.934; with 4: 0.951, 0.934, 0.948, 0.931, 0.938.

Beside the rules:

- The two failures with 2 latent traits are fits with one chain in the other state of the embedding: median `sigma_emb` 6.3 against 0.57 to 0.66 in data set 1, and 6.1 against 0.30 in data set 2. That chain's `sigma_b` sits at a different level too (0.63 against 1.10, and 0.65 against 1.80). They do not clear by 20,000 kept sweeps.
- Chain medians of `sigma_b` are 0.63 to 1.80 with 2 latent traits. With 4 they are 0.62 to 1.46 in four data sets and 0.11 to 0.17 in data set 4. The simulated truth is 0.1, so the one data set with the low state is also the one whose `sigma_b` is near the truth, and it is the one where 2 latent traits recover the measured coefficients worst (ratio 1.243, coverage 0.878).
- The package's own verdict is "ok" in 2 of 5 fits with 2 latent traits (data sets 3 and 5) and in 3 of 5 with 4 (data sets 2, 3 and 5).

What this does not show:

- The truth has 4 latent traits, so a fit with 2 has fewer than the truth. Part of the rise in `sigma_b` may be the leftover taking up what the missing latent traits carried. A truth with 2 latent traits, fitted with 2 and with 4, would separate the two. It has not been run.
- Only 1 data set of 5 shows the low state with 4 latent traits, so reading 1 rests on that one.

The summary in full:

```
gt2 d1: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 0.628 1.100 1.097 1.098
gt2 d1: measured R-hat max 1.260 (10 over 1.1) at 5,000 kept; 1.262 (11 over 1.1) at all 20000; verdict review; RMSE embedding 0.6911, measured 0.3738; coverage 0.944, width 1.433; median sigma_emb by chain 6.31 0.57 0.61 0.66; 114 minutes
gt2 d2: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 1.798 1.804 1.802 0.651
gt2 d2: measured R-hat max 1.730 (65 over 1.1) at 5,000 kept; 1.715 (60 over 1.1) at all 20000; verdict review; RMSE embedding 0.6458, measured 0.4807; coverage 0.944, width 1.723; median sigma_emb by chain 0.30 0.30 0.31 6.12; 114 minutes
gt2 d3: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 1.664 1.672 1.666 1.672
gt2 d3: measured R-hat max 1.013 (0 over 1.1) at 5,000 kept; 1.002 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.7433, measured 0.4664; coverage 0.955, width 1.672; median sigma_emb by chain 0.64 0.67 0.66 0.66; 114 minutes
gt2 d4: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 0.654 0.655 0.655 0.658
gt2 d4: measured R-hat max 1.007 (0 over 1.1) at 5,000 kept; 1.002 (0 over 1.1) at all 20000; verdict review; RMSE embedding 0.8040, measured 0.4965; coverage 0.878, width 1.434; median sigma_emb by chain 5.78 5.94 6.15 5.82; 114 minutes
gt2 d5: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 1.349 1.359 1.350 1.349
gt2 d5: measured R-hat max 1.009 (0 over 1.1) at 5,000 kept; 1.002 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.7915, measured 0.4587; coverage 0.934, width 1.594; median sigma_emb by chain 0.69 0.66 0.74 0.70; 113 minutes
gt4 d1: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 0.625 0.623 0.627 0.897
gt4 d1: measured R-hat max 1.155 (8 over 1.1) at 5,000 kept; 1.240 (7 over 1.1) at all 20000; verdict review; RMSE embedding 0.6699, measured 0.3538; coverage 0.951, width 1.386; median sigma_emb by chain 0.55 0.53 0.56 0.07; 136 minutes
gt4 d2: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 0.632 0.639 0.633 0.635
gt4 d2: measured R-hat max 1.022 (0 over 1.1) at 5,000 kept; 1.010 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.6540, measured 0.4697; coverage 0.934, width 1.570; median sigma_emb by chain 0.32 0.32 0.32 0.30; 135 minutes
gt4 d3: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 1.455 1.458 1.459 1.460
gt4 d3: measured R-hat max 1.072 (0 over 1.1) at 5,000 kept; 1.016 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.7255, measured 0.4661; coverage 0.948, width 1.683; median sigma_emb by chain 0.17 0.17 0.19 0.17; 135 minutes
gt4 d4: low share by chain 0.114 0.048 0.080 0.077; longest low stay 700 636 1250 336; long stays 4; median sigma_b 0.106 0.106 0.169 0.126
gt4 d4: measured R-hat max 1.566 (51 over 1.1) at 5,000 kept; 1.106 (1 over 1.1) at all 20000; verdict review; RMSE embedding 0.7693, measured 0.3996; coverage 0.931, width 1.387; median sigma_emb by chain 5.38 5.61 4.11 5.28; 134 minutes
gt4 d5: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 1.076 1.076 1.077 1.080
gt4 d5: measured R-hat max 1.225 (9 over 1.1) at 5,000 kept; 1.072 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.7027, measured 0.4421; coverage 0.938, width 1.595; median sigma_emb by chain 0.06 0.04 0.05 0.05; 135 minutes
```

## Fewer latent traits on a truth the fit can match, 26 November 2026

The first test's account of its truth was wrong, and this section corrects it. The simulator draws the embedding's blends from species responses of their own, apart from the latent traits of the measured block, and a fit has one set of latent traits for both. So the default truth (4 latent traits and 3 blends) has 7 latent dimensions. Neither arm of the first test matched it, and the 4-trait arm was not "the truth's number" as the section above says. The first test could therefore not tell "fewer latent traits identify `sigma_b`" from "fewer latent traits leave more for the leftover to carry".

Doug asked for the fairer test. It uses a small truth: 1 latent trait and 1 blend, 2 latent dimensions. A fit with 2 latent traits matches it and a fit with 4 has 2 to spare.

What was run:

- `lattrait.R --truth small`, with its reading rules written into the header before any run, and `lattrait_summarise.R --prefix small-`. The script was not edited after the runs began (its checksum at the end begins 14fe30d3, as at the start).
- Everything else as in the first test: 256 sites, 48 species, 64 embedding columns, `sigma_b` = 0.1, the same data and fit seeds, a library without the move built from a2a64a1, 4 chains of 5,000 + 20,000, 5 workers. The data sets differ from the first test's because the truth differs. Both arms were fitted to the same five data sets (the saved true coefficients are identical between arms).
- Arm `small-gt2` ran 13:45 to 15:48 GMT (121 to 122 minutes per fit) and arm `small-gt4` 15:48 to 17:38 GMT (108 to 109 minutes per fit), on AC power throughout. The arms ran one after the other, so their times are not a fair comparison.
- A check on the saved truth: the first singular value of the true embedding coefficients carries 94% to 98% of their sum of squares across species in all five data sets, against three values sharing it in the first test. So the truth has one blend, as intended.

The readings:

1. Low stays: **does not support the suggestion**. The 2-trait arm has 17 long stays, in data sets 2, 3 and 5 (longest by data set 0, 7,778, 1,034, 179, 1,060). The 4-trait arm has 40, in all five data sets (longest 1,461, 3,000, 2,771, 1,891, 2,403).
2. Measured coefficients: 3 fits fail with 2 latent traits (data sets 2, 4 and 5) and 5 with 4. At all 20,000 kept sweeps the coefficients over 1.1 are 13 in data set 2 with 2 latent traits, and 5 and 10 in data sets 2 and 4 with 4.
3. Recovery, error with 2 latent traits over error with 4. Embedding contribution: 0.965, 0.969, 0.986, 0.972, 0.985, which reads **no difference**. Measured coefficients: 0.879, 0.957, 0.938, 0.921, 0.969, which reads **no difference** (one ratio under 0.91, and the rule needs two).

Beside the rules:

- Share of kept sweeps under 0.02, by data set: 0.000, 0.303, 0.098, 0.010, 0.030 with 2 latent traits (mean 0.088) and 0.093, 0.347, 0.209, 0.321, 0.153 with 4 (mean 0.224). So spare latent traits make the low state more common, in every data set, which is the direction of Alex's account. Matching the truth does not remove it.
- Chain medians of `sigma_b` (truth 0.1). With 2 latent traits: 0.13 to 0.14, 0.01 to 0.05, 0.07 to 0.09, 0.09 to 0.10, 0.11. With 4: 0.06 to 0.08, 0.03, 0.04 to 0.06, 0.03 to 0.04, 0.04 to 0.08. With 2 the level is near the truth in four data sets of five. With 4 it is under the truth in all five.
- In the first test the medians were 0.6 to 1.8 against the same truth of 0.1, in both arms bar one data set. That is the leftover carrying latent dimensions the fit lacked. It is why the low state was rare there, and why 2 latent traits appeared to remove it.
- Data set 2 is the worst in both arms. With 2 latent traits one chain is low for 54% of its kept sweeps, with one stay of 7,778.
- All ten error ratios are under 1, so 2 latent traits recover slightly better on this truth, by 1% to 4% for the embedding and 3% to 12% for the measured coefficients. The rule reads this as no difference.
- Coverage of the measured coefficients' 95% intervals with 2: 0.965, 0.875, 0.920, 0.951, 0.944; with 4: 0.976, 0.924, 0.938, 0.962, 0.958. The intervals are wider with 4 (mean width 0.83 to 0.96 against 0.71 to 0.85).
- The package's own verdict is "ok" in 4 of 5 fits with 2 latent traits (all but data set 2) and in 2 of 5 with 4 (data sets 1 and 5).

What the two tests show together:

- Fewer latent traits do not stop the low state. On a truth the fit matches, the chains still spend 9% of their kept sweeps under 0.02 and have stays of over 1,000 sweeps.
- Spare latent traits make it worse, as Alex said they would.
- The first test's "long stays gone" was produced by under-fitting and is withdrawn as support for the suggestion.

What this does not show:

- Whether the move of PR 46 shortens the stays on this truth. The move was tested only on the default truth, where the low state is rare. These ten fits have no move.
- Anything beyond five data sets of one design.

The summary in full:

```
gt2 d1: low share by chain 0.000 0.000 0.000 0.000; longest low stay 0 0 0 0; long stays 0; median sigma_b 0.136 0.136 0.135 0.134
gt2 d1: measured R-hat max 1.046 (0 over 1.1) at 5,000 kept; 1.017 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.4232, measured 0.1868; coverage 0.965, width 0.766; median sigma_emb by chain 0.16 0.18 0.17 0.16; 122 minutes
gt2 d2: low share by chain 0.178 0.288 0.208 0.538; longest low stay 1170 1430 1205 7778; long stays 9; median sigma_b 0.047 0.041 0.053 0.012
gt2 d2: measured R-hat max 2.514 (222 over 1.1) at 5,000 kept; 1.172 (13 over 1.1) at all 20000; verdict review; RMSE embedding 0.4812, measured 0.2301; coverage 0.875, width 0.706; median sigma_emb by chain 0.53 0.57 0.58 0.62; 121 minutes
gt2 d3: low share by chain 0.033 0.152 0.113 0.094; longest low stay 184 1034 561 658; long stays 6; median sigma_b 0.075 0.077 0.072 0.085
gt2 d3: measured R-hat max 1.071 (0 over 1.1) at 5,000 kept; 1.024 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.3581, measured 0.2364; coverage 0.920, width 0.852; median sigma_emb by chain 0.71 0.72 0.72 0.66; 122 minutes
gt2 d4: low share by chain 0.011 0.006 0.025 0.000; longest low stay 83 30 179 0; long stays 0; median sigma_b 0.094 0.098 0.097 0.104
gt2 d4: measured R-hat max 1.135 (2 over 1.1) at 5,000 kept; 1.019 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.4294, measured 0.1862; coverage 0.951, width 0.761; median sigma_emb by chain 0.15 0.17 0.15 0.17; 121 minutes
gt2 d5: low share by chain 0.063 0.000 0.000 0.055; longest low stay 1060 0 0 572; long stays 2; median sigma_b 0.109 0.110 0.110 0.107
gt2 d5: measured R-hat max 1.149 (4 over 1.1) at 5,000 kept; 1.027 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.4193, measured 0.1974; coverage 0.944, width 0.743; median sigma_emb by chain 0.23 0.24 0.23 0.22; 121 minutes
gt4 d1: low share by chain 0.013 0.142 0.140 0.077; longest low stay 87 867 1461 353; long stays 3; median sigma_b 0.084 0.064 0.068 0.067
gt4 d1: measured R-hat max 1.480 (61 over 1.1) at 5,000 kept; 1.084 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.4386, measured 0.2125; coverage 0.976, width 0.854; median sigma_emb by chain 0.10 0.09 0.10 0.08; 109 minutes
gt4 d2: low share by chain 0.353 0.262 0.435 0.337; longest low stay 909 956 3000 1689; long stays 12; median sigma_b 0.028 0.033 0.025 0.031
gt4 d2: measured R-hat max 1.650 (139 over 1.1) at 5,000 kept; 1.156 (5 over 1.1) at all 20000; verdict review; RMSE embedding 0.4968, measured 0.2404; coverage 0.924, width 0.840; median sigma_emb by chain 0.34 0.41 0.22 0.33; 108 minutes
gt4 d3: low share by chain 0.292 0.180 0.134 0.231; longest low stay 1858 1107 1110 2771; long stays 9; median sigma_b 0.042 0.049 0.055 0.049
gt4 d3: measured R-hat max 1.494 (135 over 1.1) at 5,000 kept; 1.098 (0 over 1.1) at all 20000; verdict review; RMSE embedding 0.3630, measured 0.2520; coverage 0.938, width 0.962; median sigma_emb by chain 0.54 0.53 0.46 0.43; 108 minutes
gt4 d4: low share by chain 0.370 0.295 0.309 0.308; longest low stay 1525 1290 1891 1825; long stays 12; median sigma_b 0.029 0.035 0.034 0.032
gt4 d4: measured R-hat max 1.818 (170 over 1.1) at 5,000 kept; 1.185 (10 over 1.1) at all 20000; verdict review; RMSE embedding 0.4419, measured 0.2023; coverage 0.962, width 0.861; median sigma_emb by chain 0.10 0.11 0.09 0.12; 109 minutes
gt4 d5: low share by chain 0.145 0.213 0.033 0.221; longest low stay 761 612 268 2403; long stays 4; median sigma_b 0.063 0.044 0.080 0.064
gt4 d5: measured R-hat max 1.292 (52 over 1.1) at 5,000 kept; 1.070 (0 over 1.1) at all 20000; verdict ok; RMSE embedding 0.4256, measured 0.2036; coverage 0.958, width 0.828; median sigma_emb by chain 0.13 0.14 0.14 0.16; 109 minutes
```
