# Simulation of the sigma_b move, 22 November 2026

This records steps 1 and 2 of the evidence in `dev/superpowers/specs/2026-11-22-sigma-b-move-design.md`, run as Task 7 of `dev/superpowers/plans/2026-11-22-sigma-b-move.md`. Step 3, the replay by the AlphaEco project, is not part of this record.

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
