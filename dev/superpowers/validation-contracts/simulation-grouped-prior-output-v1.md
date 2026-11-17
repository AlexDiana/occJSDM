# Simulation study of the grouped covariate prior: output contract, version 1

Pairs with `simulation-grouped-prior-design-v1.md`. This is what the dataset section hands back to occJSDM. The real-data output contract's rules about location, the run record and what must not appear apply here too.

## Where it lives

In the dataset repository, in a study folder under its own outputs and analysis folders. occJSDM is told the path and reads it in place.

## What is handed back

1. A run record, `run_record.json`, with the keys of the real-data output contract and these additions.
2. Summary tables as CSV files, one row per truth, arm or comparison, with no site identifier.
3. A written result, `RESULT.md`.

Per-fit extracts and logs stay in the dataset repository.

## Run record additions

- `design`: the number of sites and species S, the measured covariates, the embedding set and its columns, and the data checksums.
- `truths`: for each truth, the strengths used (embedding and measured variance, leftover share, blends, active axes, management correlation), the realised variances, and the measured-only held-out AUC found by the pre-study.
- `split`: the block side, the seeds, the number of species dropped for rarity per replicate.
- `replicates`: the number per cell in the first pass and in the confirmation, and which cells were confirmed.
- `pilot`: the pre-study's time per fit by arm and the projection that was made from it.

## Summary tables

- `prediction`: per truth and arm, held-out log score, and the paired differences grouped against components and against measured only, overall and by prevalence band, with intervals and the verdict (better, not worse, worse).
- `rarity`: signed and absolute error of predicted prevalence by arm, truth and band, and the grouped-to-components ratio in the lowest band.
- `management_coverage`: for truth 4, coverage of the 95% interval at the three levels (species, community average, trait slope), by arm and by true effect (zero or non-zero).
- `management_false_positive`: under a true effect of zero, the share of species-replicates whose interval excludes zero.
- `convergence`: per cell, the share of fits failing, on the scales and on the management coefficients, and criteria on all fits against converged fits.
- `coefficients`: coverage and bias of the other measured coefficients.
- `ess`: effective sample sizes of the same quantities.
- `timing`: wall time per fit by arm.
- `replicates`: one row per replicate, truth and arm with its log score and flags, optional, with no site identifier.

## Written result

`RESULT.md`, one line per paragraph. For each criterion: pass or fail, in which truths and cells, and the interval behind it. Then convergence, then the cost and what a larger study should assume. A failure states the smallest synthetic check that reproduces it, with no dataset content. The result says plainly that the findings are conditional on the fixed design of the real sites and covariates.

## What occJSDM does on receipt

It checks the run record for the package commit, reads the tables, and decides whether a failed criterion is a package defect, a default to change (for example the axis scale or the tail of the prior), or a limit of the design. A failure at 100 replicates halts part 2 of the package design, as the package spec says.
