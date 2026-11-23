# Simulation study of the grouped covariate prior: output contract, version 2

Status: written by occJSDM (Claude, for Doug) on 21 November 2026. Not issued. It comes into force with the design contract it pairs with, `simulation-grouped-prior-design-v2.md`, when Doug says so.

This is what the dataset section hands back to occJSDM. The real-data output contract's rules about location, the run record and what must not appear apply here too. Against version 1 it adds the tables and fields that the new convergence rule and the two new criteria need. One thing that version 1 asked for changes: the measured-only AUC is no longer found by a pre-study, and is computed in the run itself. One table that was optional is now required (`replicates`).

## Where it lives

In the dataset repository, in a study folder under its own outputs and analysis folders. occJSDM is told the path and reads it in place.

## What is handed back

1. A run record, `run_record.json`, with the keys of the real-data output contract and the additions below.
2. Summary tables as CSV files, one row per truth, arm or comparison, with no site identifier.
3. A written result, `RESULT.md`.

Per-fit extracts and logs stay in the dataset repository.

## Run record additions

- `design`: the number of sites and species S, the measured covariates, the embedding set and its columns, and the data checksums.
- `truths`: for each truth, the strengths used (embedding and measured variance, leftover share, blends, active axes, management correlation), the realised variances, a statement that they are the strengths of the version 1 run, and the measured-only held-out AUC of the truth in this run.
- `split`: the block side, the seeds, the number of species dropped for rarity per replicate.
- `seeds`: a statement that no seed was used in an earlier run, and how that was checked.
- `replicates`: the number per cell in the first pass and in the confirmation, and which cells were confirmed.
- `pilot`: the time per fit by arm with 4 chains, and the projection that was made from it.
- `chains`: the number of chains by arm, and whether they ran one after another or on separate cores.
- `library`: the package commit the grouped fits were made with, and the evidence that the move on `sigma_b` ran in every grouped binary fit.
- `rules`: both rule versions by name, the date each was written and the date the fits were made, and the limits used (1.1, 0.30, 0.90, the package's small-effect line, the 10% limit for a cell, and the tolerance of 0.05 for the gradients criterion).
- `package_check`: the evidence that the study's gradients calculation reproduces `returnEmbeddingGradients()` on the grouped arm.

## Summary tables

Every criterion table is given on three fit sets, each with its number of pairs: all fits, the fits converged under the version 2 rule, and the fits converged under the version 1 rule. A pair belongs to a converged set if its grouped fit converged. Beside each comparison, on each fit set, every criterion table gives the number of pairs in which the other arm's fit failed its own convergence check.

- `prediction`: per truth and arm, held-out log score, and the paired differences grouped against components and against measured only, overall and by prevalence band, with intervals and the verdict (better, not worse, worse).
- `rarity`: signed and absolute error of predicted prevalence by arm, truth and band, and the grouped-to-components ratio in the lowest band.
- `management_coverage`: for truth 4, coverage of the 95% interval at the three levels (species, community average, trait slope), by arm and by true effect (zero or non-zero).
- `management_false_positive`: under a true effect of zero, the share of species-replicates whose interval excludes zero.
- `gradients`: for truths 1 and 4, per cell and arm, the adjusted blend R-squared at the held-out sites with each arm at its own number of gradients, and the paired differences grouped against each components arm and against raw, with intervals and the verdict.
- `gradients_detail`: the items the design contract lists as reported and deciding nothing. Every arm at each number of gradients from 1 to 10, adjusted and unadjusted, at the held-out and the fitted sites. Each arm's own number of gradients, as counts of fits by number. The reference with no fit. The grouped arm by state. The share of each arm's effect held by its leading 3 gradients.
- `shared_part`: per cell, the grouped fits in the shared, mixed and leftover states, what the truth holds, the number in the right state, and the verdict (pass, fail, not scored).
- `convergence`: per cell, for the grouped arm: the fits failing version 1, by what failed (an axis scale, `sigma_emb`, `sigma_b`, a management coefficient); the fits failing version 2 on the targets, by part (a) to (e), with part (a) by level; the fits not assessed on (c) and (d); the fits whose split can be read; and the verdict of each version for the cell. For every other arm: the fits failing parts (a) and (b).
- `convergence_cross`: per cell, the four counts of version 1 pass or fail against version 2 pass or fail. Every fit that fails version 1 and passes version 2 is listed by cell and replicate number, with what failed version 1 and its largest R-hat. Nothing that version 1 flagged is left uncounted.
- `convergence_detail`: per cell, the quantities the design contract lists as reported and deciding nothing, as medians and largest values over fits. The largest difference between two chains' embedding effects is given for the assessed and the not-assessed fits separately.
- `coefficients`: coverage and bias of the other measured coefficients.
- `ess`: effective sample sizes of the same quantities.
- `timing`: wall time per fit by arm.
- `replicates`: one row per replicate, truth and arm, with its log score, the verdict of each rule and what failed it, the split label, the state and its flags, with no site identifier. Required in version 2, because the verdict of the version 1 rule on every fit must be reported. It was optional in version 1.

## Per-fit extract

These fields are new against the version 1 run. They stay in the dataset repository, and the tables above are made from them.

- R-hat of the community-average and the trait-slope management coefficient, each computed on the quantity whose coverage is judged.
- Each chain's held-out log loss.
- The largest difference between two chains' whole embedding effects at the fitted sites, on the scale of the linear predictor, whether or not the package assessed the fit.
- R-hat of each species' embedding contribution at the fitted sites, summarised as the largest value and the share over 1.1.
- The gradients measure at each number of gradients from 1 to 10, and the fit's own number.

## Written result

`RESULT.md`, one line per paragraph. The result says plainly that the findings are conditional on the fixed design of the real sites and covariates, and that version 2 was issued after the results of a version 1 run had been seen.

For every cell, in this order.

- The number of grouped fits in the shared, mixed and leftover states, beside what the truth holds, and the verdict of the shared-part criterion (passes, fails, or not scored). A failure is given in the word "fails", whatever the cell's other verdicts.
- The version 1 verdict on convergence, in version 1's own words ("fails on mixing grounds", or passes).
- The version 2 verdict on convergence. A cell that passes on the targets and fails the shared-part criterion is written as "converged on the targets; fails the shared-part criterion (n of N fits in the right state)", never as plainly "converged". Version 2's verdict is never given without version 1's beside it.
- The package's own verdict ("ok" or "review"), as a count of fits.
- The number of fits not assessed, in the words "the chains' agreement on the embedding was not checked in n of N fits".
- Each criterion: pass or fail, and the interval behind it, on the three fit sets. A comparison with the raw arm says in words that it is against fits that failed their own convergence check, where they did.

Then the cost, and what a larger study should assume. A failure states the smallest synthetic check that reproduces it, with no dataset content.

## What occJSDM does on receipt

It checks the run record for the package commit, the rule versions and the seeds, reads the tables, and decides whether a failed criterion is a package defect, a default to change (for example the axis scale or the tail of the prior), or a limit of the design. A failure on mixing grounds, or of the prediction, rarity, management or gradients criterion, at 100 replicates halts part 2 of the package design, as the package spec says. A failure of the shared-part criterion does not. It withdraws the claims about the shared part for that cell, as the design contract says.
