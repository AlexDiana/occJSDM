# A convergence check for the grouped prior: short design

19 November 2026. Branch `claude/grouped-mixing`. Status: awaiting Doug's approval. Nothing is built.

## Why

R-hat on single embedding coefficients is the wrong test for a fit made with the grouped prior. It fails good fits, because single blends and single embedding coefficients are not pinned down and their draws are heavy-tailed. It passes the dangerous fit, in which every chain ends with the leftovers carrying the species' responses and no gradients are found, because those chains agree with each other. The evidence is in `dev/simstudy/grouped-mixing-20261118/README.md`. The purpose of the prior, as Doug stated it, is to learn habitat gradients when they exist, so the check has to say whether the chains found gradients and whether they found the same ones.

## What is asked for

A way to tell a good grouped-prior fit from a bad one, run automatically at the end of a fit and available as a function.

## Approaches

- A new function only. Nothing existing changes, but every grouped fit would still end with the misleading warnings and stay silent on the dangerous case.
- A new function, and the automatic report at the end of a fit changed for grouped fits with an embedding. Doug chose this on 19 November 2026.
- Replacing R-hat throughout with a rank-normalised version. Set aside: it changes every fit's report, and by itself it does not address either failure.

## What gets built

One exported function, `returnGroupedConvergence(fitModel)`, for any fit made with `embeddingVars`, whatever its model type. It returns a list: `verdict` ("ok" or "review"), `checks` (a data frame with one row per check: what was measured, its value, the limit, whether it is flagged, and a sentence saying what a flag means), `chains` (one row per chain), and the pieces behind checks 3 and 4. `runOccJSDM()` prints the same summary at the end of such a fit and gives one warning when the verdict is "review".

Four checks.

1. Which state each chain is in. The median of `sigma_emb` in each chain. A chain with a median of 1 or more is one where the leftovers are at least as large as the shared part. Flagged when the chains are in different states, and when all of them are in that state. The second is not always wrong (species may truly respond independently), so the sentence says so, and says that with fewer than 5,000 burn-in iterations a missed gradient is the more likely reading.
2. R-hat, the package's own, on the measured coefficients and on the logs of the scales (`sigma_b`, `sigma_emb`, the axis scales). Flagged over 1.1. This is the only check of the four that is an R-hat, and it needs each chain's variance as well as its mean; the others use only each chain's posterior mean.
3. Agreement on the set of gradients. For each chain, the posterior mean of the embedding's combined effect (sites by species). The number of gradients is the smallest number of directions that hold 90 percent of that effect pooled over chains. For each pair of chains, the share of one chain's effect that lies within the other's gradients. The smallest value over pairs is reported, and flagged under a limit to be fixed from the saved fits before the code is written (provisionally 0.90).
4. The size of any disagreement in the combined effect. How far each chain's effect is from the pooled one, as a share of the effect's size. Flagged over 0.15.

Checks 3 and 4 are not assessed, and say so, when the embedding's effect is small (provisionally, when its spread across sites is under 0.5 on the scale of the linear predictor for the typical species), as in a fit where no species responds to the embedding. They need two chains or more.

The automatic report. `computeDiagnostics()` gains an optional argument naming the embedding rows of B. With it, the measured rows of B are reported and warned on as now; the embedding rows are reported without a warning and with a line saying they are not a test of convergence; and the scales are reported on the log scale. Without the argument, and for every fit made without `embeddingVars`, it does exactly what it does today.

## Amendment, 19 November 2026: beyond the fitted sites

Made while the function was being built, after Doug pointed out that different blends must extrapolate differently. Measured on saved fits (per-chain posterior means of the coefficients against the truth): in a fit whose two chains found gradients and differ by 6 percent of the effect's size at the fitted sites, they differ by 6 percent at new sites like the fitted ones, by 11 percent where the embedding's weak directions are three times larger, and by 16 percent where every direction is equally large, and the error against the truth goes from 27 to 42 percent. So checks 3 and 4 are true of the fitted sites and of sites like them, and understate the disagreement elsewhere. Two changes follow. The wording of both checks says "at the fitted sites". And a further row is reported, `disagreement_beyond`: how far each chain's embedding coefficients are from the pooled ones, as a share of their size, which is the disagreement at sites that vary equally in every embedding direction. It is reported and never flagged, because there is no evidence yet for a limit. Checking the chains at the sites a map will be drawn for is left for later: it needs the new sites put through the fit's own standardisation.

## Amendment, 19 November 2026: after the independent review

The review of the finished function found faults in two definitions above, and the full-length fits run afterwards showed a third. The record and the numbers are in `dev/simstudy/grouped-mixing-20261118/README.md`. Each of the three changes below alters something Doug approved. He confirmed all three on 19 November 2026.

- Check 3. The agreement of one chain with another is now the share of its effect within the other's gradients divided by the share within its own. As first defined, two identical chains scored only the share their gradients hold, about 0.90, which was also the limit. The limit stays at 0.90. No saved fit has chains with clearly different gradients, so no saved fit tests it.
- Check 4. The disagreement is now the largest difference between two chains' effects, as a share of the effect's size. As first defined (each chain's distance from the pooled effect, averaged) one chain apart was diluted by the chains that agree. The limit is 0.30, which is the old 0.15 carried over to the new statistic. `disagreement_beyond` is defined the same way.
- Check 2. The R-hat of the scales is in two rows. `rhat_scales` covers `sigma_emb` and `sigma_b` and is flagged over 1.1. `rhat_axis` covers the axis scales, is reported with the number over 1.1, and does not decide the verdict: in full-length fits whose chains agree on everything reported, 2 to 5 of 64 axis scales fail R-hat.

Also: the step at the end of a fit can no longer lose the fit (an error becomes a warning); draws with missing values flag the check they belong to; and nothing copies B.

## What checks 3 and 4 do not say

They compare chain means. They say whether the chains reached the same answer and how far apart they are. They do not say whether each chain has explored its own uncertainty, which is what R-hat is for, and which is why check 2 stays for the quantities reported with intervals.

## Not in this change

- An R-hat on the combined effect at a sample of sites and species. It failed 79 percent of the time in good fits on small differences. It can be added later as a reported number that does not decide the verdict.
- Any change to the report of a fit made without `embeddingVars`.
- A rank-normalised R-hat.

## Calls made here, not by Doug

- The name `returnGroupedConvergence()`, after `returnConvergenceDiagnostics()`.
- The limits: 1 for the state, 0.90 (provisional) for the gradients, 0.15 for the disagreement, 0.5 for a small effect, 90 percent for the number of gradients. The second and third will be set against the 18 saved fits whose state is already known, which is the only evidence there is and is circular to that extent; the document that records the result will say so.
- "review", not "fail", as the word for a flagged fit, because the first check can flag a fit that is right.
- Reporting the scales on the log scale in the automatic report of grouped fits.

## How it will be tested

Tests first. With made-up effects: identical chains give agreement 1 and disagreement 0; chains whose blends differ but whose gradients are the same give agreement 1; a chain with other gradients is flagged; a small effect is not assessed. The state check on made-up draws of `sigma_emb`. `computeDiagnostics()` warns on a measured coefficient and not on an embedding coefficient, and is unchanged without the new argument. A grouped binary fit returns the list and prints the summary; a fit without `embeddingVars` is refused. Then the function's two effect-based checks are run on the saved fits of the rotation and orientation tests, where which fits have a chain apart is already known.
