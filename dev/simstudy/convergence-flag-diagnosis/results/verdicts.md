# Task 3c verdicts: community 5 diagnostic fits

The frozen decision rules of [README.md](../README.md) (frozen at 183bc2f), as amended by [AMENDMENT-1.md](../AMENDMENT-1.md) (rulings R13 to R18), applied by `analyse.R` to the 40 diagnostic chains of `design-qfar_K6-sites300-05`: the extended run (16 chains) and variants (a), (b) and (c) (8 chains each), each chain 10,000 retained draws from 40,000 post-burn iterations. The verdicts are mechanical; the last section is interpretation and is marked as such. The numbers are in `extended/verdict.csv`, `variants/verdicts.csv` and the CSVs beside them.

## Checks before any rule

- The md5s AMENDMENT-1 records match the files at HEAD: `modes.R` bab7c1c3934ba2e0c83b8f26232a0aaf, `anatomy.R` bd98ff3f40cad972f96bdeca5a07e206 and `results/modes/anchor-classifier.csv` e235c2fa641eb05c36232bec0d513bc8; so does the launcher's frozen md5 of `run.R` (f9e167d8cc8879e27b32f80a4135eb26), whose functions supply the frozen seeds, schedules and variant definitions.
- All 40 fits match the frozen run, chain, seed (20261001 plus the chain), seeded RNG state, schedule (10,000 burn-in, 10,000 retained, thinning 4), input md5, library revision 707540a, launcher-frozen script hashes, `listPriors` and variant definition; variant (b) holds species 6 `theta0` at 0.038004923556 in every draw, and variant (c) holds the nine species without `OTU_6`. No fit has a package warning. Per-fit md5s are in each `provenance.csv`.
- The source input has md5 574ed4df46b9c1bd227c79d2185a7fcb, read through the remap with md5 refusal; the saved derived input of variant (c) has md5 b2c80630eedaa0c1692a8aa9535d0571 and is identical to a fresh `drop_species()` of the source input.

## Extended run: one mode only: near-truth

- Rule (AMENDMENT-1, R15): one mode only when no chain visits the other region, that is all 16 chains stay in the same anchored region; the verdict names the region. It is checked after Slow mixing (every chain visits both regions) and Separated modes (each chain in one region, each region holding a chain), neither of which applies, and before Mixed.
- Numbers: all 16 chains are in the near-truth region, each with 10,000 of 10,000 draws near-truth and none mirror; the far-from-both share per chain is 0.57% to 1.37% (R18 limit 5%), so no chain is unknown. Over all 160,000 draws the near-truth share is 1, the mirror share 0, and 0.97% of draws are far from both components.
- Beside it (not deciding): the largest all-chain rank-normalised Rhat over the eleven species-6 quantities is 1.0015 (`B0`) and the largest per-chain split-Rhat 1.0074 (chain 8, `B0`); Task 1's `chain_separation()` labels every one of the ten species `agrees`, with species 6 separation at most 0.14.
- Secondary checks agree with the primary on every chain: the `theta0` cut at 0.135 puts every chain in the near-truth region (0.01% to 0.18% of a chain's draws above the cut, 177 of 160,000 in all), and the refitted mixture finds one mode. No chain-level disagreement between the primary and either secondary.

## Variant (a), Beta(1, 100) theta0 prior: explains, flagged uninformative

- Rule (README.md with R14): the variant explains the split if every species-6 quantity agrees, that is its Task 1 label is `agrees` and its all-chain Rhat is at most 1.05. AMENDMENT-1: since the extended run shows one mode only, the verdict is flagged as uninformative about the split.
- Numbers: 11 of 11 quantities agree; largest all-chain Rhat 1.0019 (`mean_psi_original_sites`), largest separation 0.12.
- Beside it (not deciding): all 8 chains in the near-truth region (far-from-both 0.07% to 0.31%); the `theta0` cut and the refitted mixture (one mode) agree on every chain.

## Variant (b), species 6 theta0 fixed at its generating value: explains, flagged uninformative

- Rule: as for (a), with `theta0` left out as fixed by design.
- Numbers: 10 of 10 remaining quantities agree; largest all-chain Rhat 1.0015 (`mean_psi_original_sites`), largest separation 0.10.
- Beside it (not deciding): all 8 chains in the near-truth region (far-from-both 0.39% to 1.09%); the refitted mixture (with `theta0` dropped as constant) finds one mode and agrees on every chain. The `theta0` cut does not apply.

## Variant (c), species 6 removed: isolates, flagged uninformative

- Rule (R13): for each of the other nine species, the change in the pooled posterior mean of `mean_psi_original_sites` (8 variant chains minus 16 extended chains) is below 0.01 in absolute value, or within 2 Monte Carlo standard errors of the difference (the two runs' MCSEs combined in quadrature). Flagged uninformative as for (a).
- Numbers: 9 of 9 species pass, all by the 0.01 condition; the largest absolute change is 0.0082 (`OTU_1`, combined MCSE 0.0017), and `OTU_7`, the species that drifts in the pr11 fit, changes by 0.0059 (combined MCSE 0.0072). Six species (`OTU_2`, `OTU_4`, `OTU_5`, `OTU_7`, `OTU_9` and `OTU_10`) also lie within 2 MCSEs; `OTU_1`, `OTU_3` and `OTU_8` pass by 0.01 alone.
- Beside it (not deciding): the largest change over the 300 sites in a species' posterior mean occupancy probability is 0.067 (`OTU_9`, site 19).

## Posterior mass and species 6 against the generating values

- Mass shared between the modes (extended run, primary assignment): near-truth 1, mirror 0; far from both 0.97%.
- Species 6 in the near-truth region (all 16 chains) against the generating values: `theta0` mean 0.036 (truth 0.038), `B0` -0.48 (truth -0.37), collection intercept -0.66 (truth -0.82), `B_slope1` 0.60 (truth 0.70), `B_slope2` -1.78 (truth -1.87), and mean occupancy over the original sites 0.475 (truth 0.511); every generating value lies inside the pooled 95% interval. There is no mirror region to compare (`extended/region-means-vs-truth.csv`).

## Interpretation (not part of the verdicts)

Superseded: Q4 of [REPORT.md](../REPORT.md) replaces the first bullet below. It limits the variants' absence of the mirror to variants (a) and (b), lists thinning and the chains' common starting values among the differences, drops "chance alone is possible but unlikely", and counts the pr11 initial fit's 2 chains beside these 20.

- The mirror mode that held chains 2 and 4 of the saved pr11 fit did not appear in any of the 16 extended chains, nor in any chain of the variants, so this run cannot say what produces it. The diagnostic chains differ from the pr11 fit in more than their number: each runs in its own process with its own seed, and with a longer burn-in (10,000 against 6,000 iterations); thinning by 4 changes only which draws are kept. Treating all 20 chains as exchangeable, the chance that both mirror chains fall among the 4 pr11 chains is 6 in 190 (0.032, Fisher's exact test, one-sided), so chance alone is possible but unlikely; which of these differences matters, if any, is not tested here.
- Variant (a) moves species 6 `theta0` towards zero (mean 0.011 against 0.036 in the extended run, truth 0.038) and `B0` up (-0.11 against -0.48), with mean occupancy 0.52 (truth 0.51); this is the prior's effect within the near-truth region, not a change of mode.
