# Amendment 1 to the Task 3 protocol

Dated 30 September 2026, and made before any diagnostic fit was launched. The seed base 20261001 and the archive name `convergence-diagnosis-20261001` are fixed labels chosen in advance, not dates. It amends [README.md](README.md) as frozen in commit 183bc2f587f1c686d9caa1f5236f87e812fbe4fa (README.md md5 at the freeze `8b63ba629693b7adb770fa7171e96062`). Where the two differ, this file holds. README.md is otherwise unchanged apart from a one-line pointer to this file at its top, so `git diff 183bc2f -- README.md` shows only that line and a blank line after it. This is the dated amendment file README.md provides for, numbered rather than named by its date. Rulings R13 to R17 are the controller's, made after review of the frozen protocol.

Unchanged: the fits, variants, schedules, seeds, library, `run.R`, `fixed-theta0.R`, `launch.R`, and every file hashed into a fit or frozen by the launcher (the launcher's frozen md5s still match). The launch command in README.md stands.

## Why

Review found that the frozen mode assignment, the mixture refitted to each run's pooled draws (README.md, "How the rules are measured", the "Regions" and verdict bullets, lines 75 and 76 at the freeze), misassigns draws when one mode holds only one or two of 16 chains. On pseudo-chains cut from the saved pr11 fit, 15 near-truth and 1 mirror chain gave a wrong second component (all 15 majority chains then visited both regions, up to 7.6% of a chain's draws off-region), 1 and 15 gave a single mode, and 14 and 2, or 7 and 1 of 8, failed as well. The causes are EM stuck in a poor solution, and a mixture likelihood that prefers splitting the majority mode's non-Gaussian shape to separating the minority chains. Review also found that the single-mode case named no mode, and three readings of the frozen rules were settled at the same time.

## R17: the primary mode assignment is anchored, not refitted

This supersedes the "Regions" bullet of README.md (line 75 at the freeze) and the use of the refitted mixture in the verdict bullet (line 76).

- The primary assignment uses a frozen two-component classifier calibrated on the saved pr11 fit of `design-qfar_K6-sites300-05` (fit md5 `07b2d9eb4356b7973f87485a0f8bcd09`), species 6: component 1, **near-truth**, is the mean and covariance of chains 1 and 3; component 2, **mirror**, is the same for chains 2 and 4; the weights are equal.
- It uses only quantities that stay free in every variant: species 6's two environmental slopes (`B_slope1`, `B_slope2`), its collection intercept (`beta_theta_intercept`) and `mean_psi_original_sites`. It does not use `theta0` (fixed in variant b and squeezed by the prior in variant a) or `B0`; both are reported for each mode, as chain means by region.
- The classifier is `results/modes/anchor-classifier.csv` (md5 `e235c2fa641eb05c36232bec0d513bc8`): each component's weight, 24,000 calibration draws, four means and 4 x 4 covariance, each number as a decimal and as an exact hexadecimal float, written by `verify.R --mode=anchor-calibration`. `read_anchor()` in `modes.R` refuses the file if its md5 differs. Component means (near-truth, then mirror): `B_slope1` 0.603 and -1.266, `B_slope2` -1.793 and 2.018, `beta_theta_intercept` -0.664 and -2.832, `mean_psi_original_sites` 0.475 and 0.448; the squared Mahalanobis distance between the means is 74.9.
- `assign_anchored()` gives each retained draw the component of higher Gaussian density (equal weights); `anchored_regions()` gives each chain's share of near-truth and mirror draws and its region. As before, a chain visits a region if at least 1% of its retained draws (100 of 10,000) carry it, and stays in one region if the other holds less than 1%.
- Validation on the saved fit (`results/modes/anchor-validation.csv`): 1 of 48,000 draws is labelled against its chain's region in sample; calibrated on chains 1 and 2 and applied to chains 3 and 4, or the reverse, or on the first halves of every chain and applied to the second halves, at most 1 draw per chain is labelled against its region. On pseudo-chains of the same fit (`results/modes/imbalance-check.csv`) every chain is assigned correctly in all seven splits (15 and 1, 14 and 2, 8 and 8, 2 and 14, 1 and 15 of 16; 7 and 1, 1 and 7 of 8), with at most 0.07% of a chain's draws off-region. `test-modes.R` holds these splits, a set with `theta0` held constant as in variant (b), and a classifier calibrated on chains 1 and 2 applied to chains 3 and 4, as regression tests.
- The share of each chain's draws far from both components is computed with the primary assignment and, under R18 below, can place a chain in an unknown region.
- Secondary checks, reported alongside the primary assignment: (i) a `theta0` cut at 0.135, the valley of the pooled pr11 `theta0` density (0.1351), mirror above it, applied only where `theta0` is free (the extended run and variant a); on the saved fit 0.03% to 0.09% of a chain's draws fall on the wrong side, and it assigns all seven pseudo-chain splits correctly. (ii) The refitted mixture of README.md, now also started from chain partitions (each chain against the rest, and the chains split at the largest gap in their chain means of each quantity); this fixes the 14 and 2, 2 and 14, 7 and 1 and 1 and 7 splits but not 15 and 1 or 1 and 15. `compare_assignments()` gives, per chain, each method's region and share of mirror draws and whether it agrees with the primary; every chain where a secondary disagrees with the primary is reported. The verdict uses the primary only.

## R18: a region unlike both anchors

Added on 30 September 2026 after re-review, still before any diagnostic fit. The anchored classifier always assigns the component of higher density, so a chain in a region unlike both anchors would otherwise still be called near-truth or mirror: in a probe, a near-truth chain with occupancy shifted up by 0.3 had 68% of its draws far from both components and was classed 100% near-truth.

- A draw is far from both components when its squared Mahalanobis distance to each component (using that component's mean and covariance) exceeds the 99.9% quantile of the chi-squared distribution with 4 degrees of freedom, one per anchor quantity (18.47). This definition is unchanged from `assign_anchored()` as first written.
- A chain with more than 5% of its retained draws far from both components (`UNKNOWN_SHARE` in `modes.R`) is in an **unknown** region, whatever its near-truth and mirror shares. `anchored_regions()` then gives its region as `unknown`, `region_pattern()` returns `other`, so the extended-run verdict is Mixed, and `region_note()` names the chain or chains in the report.
- Reference, from the saved pr11 fit the anchor was calibrated on: 0.36% to 0.94% of each whole 12,000-draw chain is far from both components (0.27% to 1.03% if thinned by 4 like the diagnostic chains). The share varies with the length of the stretch examined, because far draws come in runs of up to 23 iterations: over consecutive blocks of 1,500 draws the largest share is 5.3% (1 block of 32 above 5%), over blocks of 3,000 draws 2.7%, 6,000 draws 1.5% and 12,000 draws 0.9% (`results/modes/anchor-validation.csv`). A diagnostic chain spans 40,000 post-burn iterations, so a chain from either known region is expected well below 5%.
- The far-from-both share of every chain is written to `mode-mass.csv` (column `share_far`) and carried by `compare_assignments()` (column `primary_share_far`). The primary, anchored assignment is the one written to `mode-mass.csv` and the one that colours the 16-chain figure (`plot_chain_strips()`, where a chain in an unknown region is drawn in grey with its far-from-both share); the secondary assignments are reported beside it.
- In the pseudo-chain checks the 1,500-draw block above 5% is classed unknown (splits 2 and 14, and 1 and 15, of 16; `results/modes/imbalance-check.csv`); every chain's majority component is still correct, and the regression tests check the classification of such short blocks separately from the unknown rule.

## The extended-run verdict (replacing README.md lines 76 and 77 at the freeze)

- **Slow mixing**: every one of the 16 chains visits both anchored regions, and the rank-normalised all-chain Rhat (`rhat_all_chains`) is at most 1.01 for all eleven species-6 quantities.
- **Separated modes**: every chain stays in one anchored region, each region holds at least one chain, and every chain's own split-Rhat (`split_rhat_within_chain`) is at most 1.05 for all eleven species-6 quantities.
- **One mode only** (R15): no chain of the extended run visits the other region, that is all 16 chains stay in the same anchored region. The verdict names that region: "one mode only: near-truth" or "one mode only: mirror". This replaces README.md line 77, which gave Mixed for this case.
- **Mixed**: anything else, including any chain in an unknown region (R18).
- `region_pattern()` in `modes.R` reads the region part of these rules from `anchored_regions()`: every chain visits both, each chain in one region, one mode only (named), or other.

## R14: reading `agrees`

A species-6 quantity agrees if its Task 1 label (from `chain_separation()` on the variant's chains) is `agrees` and its all-chain rank-normalised Rhat is at most 1.05. This confirms the reading in README.md line 78 at the freeze. A quantity held fixed by design (`theta0` in variant b) is still left out. Variant (a) or (b) explains the split if every other species-6 quantity agrees. The anchored regions of the variant's chains, and for variant (a) the `theta0` cut and for both the refitted mixture, are reported beside the verdict but do not decide it.

## R13: variant (c)

Variant (c) isolates the split if, for each of the other nine species, the change in its occupancy mean (the pooled posterior mean of `mean_psi_original_sites` over the 8 variant (c) chains minus that over the 16 extended-run chains) is less than 0.01 in absolute value, or lies within 2 Monte Carlo standard errors of the difference, the two runs' `posterior::mcse_mean` combined in quadrature. Both conditions, the change and its standard error are reported for each species. This replaces the threshold-only rule of README.md line 79 at the freeze. The reason: species 7 drifts in the pr11 fit (chain means 0.477 to 0.585), so sampling noise alone could exceed 0.01.

## R16: the pre-freeze change of mode method

The change of `assign_modes()` made before the freeze (from the logit scale with one EM start to the quantities' own scales with several starts, `b889208`), tested only on the already-known saved pr11 fit and disclosed in README.md, is accepted. Under R17 that method is a secondary check.

## Variant verdicts when the extended run shows one mode

If the extended-run verdict is one mode only, the variant verdicts are still computed and reported, but flagged as uninformative about the split, since there is then no split in the extended run for a variant to explain or isolate.

## Files Task 3c must check before applying the rules

Task 3c stops before applying any rule unless these md5s match (updated for R18): `modes.R` `bab7c1c3934ba2e0c83b8f26232a0aaf`, `anatomy.R` `bd98ff3f40cad972f96bdeca5a07e206`, and `results/modes/anchor-classifier.csv` `e235c2fa641eb05c36232bec0d513bc8`. The launcher already guards `run.R`, `fixed-theta0.R`, the library fingerprint and the clone archive.

## Recovery: a stale derived-input lock

The launcher refuses to start only on locks under `ARCHIVE/fits`. Each variant (c) process also takes `ARCHIVE/inputs/qfar_K6-sites300-05-without-species6.rds.lock` for the few seconds it spends checking the saved derived input (which already exists, written by the pilot). If a variant (c) process is killed while holding it, every later variant (c) chain waits 600 seconds for it and then fails. Before launching, and before any relaunch, check with `ls -d ARCHIVE/inputs/*.lock`; if one exists and `pgrep -fl convergence-flag-diagnosis/run.R` prints nothing, remove it with `rmdir` and then launch.

## Disk

Nothing else should write to the volume while the diagnostic fits run. When R18 was added about 76 GiB were free (12 GiB when this amendment was first written), and the 40 fits need about 2.3 GB.
