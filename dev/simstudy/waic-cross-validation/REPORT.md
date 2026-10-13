# Site cross-validation of the Lesson 4 factor-count comparison: step 1 report

13 October 2026. The design and decision rule were fixed in [DESIGN.md](DESIGN.md) before any fit; nothing in it was changed afterwards.

## Answer

On Lesson 4's survey, one and two site factors predict new sites equally well, and WAIC gives the same answer as direct cross-validation. Under the pre-specified rule both verdicts are "no meaningful difference", the cross-validation verdict survives its Monte Carlo error and every single-chain recomputation, and so step 2 (full leave-one-site-out) is not needed.

- Ten-fold cross-validation: elpd -11,567.54 for one factor and -11,566.95 for two. The difference (one minus two) is -0.59, with paired across-site SE 0.27 and Monte Carlo SE 0.20. On the deviance scale that WAIC uses, two factors are ahead by 1.18.
- WAIC from the full-data refits, all 24,000 draws: elpd_waic -11,560.57 and -11,559.91, difference -0.66, paired SE 0.30, Monte Carlo SE 0.21. On the deviance scale, two factors are ahead by 1.31, matching the 1.31 Lesson 4 quotes from 4,000 draws of the archived fits.
- Both methods lean slightly towards two factors, the generating count, by about two paired SEs. The difference is far below the threshold of 4 elpd (8 on the WAIC scale) under which, following the cross-validation FAQ, a difference is small whatever its SE.
- Single chains: recomputing cross-validation from each chain's draws alone gives differences of -0.83, +0.06, -0.54 and -1.04, all "no meaningful difference". Moving the full difference by two Monte Carlo SEs either way does not change the verdict either.

## The three questions

1. Can users be told to rank models by WAIC? On this survey, WAIC's total and its difference agree closely with direct cross-validation, although 99 of 100 sites in each model trigger WAIC's reliability warning (pointwise penalties 0.39 to 1.63). So the warning did not signal a misleading comparison here. This is one survey with one pair of models: it shows that the corrected `extractWAIC()` gave the right answer here, not that it is reliable in general.
2. Lesson 4's answer: no meaningful difference in predicting new sites between one and two factors. The data do not choose a factor count, and the slight lean towards two is not evidence of the generating count.
3. The held-out scoring code in `cv-lib.R` works and is tested; it is a basis for a user-facing site-level cross-validation function if one is wanted.

## Site-level agreement

- Per-site expected log predictive densities from the two methods correlate at 0.9999 in each model. WAIC is higher than cross-validation by 0.07 per site on average (about 7 in total per model), the optimism expected when cross-validation trains on 90 sites rather than 100.
- Per-site paired differences agree less well: their correlation between methods is 0.48. Both methods favour two factors at about the same number of sites (56 by cross-validation, 55 by WAIC), and the largest single-site difference is 0.10 by cross-validation and 0.13 by WAIC. WAIC tracks the totals and the site-to-site variation in fit, but not which individual sites favour which model.

## Checks

- Convergence: all 22 refits pass the strict screen (Rhat at most 1.05, bulk and tail ESS at least 400) on every occupancy, collection and detection parameter and on the fitted sites' occupancy probabilities: largest Rhat 1.012, smallest bulk ESS 568, smallest tail ESS 620. No fit raised a warning.
- Quadrature: recomputing every score on 1,000 evenly spaced draws with stricter quadrature changed site log likelihoods by at most 2.8e-6 (tolerance 1e-4).
- Monte Carlo error of held-out scores: at most 0.052 per fold of ten sites.
- Reproducibility: a full-length rehearsal (an accidentally full-schedule dry run in a scratch directory) and the run of record gave identical posterior results in all 22 fits and identical result files.
- Archived fits: the full-data refits under R 4.6.1 and occJSDM `c102745` match the archived Lesson 4 fits (R 4.5.0, `a56f554`, same seeds) to within 5.3e-15 in `B0` and 1e-15 in the factor SD, and exactly in `theta0` and `p`. The package changes since `a56f554` did not change these fits, and R 4.6.1 changes only their last bits.
- Scorer tests: `test-cv-lib.R` passes 43 expectations, including exact reproduction of `computeSiteWAIC()` on training sites and agreement to 1e-6 with a brute-force likelihood on withheld sites.

## Correction to DESIGN.md

DESIGN.md says single-threaded fits are "the only setting at which fits are bit-reproducible", citing AGENTS.md. That finding dates from August. `?runOccJSDM` now says random draws run serially on the main thread whatever the RcppParallel thread count, so the thread count no longer changes the sampling stream. Running single-threaded remained the right choice for eight fits side by side; the stated reason was out of date. DESIGN.md is left as written, since it was fixed before the run.

## Limits

One survey (100 sites, 10 species, three samples, two primers, six PCRs), one pair of models and one simulated community. Ten-fold cross-validation trains on 90 sites, which lowers both models' scores slightly; the paired difference is much less affected. Spatial models are outside the scope: the site score rejects them.

## Run

- Software: R 4.6.1 and occJSDM `c102745`, installed from a `git archive` export into the study library; scripts at `819f9a4`.
- 22 fits of four chains, 3,000 burn-in and 6,000 retained iterations, one thread each, eight at a time: 208 to 259 seconds per fit. Fits finished in 12 minutes and scores in 18 more; the two-factor full-data WAIC took 18 minutes alone.
- Raw fits, scores, setup and logs: `dev/simstudy/results/waic-cv-20261013/` (ignored). Compact results here, in `results/`: `decision.txt`, `totals.csv`, `per-chain.csv`, `agreement.csv`, `sites.csv`, `fits.csv` and `archive-comparison.csv`.

Reproduce from the repository root, with `out` a new directory holding `library/` (this package installed from a `git archive` export) and `source-commit.txt`:

``` sh
R_LIBS="$out/library" Rscript dev/simstudy/waic-cross-validation/test-cv-lib.R
R_LIBS="$out/library" Rscript dev/simstudy/waic-cross-validation/prepare.R "$out"
mkdir -p "$out/scripts" && cp dev/simstudy/waic-cross-validation/*.R dev/simstudy/waic-cross-validation/launch.sh "$out/scripts/"
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' \
  sh "$out/scripts/launch.sh" "$out" dev/simstudy/waic-cross-validation/results \
  > "$out/logs/launch.log" 2>&1 < /dev/null &
```
