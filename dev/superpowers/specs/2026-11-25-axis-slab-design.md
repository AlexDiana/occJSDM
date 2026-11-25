# A slab for the embedding scales of the grouped prior

Date: 25 November 2026. Status: design written by Claude and revised the same day after a reviewer read it as a builder would (what changed is listed under *Decisions, and who made them*). Doug approved it on 25 November 2026, with the cap's default at 2, and asked for the plan, which is `dev/superpowers/plans/2026-11-25-axis-slab.md`. His approval covers the design and the writing of the plan, not the build. Nothing is built. Branch `claude/axis-slab`, cut from the head of draft PR 46 (`claude/sigma-b-move`, a03a1b8), not pushed.

## Purpose

Remove the wide upper tail of the grouped prior on the embedding rows of `B`, without moving the rest of the prior.

Under the present prior the embedding's contribution to a species' linear predictor has a prior median of about 1.5 logits and a 95% point of about 27. One prior draw in twenty gives occupancy probabilities that are almost all 0 or 1. Alex Diana suggested the regularised horseshoe on 25 November 2026 (comment on PR 41) and said the problem he means is this one: the upper tail, and with it the risk that a coefficient the data hardly inform runs to a very large value. Doug decided the same day to design it now, before a validation run is issued against the present prior, and chose the form set out here.

## The evidence so far

A prior-only simulation, with nothing fitted: 250 sites, 8 latent traits, 2,000 prior draws of one species' coefficients. The number reported is the standard deviation across sites of the embedding's contribution to that species' linear predictor, in logits. The scripts and their output are in `dev/simstudy/axis-slab-20261125/scratch/`. They are scratch code written before this design and they copy the prior instead of calling the package.

At 64 embedding columns:

- Present prior (`0.7 / K`, no slab): median 1.45, 95% point 27, 99% point 158.
- A cap of 2 on the axis scales alone: median 1.4, 95% point 11.5, 99% point 49. This fails. The embedding leftover scale `sigma_emb` is also half-Cauchy, and no cap on an axis scale reaches it.
- A cap of 2 on the axis scales and on the leftover's standard deviation (this design): median 1.3, 95% point 5.7, 99% point 9.5.
- A wide Normal factor on every coefficient part (the other form that works, not chosen): median 0.5, 95% point 1.4 at the present axis scale.
- Half-t with 3 degrees of freedom on the axis scales (the fallback named on 12 November): median 1.7, 95% point 13, 99% point 92.

At 128 columns the first three rows read 1.33 and 27.8, 1.26 and 11.0, and 1.26 and 5.7 (median and 95% point).

## Goals

- The prior's 95% point comes down to single figures, with the median still between 1 and 2 at 64 and at 128 columns, and with `0.7 / K` still the default axis scale.
- No coefficient part of an embedding row can have a prior standard deviation above the cap.
- With an infinite cap a grouped fit reproduces the present prior and sampler exactly, so the two can be compared on the same data and the evidence gathered for PRs 44 to 46 can still be reproduced.
- Fits under the fixed prior are unchanged bit for bit.

## Not goals

- The measured covariates. Their prior (`sigma_b` half-Cauchy, axis scale fixed at 1) stays as it is.
- The slow mixing of the axis scales of inactive columns, which sit near zero. A cap acts on large scales only.
- The division of a species' response between the shared part and the leftover, and the low state of `sigma_b`. Those are the questions of PRs 44 and 46.
- A prior on the cap. It is one fixed number.
- A scale per coefficient. The axis scale stays one number per embedding column.

## What must not happen

- A grouped fit with an infinite cap, or any fit under the fixed prior, drawing different numbers from before.
- A true response well inside the cap being shrunk measurably more than it is now.
- The cap being reported as a fix for a problem it does not address. The design says which problems those are, above.

## Approaches

Four were considered. Doug chose the second on 25 November 2026, on Claude's recommendation.

1. *A cap on the axis scales only.* The regularised horseshoe of Piironen and Vehtari (2017) applied to `lambda_k`, as first sketched to Alex. It fails on the numbers above: the 95% point stays near 11.
2. *A cap on the axis scales and on the leftover's standard deviation.* Chosen. The same construction, extended so that the cap also holds the leftover. The body of the prior stays where it is and only the tail is cut. Nothing depends on the number of species. The cost is that the draws of the axis scales and of `sigma_emb` lose their closed form, and the two embedding scale moves of PR 44 need new acceptance terms.
3. *A wide Normal factor on every coefficient part.* The construction of Nishimura and Suchard (arXiv 1911.02160), who keep every draw in closed form and prove convergence rates for the Polya-Gamma sampler in logistic regression. Not chosen, because here one scale is shared by many coefficients, and then the construction moves the body of the prior as well: the median falls to 0.5, restoring it needs an axis scale near `4 / K` (about six times weaker pull of inactive columns towards zero), and the implied prior on `sigma_emb` tightens as species are added (median 0.28, 0.17 and 0.06 at 30, 100 and 1,000 species). It is the fallback if the new draws of approach 2 mix badly.
4. *A half-t on the axis scales.* The fallback named on 12 November 2026. Not chosen: its tail is still wide, and it is not what Alex suggested.

## The prior, in words

There is one new number, the cap, written `c` below and set by `slab_scale` in `listPriors`. Its default is 2. It is in logits per standard deviation of an embedding column, because embedding columns are standardised.

- The axis scale `lambda_k` and the leftover scale `sigma_emb` keep their half-Cauchy priors, at `0.7 / K` and 1.
- What multiplies row k of `B` is no longer `lambda_k` but a capped version of it, the *effective axis scale*. For a small `lambda_k` the two are equal. As `lambda_k` grows the effective scale levels off at the cap.
- The community average, the trait effects and the latent-trait effects of column k use the effective axis scale where they used `lambda_k`.
- The leftover's standard deviation was `sigma_emb` times `lambda_k`. It is now that product, capped in the same way.

So well below the cap nothing changes, and every part of a row has a Normal prior with standard deviation at most `c`.

The cap is on each part, not on the whole coefficient. A species' coefficient on a column is the sum of the community average, the latent-trait part and the leftover (and a trait part when there are traits), so with every part at the cap its prior standard deviation is about `c` times the square root of 3, which is 3.5 at the default. A coefficient of 3.5 on a standardised column multiplies the odds by about 30 for one standard deviation of that column.

The cap is not far above the truths this project simulates. The reviewer's replica of the simulator's default embedding truth (3 blends, 10 active columns, 48 species) gave a largest single coefficient with a median of 3.5, and a strongest column whose root mean square across species had a median of 1.36 and was over 2 in about 2% of data sets. A column can reach such a spread under the cap: an effective scale of 1.36 needs an uncapped `lambda_k` of about 1.85, which costs little under a half-Cauchy. The cost rises steeply only as the effective scale nears the cap itself. Whether this shrinks real responses is what evidence step 3 measures, and it is why the default of 2 was put to Doug as a decision, recorded below.

## Mathematics

Notation as in the grouped prior's design: for embedding column k, `mu_k` and `g_k` are column k of `G`, `c_k` is column k of `C`, `e_ks` is the leftover of species s, q is the number of latent traits, S the number of species, a the axis scale and b the scale of the prior on `sigma_emb`.

**The prior.** Unchanged:

``` math
\lambda_k \sim C^{+}(a), \qquad \sigma_{\mathrm{emb}} \sim C^{+}(b)
```

New, with precisions adding:

``` math
\frac{1}{\tilde\lambda_k^{2}} = \frac{1}{\lambda_k^{2}} + \frac{1}{c^{2}}, \qquad \frac{1}{\omega_k^{2}} = \frac{1}{\sigma_{\mathrm{emb}}^{2}\lambda_k^{2}} + \frac{1}{c^{2}}
```

``` math
\mu_k,\ g_{hk} \sim N(0,\ \tilde\lambda_k^{2}), \qquad c_{jk} \sim N(0,\ \tilde\lambda_k^{2} / q), \qquad e_{ks} \sim N(0,\ \omega_k^{2})
```

The first of these is the regularised horseshoe of Piironen and Vehtari (2017, equation 2.8) with their global scale folded into a. The second is Claude's extension of it to the leftover. As c goes to infinity both reduce to the present prior. The leftover's cap is written on the product `sigma_emb * lambda_k` and not on the effective axis scale, so the leftover's standard deviation is at most c whatever `sigma_emb` is.

**The draw of an axis scale.** Write u for the log of `lambda_k`, and

``` math
P_g(u) = e^{-2u} + c^{-2}, \qquad P_e(u) = \frac{e^{-2u}}{\sigma_{\mathrm{emb}}^{2}} + c^{-2}
```

for the two precisions. Given everything else, the log density of u is, up to a constant,

``` math
f_k(u) = u - \log\!\left(1 + \frac{e^{2u}}{a^{2}}\right) + \frac{g_b + q}{2}\log P_g - \frac{P_g}{2}\left(\mathit{ss}_G + \mathit{ss}_C\right) + \frac{S}{2}\log P_e - \frac{P_e}{2}\,\mathit{ss}_R
```

- `ss_G` is the sum of squares of column k of `G` (its `g_b` rows, the community row included).
- `ss_C` is the sum of squares of column k of `C`, each divided by `c_base[k]^2`.
- `ss_R` is the sum of squares of column k of the leftovers.
- The leading u is the Jacobian of the log.

These are the three sums the present draw already forms. The density is proper: it falls as `-u` for large u (the half-Cauchy tail, since both precisions level off at the cap) and much faster for small u.

**The draw of the leftover scale.** Write v for the log of `sigma_emb`. Given everything else,

``` math
f_\sigma(v) = v - \log\!\left(1 + \frac{e^{2v}}{b^{2}}\right) + \sum_{k\,\in\,\text{embedding}} \left[ \frac{S}{2}\log P_{e,k}(v) - \frac{P_{e,k}(v)}{2}\,\mathit{ss}_{R,k} \right], \qquad P_{e,k}(v) = \frac{e^{-2v}}{\lambda_k^{2}} + c^{-2}
```

**How they are drawn.** Neither density is a standard one. Each is drawn by univariate slice sampling (Neal 2003) on the log scale, with stepping out and shrinkage, started from the current value with an initial width of 1. A slice draw leaves the conditional exactly invariant and needs no tuning. The order in the sweep is the present one: `sigma_b` (unchanged, closed form), then `sigma_emb`, then each axis scale given the new `sigma_emb`. During the staged burn-in the draw of `sigma_emb` is skipped, as now.

**The axis scale move.** It proposes `lambda_k' = lambda_k * exp(delta)`. At present it multiplies the whole of row k by `exp(delta)`. Under the cap it multiplies each part by the ratio of its own prior standard deviations, new over old:

``` math
\rho_g = \tilde\lambda_k' / \tilde\lambda_k, \qquad \rho_e = \omega_k' / \omega_k
```

Column k of `G` and of `C` is multiplied by `rho_g`, and column k of the leftovers by `rho_e`. Each part keeps its value in units of its own prior standard deviation, so the prior density of the scaled parts at the new scale, times the Jacobian of the scaling, equals the prior density of the old parts at the old scale, exactly as in the present move. Row k of `B` moves by

``` math
d_s = (\rho_g - 1)\,(B_{ks} - e_{ks}) + (\rho_e - 1)\,e_{ks}
```

and the log acceptance ratio is

``` math
\mathit{lin} - \tfrac{1}{2}\mathit{quad} + \delta - \log\!\left(1 + \lambda_k'^{2}/a^{2}\right) + \log\!\left(1 + \lambda_k^{2}/a^{2}\right)
```

``` math
\mathit{lin} = \sum_{i,s} (\kappa_{is} - \Omega_{is}\eta_{is})\, x_{ik}\, d_s, \qquad \mathit{quad} = \sum_{i,s} \Omega_{is}\, x_{ik}^{2}\, d_s^{2}
```

with the likelihood conditioned on the Polya-Gamma weights as now. With an infinite cap `rho_g` and `rho_e` both equal `exp(delta)` and this is the present move.

**The leftover scale move.** It proposes `sigma_emb' = sigma_emb * exp(delta)` and multiplies column k of the embedding leftovers by `rho_k = omega_k' / omega_k`, one factor per column. With `D` the sites by species matrix `X_emb diag(rho - 1) R_emb'`, the log acceptance ratio is

``` math
\sum (\kappa - \Omega\eta)\, D - \tfrac{1}{2}\sum \Omega\, D^{2} + \delta - \log\!\left(1 + \sigma_{\mathrm{emb}}'^{2}/b^{2}\right) + \log\!\left(1 + \sigma_{\mathrm{emb}}^{2}/b^{2}\right)
```

With an infinite cap every `rho_k` equals `exp(delta)` and this is the present move. The `sigma_b` move touches measured columns only and does not change.

## Where it goes

All in R. The C++ core already takes one prior standard deviation per covariate and does not change.

- `read_coef_prior()` reads `slab_scale`: a positive number or `Inf`, default 2, an error without `embeddingVars`. It is stored in the prior's record.
- `coef_prior_sds()` returns the capped standard deviations for the embedding columns: `omega_k` for the leftovers, the effective axis scale for `G`, and `c_base` times it for `C`. This one change carries the cap into the joint coefficient draw, the range step, `sample_GC()` and `sample_A()`.
- `sample_coef_scales()` uses the two slice draws for `sigma_emb` and the axis scales when the cap is finite.
- `axis_scale_move()` and `sigma_emb_move()` take the forms above when the cap is finite. The axis move needs the leftovers, which it is not passed today.
- `coef_C_start()` uses the effective axis scale. At the starting values the difference is negligible.
- The chain stores the effective axis scale in `lambda_output`, where it stored `lambda_k`. Every reader of `lambda_output` (`returnCovariateScales()`, the convergence summary, `returnGroupedConvergence()`) wants the number that multiplies the row. Above the cap the uncapped `lambda_k` is not informed by the data and wanders in its Cauchy tail, so it is not stored and not checked. `sigmaemb_output` stays `sigma_emb`.
- The `binary`, `occupancy` and `two_stage` model types get the cap by default. In a `continuous` fit the cap would be in units of the response, not logits, so there the default is `Inf` and a cap applies only if the user sets one. The scale moves stay in binary fits only, as now.

## What changes for users

- One new setting, `slab_scale` in `listPriors`, default 2 (`Inf` in a `continuous` fit). `Inf` gives the prior of PR 46.
- `returnCovariateScales()` reports the effective axis scale in its `axis` rows.
- The help page `?"grouped-prior"` gains the setting, a sentence on what the cap does, and a corrected *Status* paragraph on the prior's tail. `NEWS.md` gains an entry under the development heading.
- Nothing changes for a fit without `embeddingVars`.

## Tests

In a new file `tests/testthat/test-axis-slab.R`, and by extending the existing tests of the two moves.

- The capped standard deviations: equal to the present ones with an infinite cap, never above the cap, and equal to the formulas above at a few hand-computed values.
- The slice draws sample their conditionals: for a small state, many draws of an axis scale and of `sigma_emb` agree with the conditional found by numerical integration of `f_k` and `f_sigma`.
- Each move accepts on the exact ratio of the target density: the log ratio it returns equals the difference of the full log target at the proposed and current states, plus the log Jacobian of the scaling. Under the cap that Jacobian is `(g_b + q) log(rho_g) + S log(rho_e)` for the axis move and S times the sum over columns of `log(rho_k)` for the leftover scale move, where the existing tests add the number of scaled coefficients times `delta`.
- An accepted axis move changes column k of `G`, `C` and the leftovers by the two factors and nothing else. An accepted leftover scale move changes the embedding leftovers by one factor per column and nothing else.
- With no data, the scale draws together with all three moves leave the prior unchanged: the existing prior-only test, run with a finite cap, comparing the draws of the effective scales with direct draws from the prior.
- `slab_scale` is validated, and is an error without `embeddingVars` or under the fixed prior.
- Identity: with `slab_scale = Inf`, grouped fits of model types binary, occupancy and two-stage reproduce the draws of a03a1b8 exactly, and fixed-prior fits reproduce their reference fits exactly.
- `lambda_output` never exceeds the cap in a fit with a finite cap.
- The whole suite, with `NOT_CRAN=true`.

## Evidence before it is kept

In this order. The reading rules are fixed here, before any run.

1. *The prior check, with the package's own functions.* 250 sites, 8 latent traits, 2,000 prior draws, at 64 and at 128 embedding columns, default settings. It **passes** if the median is between 1 and 2 and the 95% point is under 8 at both column counts. The same script run with `slab_scale = Inf` is reported beside it and should give a 95% point near 27.
2. *A rare species that one column separates.* One simulated binary data set, small on purpose so that one species can drive a column's axis scale. It is built by a new script, not by `simulateOccJSDMData()`, with seed 3001:
    - 150 sites and 10 embedding columns, drawn as independent standard Normals and then standardised, so that column 1 is the only direction that separates.
    - 6 species. Species 1 is present at exactly the 6 sites with the largest values of column 1 and absent elsewhere. Species 2 to 6 are present independently at each site with probability 0.3 and respond to nothing.
    - The fit: model type `binary`, all 10 columns in `embeddingVars`, no measured covariate, no traits, `n_lattrait = 2`, `n_factors = 2`, 4 chains of 5,000 + 5,000, fit seed 3101, once with the default cap and once with `slab_scale = Inf`.
    - The number read is the 99% point, over the kept draws of all chains, of the absolute value of species 1's coefficient on column 1. The largest draw is reported beside it.
    - With the cap it must be under 14, which is four times the largest prior standard deviation the whole coefficient can have (3.5). Otherwise the cap **fails** at its purpose, and the first thing to look at is which part of the coefficient carries the excess.
    - Without the cap the same number is reported under the same limit: over 14 reads "the tail did harm here", and under 14 reads "no harm seen in this data set".
3. *The same simulated data with and without the cap.* The five data sets of the `sigma_b` move's simulation (256 sites, 48 species, 5 measured covariates, 64 embedding columns, 2 traits, 4 latent traits, binary, data seeds 1001 to 1005 and fit seeds 2001 to 2005), fitted from one library with `slab_scale = 2` and with `Inf`, 4 chains of 5,000 + 5,000. The simulator returns the whole true `B` and not its parts, so the truth is described by two numbers per data set, reported before fitting: the largest absolute true coefficient on an embedding column, and the largest root mean square across species of a column's true coefficients. Four readings:
    - *Recovery.* The root mean squared error of the posterior mean of the embedding's contribution to the linear predictor, against the truth, per data set. **No harm** if with the cap it is at most 10% larger than without in at least 4 data sets and at most 25% larger in every one. Otherwise **harm**.
    - *Shrinkage of the strongest columns.* In each data set, take the 5 embedding columns whose true coefficients have the largest root mean square across species, and regress the posterior means of those coefficients on their true values through the origin. **No extra shrinkage** if the slope with the cap is at least 0.95 of the slope without in at least 4 data sets. Otherwise **extra shrinkage**, which is the thing *What must not happen* names.
    - *Convergence.* A fit fails if `returnGroupedConvergence()` flags `rhat_measured`, `gradients` or `disagreement`. These are the targets of the version 2 validation contract, under which no scale and no chain state fails a fit. **No worse** if the number of failing data sets with the cap is at most the number without plus 1. The earlier record of these data sets suggests the count without the cap is not zero.
    - *Time.* Minutes per fit, with the power state recorded. **Slower** if the mean with the cap is more than 10% above the mean without.
4. *A replay on real sites,* by the AlphaEco project, as was done for the `sigma_b` move. It is not part of this build. Whether to ask for it, and moving the pin of the version 2 validation contract, are Doug's decisions after steps 1 to 3 are read.

If step 1 or 2 fails, or step 3 reads "harm" or "extra shrinkage", the cap is not kept as the default and the choice of what next (a different cap, approach 3, or the present prior) is Doug's.

## What to expect, and what not to

These are Claude's expectations from reasoning, not evidence.

- Step 1 should pass with room to spare, since the scratch simulation gave 95% points of 5.6 to 5.8.
- Step 2 may show that the present prior was not running away. An axis scale is shared by every species, so one species cannot inflate it without making the other species' leftovers implausible. This is a protection the textbook horseshoe, with a scale per coefficient, does not have. The data set of step 2 is small on purpose so that the protection is weak. If the uncapped fit stays bounded even there, the case for the cap rests on the prior check alone, and Doug and Alex should know that.
- Step 3 is a harder test than Claude first thought, because the strongest simulated columns are not far under the cap. Claude still expects little difference, since a column's spread of 1.4 is reachable at small cost, but this is the reading most likely to go against the default of 2.
- The cap will not make the axis scales of inactive columns converge faster.

## Decisions, and who made them

- Alex Diana suggested the regularised horseshoe and said it is for the prior's upper tail (25 November 2026).
- Doug decided to design it now and chose approach 2, "cap the scales" (25 November 2026).
- Claude's calls, shown to Doug in the conversation before he chose, and open to correction: the cap is fixed, default 2, and settable; it covers embedding columns only; an infinite cap reproduces the present prior and sampler; the work lives on its own branch stacked on PR 46, with a draft PR of its own.
- Claude's calls made while writing this, which Doug has not seen before: the leftover's cap is on the product of `sigma_emb` and `lambda_k`; the draws are slice draws; the two moves scale each part by the ratio of its own standard deviations; `lambda_output` stores the effective axis scale; a `continuous` fit has no cap unless one is set; and every reading rule under *Evidence before it is kept*.
- A reviewer with no part in the conversation read the first draft on 25 November 2026 as a builder would. It found no error in the mathematics, which Claude had also checked numerically against the full target density on a small random state. Its questions changed four things: the data set of step 2 is now specified and its limit raised from 10 to 14, because 10 left no margin over the prior at the cap; step 3 says how the truth is compared with the cap and gains the shrinkage reading; a failing fit in step 3 is defined; and the `continuous` model type is treated apart. It also showed that the simulated truths are nearer the cap than the first draft said, which is now stated under *The prior, in words*.

## The cap's default: decided

Doug chose 2 on 25 November 2026 ("Yes, write the plan"), after being shown the two options below. The default of 2 was Claude's proposal before the reviewer's finding about the simulated truths. The alternative was 3.

- *A cap of 2.* Prior median 1.3 and 95% point 5.7 at 64 columns (1.26 and 5.7 at 128). The whole coefficient's largest prior standard deviation is 3.5. The strongest simulated column, with a spread of about 1.4, sits at 70% of the cap.
- *A cap of 3.* Prior median 1.3 and 95% point 7.0 at 64 columns (1.3 and 7.8 at 128). The whole coefficient's largest prior standard deviation is 5.2. The strongest simulated column sits under half the cap. The pass mark of step 1 would become a 95% point under 10, and the limit of step 2 would become 21.

Claude recommended 2, because it cuts more of the tail and step 3 tests directly whether it shrinks the strongest columns. If step 3 reads "extra shrinkage", 3 is the first thing to try.

## Left to the builder

- Whether the finite-cap and infinite-cap paths share code or branch, so long as the identity test holds bit for bit.
- Whether the axis scales are slice-drawn one at a time or all at once in vector form. The requirement is that the scale draws add under 5% to the time of a sweep at 64 columns, measured and reported.
- The limit on stepping out in the slice draw, and how a density value that is not finite is treated.
- The layout of the study folder and its scripts, following `dev/simstudy/sigma-b-move-20261122/`.

## Risks

- *Two paths.* Keeping the present sampler reachable with an infinite cap means two sets of scale draws and two forms of each move. If the cap is kept, the infinite path could later be retired, which is a separate decision.
- *A real response above the cap.* A community that truly responds by more than about 2 logits per standard deviation of one column is shrunk. The setting is there for that case, and step 3 reports the largest true part.
- *The meaning of `sigma_emb` above the cap.* Below the cap it is the ratio of the leftover's scale to the shared part's. Where a column sits at the cap that reading is lost.
- *Evidence gathered without the cap.* The mixing trials of PRs 44 and 46 and the gradient comparisons of PR 45 were run under the present prior. Step 3 covers five data sets. Only a replay (step 4) would show that those results carry over.
- *The comparison in step 3 is noisy.* Chains of this model are known to settle in different states, which moves the error of a fit by more than a cap is likely to. The 25% bound is wide for that reason, and a "harm" reading needs a look at the chains before it is believed.

## References

- Piironen, J. and Vehtari, A. (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. *Electronic Journal of Statistics* 11(2): 5018-5051. The regularised horseshoe is their equation 2.8. They name logistic regression with separable data as a case where the plain horseshoe does harm.
- Nishimura, A. and Suchard, M. A. Shrinkage with shrunken shoulders: Gibbs sampling shrinkage model posteriors with guaranteed convergence rates. arXiv 1911.02160. Journal details not checked.
- Neal, R. M. (2003). Slice sampling. *Annals of Statistics* 31(3): 705-767. Cited from memory, not checked.
