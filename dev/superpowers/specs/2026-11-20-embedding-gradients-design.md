# Gradients of the embedding response, and their comparison with covariates

Date: 20 November 2026. Status: design written from a conversation with Doug on 20 November 2026 and approved by him the same day; implemented on this branch on 20 November 2026 (plan in `dev/superpowers/plans/2026-11-20-embedding-gradients.md`). Branch `claude/embedding-gradients`, stacked on `claude/grouped-prior-part-2` (draft PR 44), which is stacked on `claude/grouped-covariate-prior` (PR 41).

## Purpose

Give a user who has fitted `runOccJSDM()` with an embedding two functions for reading the fit afterwards.

- `returnEmbeddingGradients()` returns the leading gradients of the species' whole response to the embedding: a score for every site, a position for every species, and the share of the response each gradient holds.
- `compareEmbeddingGradients()` compares those gradients with other covariates known for the same sites, in the main use covariates that were not in the model, and returns what the covariates leave unexplained as "unmeasured gradients".

## Why

Doug wants to use a fit for interpretation as well as prediction: relate what the species respond to in the embedding to covariates measured afterwards, to generate hypotheses about what drives the species' distributions. The embedding columns have no individual meaning, so the only readable object is the low-dimensional set of gradients the species share.

The package already has `returnEmbeddingBlends()`, but it works from the shared part of the response only (`A %*% C` on the embedding columns). The simulations of 20 November 2026 (see Evidence) showed that the division of a species' response into shared part and leftover is poorly determined, and that at weak signal a fit can put most of the response in the leftovers although the species share gradients. The leading directions of the whole response were recovered well in the same fits. So interpretation should start from the whole response.

The second function exists because the natural mistake is to correlate one gradient with one covariate. Single gradients are defined only up to rotation, so that correlation depends on an arbitrary choice. The comparison must be made with the gradients as a set, and the function makes the correct analysis the easy one.

## The picture to keep in mind

Doug's summary of 20 November: this is a JSDM inside a JSDM. The main model already estimates hidden environmental covariates (`n_factors`) from what the species share after the covariates are removed. The embedding gives more of them.

The two kinds differ in one respect that matters.

- The hidden environmental covariates are an unconstrained ordination: each site's score is a free parameter, learned only from the species at that site. They cannot be computed at a site with no species data.
- The embedding gradients are a constrained ordination: each site's score is a weighted sum of that site's embedding columns. They can be computed, and mapped, wherever the embedding exists.

A model with both is what the ordination literature calls concurrent ordination (constrained and unconstrained axes fitted together). Two consequences follow. First, a hidden environmental covariate can partly duplicate an embedding gradient, most likely when the prior shrinks the embedding response and the factors absorb what is left. Second, the comparison function can test this, because the factors' site scores are just another table of per-site numbers: pass them in as the covariates. No extra function is needed.

## Covariates that were in the model

Doug's question of 20 November: are the embedding gradients not fitted to what is left after the measured covariates are removed, so that comparing them with those covariates is empty?

They are not fitted to a residual. The measured covariates and the embedding columns are fitted together, so each coefficient is a partial effect: an embedding coefficient is the effect of that column with the measured covariates held fixed. That removes the measured covariates from the coefficients, not from the site scores. A gradient's site score is a weighted sum of embedding columns, and embedding columns are correlated with the measured covariates (the satellite sees elevation), so a gradient can be strongly correlated with a covariate that was in the model. Only the hidden environmental covariates come close to the residual picture, and even they are fitted jointly.

So the two kinds of covariate are read differently.

- **Not in the model (the main use).** A match suggests the gradient may be standing in for that covariate. This is the hypothesis-generating use.
- **In the model.** A match shows overlap between the covariate and what the embedding carries: species respond, beyond that covariate's own coefficient, along a direction of habitat correlated with it. It can mean a humped response to the covariate, or a correlated habitat feature. It also warns that the division of the response between the covariate and the embedding is uncertain, as in "Reading a management covariate".

To ask what a new covariate adds beyond those in the model, put both in the table and read the new covariate's `unique` value. To get gradients that are free of the model's covariates, pass those covariates and read the unmeasured gradients. No extra argument is needed for either.

## Goals

- Return gradients computed from the whole embedding response, per chain and pooled, with the number of gradients agreeing with what `returnGroupedConvergence()` reports.
- Make every number that is presented as a finding independent of how the gradients are rotated.
- Match sites by name, never by row order.
- Give the unmeasured gradients as site scores that can be mapped and species positions that can be read.
- Give Alex a mathematical statement of each quantity (below), and users a help section in plain language.

## Non-goals

- Uncertainty intervals for single gradients from per-draw decompositions. Matching gradients across draws is the rotation problem again, and it has not been validated.
- p-values or any other test. Sites are spatially autocorrelated, and the purpose is hypothesis generation.
- An argument for scoring new sites. The weights are returned, so scores at new sites can be computed by hand; a `newdata` argument can follow if maps of gradients are wanted.
- Plotting functions.
- Any change to `returnEmbeddingBlends()` beyond one sentence in its help pointing to the new function.
- Nonlinear comparisons beyond optional squared terms.
- A simulation study of the comparison function. Under Doug's rule of 17 November 2026 such a study would be designed here and run in the AlphaEco project; whether to ask for one is his decision.

## Anti-goals

These would make the work a failure even if the code ran.

- A result that depends on the rotation of the gradients and is not labelled as depending on it.
- Covariates silently matched to the wrong sites.
- Help text or output that invites reading one gradient as one driver, or a match as evidence of a cause.
- Any change to what `returnGroupedConvergence()` returns, to the sampler, to the prior or to a default.

## Constraints

- No change to the sampler, the prior, any default, or the output of `returnGroupedConvergence()`. The existing tests must pass unmodified.
- `B_output` is the largest object a fit holds. Chain means of its embedding rows must be computed without copying the whole array.
- The gradients must be computed from the same effect array the convergence check uses, through shared code, so the two cannot drift apart.
- Base R only inside the functions; no new dependency.
- Project rules as in `AGENTS.md` and `CLAUDE.md`: install from a git export, never `devtools::load_all()`, R 4.6, markdown written one paragraph per line with no pipe tables and no em-dashes.
- No dataset content in the repository. Tests use constructed or simulated data only.

## Approaches considered

- **From posterior means (chosen).** Each chain's posterior mean of the embedding coefficients gives that chain's effect at the fitted sites; its singular value decomposition gives the gradients. The comparison is ordinary least squares. This is the computation the convergence check already makes and the one the simulations of 20 November validated. Per-chain and pooled values are both returned, so disagreement between chains is visible.
- **A decomposition per posterior draw, with intervals (rejected).** It gives intervals for single gradients, but those are not identified, so the draws must first be aligned to one another, which is the rotation problem again. It is slow with thousands of species and has not been validated.
- **Covariates inside the model as constraints on the gradients (rejected).** That is a different model, not a way of reading this one, and the covariates of interest are measured after the fit.

Doug chose the first on 20 November 2026, and chose to include the comparison as a second function.

## Mathematical description

Notation. The subscript `emb` marks the embedding block, as in the package's code (`X_emb`, `sigma_emb`, `is_emb`). `X_e` is avoided because in the models that motivated this work it means the measured environmental covariates, and `X_a` because it names one embedding, AlphaEarth. There are `n` fitted sites, `S` species, `K` embedding columns, `T` chains and `M` kept draws per chain. `X_emb` is the `n` by `K` block of embedding columns of the standardised occupancy design matrix `X_psi`, so each column has mean 0 over the fitted sites. `B_emb` is the `K` by `S` block of embedding rows of the occupancy coefficients `B`. `'` is transpose and `||.||` the Frobenius norm.

The linear predictor of species `s` at site `i` contains the covariate term `x_i' b_s` and the hidden-environmental-covariate term `u_i' l_s`. Under the grouped prior `B_emb` is a centre (community average, traits, and blends `A C`) plus a leftover; the exact parametrisation is in the grouped-prior design of 11 November 2026. The functions here use `B_emb` only, never the centre or the leftover separately.

### Gradients

```text
Bbar_t  = (1/M) sum_m B_emb[t, m]          K x S   chain t's posterior mean
E_t     = X_emb Bbar_t                     n x S   chain t's embedding effect
E       = (1/T) sum_t E_t                  n x S   pooled effect
E       = U D V'                           thin SVD, d_1 >= d_2 >= ... >= 0

r       = smallest r with  sum_{k<=r} d_k^2 / sum_k d_k^2  >=  0.90     (unless the user sets r)

Z       = sqrt(n) U_r                      n x r   site scores
L       = V_r D_r / sqrt(n)                S x r   species loadings
W       = Bbar V_r D_r^{-1} sqrt(n)        K x r   weights,  Bbar = (1/T) sum_t Bbar_t

share_k = d_k^2 / sum_j d_j^2
```

Properties.

- `Z L' = U_r D_r V_r'`, the best rank-`r` approximation of `E`.
- `Z = X_emb W`, so the site scores are weighted sums of the embedding columns.
- Each column of `Z` has mean 0 and mean square 1 over the fitted sites, because the columns of `X_emb`, and so of `E`, have mean 0. A species' loading is therefore the change in its linear predictor per one standard deviation of the gradient.
- These properties need `d_r` above 0. Because the columns of `X_emb` have mean 0, the rank of `E` is at most the smallest of `n - 1`, `S` and `K`, and less when embedding columns are collinear, so `r` is limited to the numerical rank of `E`.
- `0.90` is `GROUPED_LIMITS$hold`, and `E_t` and `E` are the arrays `grouped_convergence()` builds, so whenever the convergence check assesses the effect, the default `r` equals its `n_gradients`.
- For any orthogonal `r` by `r` matrix `Q`, `(Z Q)(L Q)' = Z L'`. The pair `(Z, L)` is one representative of an equivalence class; the function returns the principal-axis representative, with signs fixed.

Agreement between chains is the convergence check's own measure. With `P_t` the projection onto the top `r` left singular vectors of `E_t`,

```text
agreement[a, b] = || P_b E_a ||^2 / || P_a E_a ||^2 ,      agreement = min over a != b
```

It is 1 when two chains have the same gradients and cannot exceed 1.

### Comparison with covariates

Let `I` be the `m` sites used (in the fit, present in the covariate table, no missing covariate value). Subscript `I` selects those rows. `F` is the covariate design matrix with an intercept column, `H` the projection onto its column space (`F (F'F)^{-1} F'` when `F` has full column rank; computed by QR so that collinear covariates are allowed), `q` the rank of `F` minus 1, and `J` the centring matrix. `E_r = Z L'` is the rank-`r` effect.

Check 1, covariate from gradients. For a covariate `y` (an `m`-vector) and `k` from 1 to `r`, with `H_k` the hat matrix of an intercept and the first `k` columns of `Z_I`:

```text
R2_k(y)     = 1 - || (I - H_k) y ||^2 / || J y ||^2
adjR2_k(y)  = 1 - (1 - R2_k) (m - 1) / (m - k - 1)
```

At `k = r` this is unchanged by any rotation of `Z`.

Check 2, gradients from covariates. This is the redundancy-analysis statistic of the rank-`r` effect on the covariates:

```text
explained     = 1 - || (I - H) E_r[I, ] ||^2 / || J E_r[I, ] ||^2
adj explained = 1 - (1 - explained) (m - 1) / (m - q - 1)
```

It depends on `(Z, L)` only through `Z L'`, so it is unchanged by rotation. It is computed without forming `E_r`. Let `N = (L'L)^{1/2}` be the symmetric square root and `G = Z_I N`. Then `O = L N^{-1}` has orthonormal columns and `E_r[I, ] = G O'`, so the same ratio holds with `G` in place of `E_r[I, ]`. This holds for any representative `(Z, L)` with `L` of full column rank. For the principal-axis representative the function returns, the columns of `L` are orthogonal, `N` is diagonal with entries `g_k = || L[, k] ||`, and `G[, k] = g_k Z_I[, k]`.

- Per covariate `j`: `alone_j` is `explained` with `F` reduced to the intercept and covariate `j`'s columns; `unique_j` is `explained` minus `explained` with covariate `j`'s columns removed.
- Per gradient `k`, for the principal-axis representative only: `R2_k`, the R-squared of `G[, k]` on `F`, and `share_k = || J G[, k] ||^2 / || J G ||^2`, the gradient's share of the rank-`r` effect at the sites used. Then `explained = sum_k share_k R2_k`. These depend on the principal-axis convention and are labelled so.

Check 3, unmeasured gradients. Decompose what the covariates leave:

```text
(I - H) G  = P Lambda Q'                   thin SVD
Z_u        = sqrt(m) P                     m x r   unmeasured site scores
L_u        = O Q Lambda / sqrt(m)          S x r   unmeasured species loadings,  O = L (L'L)^{-1/2}
share_u[j] = Lambda_j^2 / || J G ||^2
```

Properties: `Z_u L_u' = (I - H) E_r[I, ]`, and `sum_j share_u[j] = 1 - explained`. `(I - H) G` has rank at most the smaller of `r` and `m - q - 1`. For every `j` with `Lambda_j` above 0, column `j` of `Z_u` has mean 0 and is uncorrelated with every column of `F`; the singular vectors of zero singular values are arbitrary, so those components are dropped and only the ones with `Lambda_j` numerically above 0 are returned. This is the principal component analysis of the residual of the redundancy analysis in check 2.

Per chain, checks 1 and 2 are repeated with chain `t`'s own `(Z, L)` from `E_t`. Check 3 is pooled only, because single unmeasured gradients cannot be matched across chains.

## Design: `returnEmbeddingGradients(fitModel, n_gradients = NULL)`

- **Accepts** any fit for which `returnEmbeddingShare()` works: the grouped prior with at least one embedding column. It does not need latent traits. Otherwise it stops with a message in the style of the neighbouring functions.
- **`n_gradients`**: `NULL` for the 90% rule, or one whole number from 1 to the numerical rank of the pooled effect, which is at most the smallest of `n - 1`, `S` and `K`. A larger value stops with a message giving the rank.
- **When the convergence check does not assess the effect.** The check declines when the effect is small, with one chain, or with non-finite values. This function behaves differently on purpose (Claude's call): with a small effect it still returns the gradients by the same 90% rule, with the warning below, because the user asked for them; with one chain it returns them with `agreement` as `NA`; with non-finite values it stops. So the equality with the check's `n_gradients` and `agreement` holds whenever the check assesses the effect, and the help says so.
- **Returns** a list.
    - `n_gradients`, and `rule`, which is `"hold 0.90"` or `"set"`.
    - `site_scores`: `Z`, with site names from `infos$siteNames` as row names and `gradient1`, `gradient2`, and so on as column names.
    - `species_loadings`: `L`, with species names as row names.
    - `weights`: `W`, with the embedding column names as row names.
    - `share`: a data frame with columns `gradient`, `share` and `cumulative`, with one row for every direction up to the numerical rank of the pooled effect, not only the `r` returned, so the reader can see what one more gradient would hold.
    - `size`: the convergence check's size of the effect (the median over species of the standard deviation across sites of `E`).
    - `agreement` and `agreement_pairs`: as defined above, at the `r` in use, whether that came from the rule or from the user; `NA` with one chain. With a user-set `r` the value can differ from the convergence check's, which uses its own `r`, and the warning below then names the `r` it refers to.
    - `by_chain`: one entry per chain with that chain's `site_scores`, `species_loadings` and `share`, from `E_t` with the same `r`. A chain whose own effect has rank below `r` gets `NULL`, with a warning.
- **Signs.** Each gradient is flipped so that its weight of largest absolute value is positive, the convention `returnEmbeddingBlends()` uses. Scores, loadings and weights flip together. Per-chain gradients are signed by the same rule on their own weights and are not matched to the pooled ones.
- **Warnings.** One when `size` is below `GROUPED_LIMITS$small_effect`, saying the effect is small and its gradients are mostly noise. One when `agreement` is below `GROUPED_LIMITS$gradients`, saying the chains found different gradients and pointing to `returnGroupedConvergence()`. There is no warning for the `state` or R-hat rows: the evidence is that those disagreements do not reach the gradients, and the help says so.

## Design: `compareEmbeddingGradients(gradients, covariates, site_col = "Site", squared = FALSE)`

- **`gradients`**: the list from `returnEmbeddingGradients()`.
- **`covariates`**: a data frame with one row per site. The column named by `site_col` holds site names, compared as character with the row names of `gradients$site_scores`. Rows are matched by name in any order. Every other column is a covariate.
- **Site handling.** A missing `site_col` or a repeated site name is an error. Rows whose site is not in the fit are ignored, fitted sites absent from the table are left out, and sites with a missing value in any covariate are left out; each of these is reported in a message with its count. All three checks use the same sites. It is an error if no site matches, or if `m` is less than `q + 3` or less than `r + 2`.
- **Covariate handling.** Numeric and logical columns are used as numbers. Factor and character columns are categorical and enter `F` as indicator columns. A covariate that is constant over the sites used is an error naming it. With `squared = TRUE`, every numeric covariate with more than two distinct values also contributes the square of its centred values to `F`; the square counts as part of that covariate in `alone` and `unique`.
- **Returns** a list.
    - `covariate_from_gradients` (check 1): a data frame with columns `covariate`, `chain` (`"pooled"` or the chain number), `n_gradients` (1 to `r`), `r_squared` and `adj_r_squared`. A categorical covariate has one row set per level indicator, or one for a two-level covariate. Squared terms do not enter this check.
    - `gradients_from_covariates` (check 2): a list of three data frames. `overall` has columns `chain`, `explained` and `adj_explained`. `by_covariate` (pooled) has columns `covariate`, `alone` and `unique`. `by_gradient` (pooled) has columns `gradient`, `share`, `r_squared` and `adj_r_squared`, where `share` is the gradient's share of the rank-`r` effect at the sites used (`share_k` above), not the `share` in the gradients list, so that the shares times the R-squared values sum to `explained`.
    - `unmeasured` (check 3): a list with `site_scores` (`Z_u`, rows named by the sites used, columns `unmeasured1` and so on), `species_loadings` (`L_u`) and `share` (columns `gradient`, `share`, `cumulative`, as shares of the rank-`r` effect at the sites used). Components with no share are dropped, so there can be fewer than `r`, and none when the covariates account for everything; the matrices then have no columns and `share` no rows.
    - `sites`: the number used and the names left out, by reason.
- **No tests, no intervals.** The help explains why.

## Help text

- A new section, "Interpreting the gradients", on the grouped-prior help page after "Reading a management covariate". Its text is the appendix at the end of this document, which is the draft Doug saw on 20 November 2026 revised to name the two functions and to add the constrained and unconstrained picture above. The builder turns the appendix into roxygen and may adjust wording to fit the page, but keeps every caution.
- It states one use explicitly: pass the hidden environmental covariates' site scores as `covariates`, from the same fit to see whether a factor duplicates an embedding gradient, or from a fit without the embedding to see which hidden environmental covariates the satellite sees. `returnOrdinationScores()` returns an array of quantiles by sites by factors, not a per-site table, so the help gives the recipe: take the median slice and make a data frame with a `Site` column from `infos$siteNames`. Factors are themselves defined only up to rotation, so the reader is told to read `overall`, not single factors.
- The list of summarising functions on that page gains the two new ones.
- The help of `returnEmbeddingBlends()` gains one sentence: it works from the shared part only, and `returnEmbeddingGradients()` works from the whole response.

## Tests

On constructed fit objects with a known low-rank effect (no sampling):

- `site_scores %*% t(species_loadings)` equals the best rank-`r` approximation of the pooled effect, and `X_emb %*% weights` equals `site_scores`.
- The 90% rule picks the constructed rank, and `n_gradients` overrides it; invalid values stop.
- Identical chains give `agreement` 1; chains built with different gradients give a value below the limit and the warning; a tiny effect gives the small-effect warning.
- Signs follow the stated rule.

Against the convergence check:

- On an effect the check assesses, `n_gradients` and `agreement` equal `returnGroupedConvergence()`'s `effect$n_gradients` and `effect$agreement`. The cached grouped fit the existing tests use has 2 chains of 20 + 20 iterations and may have a small effect; if the check does not assess it, this test uses a constructed fit object large enough to be assessed.
- On a small effect the function still returns gradients, with the warning, while the check reports "not assessed".
- A user-set `n_gradients` above the rank stops.

For the comparison, on constructed gradients:

- A covariate that is an exact linear combination of the site scores has R-squared 1 at `k = r`; an unrelated one has adjusted R-squared near 0.
- Covariates that span the gradients give `explained` 1 and no unmeasured gradient. Covariates that span one of two gradients leave exactly one unmeasured gradient, equal up to sign to the other, with the expected share.
- `sum(unmeasured$share$share)` equals `1 - explained`, the returned unmeasured site scores are uncorrelated with every covariate, and `by_gradient`'s shares times its R-squared values sum to `explained`.
- Rotating `site_scores` and `species_loadings` by the same orthogonal matrix leaves `overall`, `by_covariate`, the unmeasured shares and check 1 at `k = r` unchanged.
- Shuffling the rows of `covariates` changes nothing. Extra rows, absent sites and missing values are left out with the stated messages and counts. A missing `site_col`, repeated names, a constant covariate and too few sites each stop.
- Categorical covariates and `squared = TRUE` behave as described.

The existing suite passes unmodified.

## Evidence

From the AlphaEco simulation session, 20 November 2026, on occJSDM 2b9b890; 20 grouped fits (truths t1 at strong signal and t4z at weak signal, 10 replicates each, 2 chains of 5,000 + 5,000). The truth has three shared gradients.

- The top three site directions of the fitted whole effect predict the true gradients' site scores with R-squared 0.95 in t1 and 0.81 in t4z (0.96 and 0.79 at held-out sites). The embedding's own principal components reach 0.43 and 0.36.
- The nine t4z fits whose chains all ended in "leftovers dominate" still reach 0.80.
- The shared part alone matches the true blends (mean squared cosine, sites then species) at 0.94 and 0.76 in t1, but 0.55 and 0.29 in those nine fits.
- The two chains' top three site subspaces agree almost exactly (median 1.00, lowest 0.95), including fits whose chains divide the effect very differently.
- The species side is weaker than the site side: whole effect against truth 0.84 in t1 and 0.60 in t4z.
- The 90% rule gave 2 or 3 gradients in t1 and 3 or 4 in t4z.

Limits of that evidence: it does not show the grouped prior is needed for this, since the comparison was with the embedding's own components and not with fits under other priors; the simulated truths use 10 of 64 embedding columns, which earlier tests showed flatters this prior; one set of sites, one embedding, 48 species per fit. The comparison function itself has no simulation evidence; its tests are algebraic.

## Left to the builder

- Where the code lives (a new file is suggested) and how the effect array is shared with `grouped_convergence()`, provided that function's results do not change.
- The wording of messages and errors, in the style of the neighbouring functions.
- The sign rule for the unmeasured gradients, provided it is deterministic and documented.
- Version, `NEWS` and `TODO.md` entries, following `AGENTS.md`.

## Decisions that are Doug's

- Approval of this design.
- The two function names, which are Claude's proposals.
- Whether a simulation study of the comparison function is requested from the AlphaEco project.

Calls Claude made after a fresh reader checked this document on 20 November 2026, for Doug to confirm or change:

- With a small effect the gradients function still returns gradients, with a warning, although the convergence check declines to assess such a fit.
- With a user-set number of gradients, chain agreement is computed at that number, so it can differ from the convergence check's.
- Unmeasured gradients that hold nothing are dropped, so there can be fewer than the number of gradients, or none.
- In the per-gradient table of the comparison, a gradient's share is its share at the sites used.

## Appendix: text of the help section "Interpreting the gradients"

The embedding columns have no meaning one by one. What can be read is the small set of gradients that the species share in their responses to the embedding. Two gradients mark out a plane and three a block. Any rotation of them marks out the same plane or block and fits the data equally well, so the data determine the plane, not the particular directions a function returns. Read the gradients as a set, never one by one.

`returnEmbeddingGradients()` takes the gradients from each species' whole response to the embedding, the shared part plus the leftover. How a response is divided between those two is poorly determined, and at weak signal a fit can put most of it in the leftovers although the species share gradients. The leading directions of the whole response are far more stable. In simulations they recovered the true gradients across sites, and chains agreed on them, even in fits where the division was wrong and the scales failed their convergence checks. `returnEmbeddingBlends()` works from the shared part only, so rely on it only when the `state` row of `returnGroupedConvergence()` is 0.

These gradients are hidden environmental covariates of a second kind. Those set by `n_factors` are free scores for each site, learned only from the species found there. The embedding gradients are weighted sums of a site's embedding columns, so they can be computed and mapped wherever the embedding exists.

`compareEmbeddingGradients()` relates the gradients to covariates measured at the same sites. It makes three checks.

- Covariate from gradients: how well do the site scores of the first few gradients, taken together, predict each covariate?
- Gradients from covariates: how much of the gradients do all the covariates together account for, and how much does each covariate account for alone and uniquely?
- Unmeasured gradients: the leading directions of what the covariates leave unexplained. This is habitat variation that the embedding records and the species respond to, but that the covariates do not describe, which includes a response to a measured covariate that is not a straight line. Each has a score at every site used, which can be mapped for those sites, and a position per species.

The covariates can also be the site scores of the hidden environmental covariates, from this fit or from a fit without the embedding. The comparison then shows whether a factor duplicates an embedding gradient, or which hidden environmental covariates the satellite sees. Take the median slice of `returnOrdinationScores()` and add a `Site` column. Factors are also defined only up to rotation, so read the overall share, not single factors.

Cautions.

- The number of gradients is uncertain. Check that a conclusion holds with one more and one fewer.
- A match is a hypothesis about a driver, not evidence of one. Anything correlated with the covariate across your sites would match as well. The function gives no test and no interval, because sites near one another are not independent.
- The number of gradients is not the number of drivers. One covariate with a humped response can take up two gradients, and one gradient can be a mixture of several covariates. Set `squared = TRUE` when responses may be humped.
- Species' positions on the gradients are less well determined than the sites' scores, most of all at weak signal.
- A covariate that was in the model can still match the gradients, because the embedding columns are correlated with it. The model holds it fixed when estimating the embedding coefficients, but that does not make the gradients' site scores unrelated to it. Such a match shows overlap between the covariate and the embedding, not a new driver. To ask what a new covariate adds beyond the model's covariates, give both and read `unique`.
