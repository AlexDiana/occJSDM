# Embedding gradients Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `returnEmbeddingGradients()` and `compareEmbeddingGradients()`, which read a grouped-prior fit's embedding response as a set of gradients and compare them with per-site covariates.

**Architecture:** Both functions are post-fit summaries in a new file `R/embedding-gradients.R`. The gradients come from the same chain means and effect array the convergence check builds, through two small helpers extracted from `R/grouped-convergence.R` with no change in that check's results. The comparison is ordinary least squares done by QR, on the gradients list alone.

**Tech Stack:** Base R (`svd`, `qr`, `stats::model.matrix`), roxygen2 8.1.0, testthat edition 3.

**Spec:** `dev/superpowers/specs/2026-11-20-embedding-gradients-design.md`. Read it first, in particular "Mathematical description" and the two "Design" sections; this plan does not repeat the formulas.

## Global Constraints

- Branch `claude/embedding-gradients`, worktree `~/src/occJSDM/.worktrees/embedding-gradients`. Do not touch `claude/grouped-covariate-prior` (PR 41) or `claude/validation-contracts`.
- No change to the sampler, the prior, any default, or anything `returnGroupedConvergence()` returns. No existing test is edited; new tests may be appended to an existing test file.
- `B_output` is never copied whole. Chain means come from `blockRhatESS()`, as in the convergence check.
- Base R only inside the new functions. No new dependency.
- Never `devtools::load_all()` and never `R CMD build` in a checkout. R is 4.6. Install and test with the loop below, with `SP` set to a scratch directory outside the repository.
- Notation in comments and help: the embedding block is `X_emb` and its coefficients `B_emb`, never `X_e` or `X_a`.
- Markdown (`NEWS.md`, this plan, the spec): one paragraph per line, no pipe tables, no em-dashes, and in prose outside code write `-\>`, `\~`, `\<`, `\>`, `\|`. Re-read a markdown file before and after editing it.
- Help text is plain language in the style of the `grouped-prior` page. No dataset content anywhere; tests use constructed or simulated data.
- Commit messages end with `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>`. Do not push.

The loop, referred to below as "install", "document" and "test FILTER":

```bash
# install
mkdir -p $SP/build $SP/lib
git ls-files -co --exclude-standard -z | rsync -a --from0 --files-from=- ./ $SP/build/
R_PROFILE_USER=/dev/null R CMD INSTALL --no-docs --library=$SP/lib $SP/build > $SP/install.log 2>&1 || tail -30 $SP/install.log

# document (after install; then install again so the new NAMESPACE is in the library)
R_PROFILE_USER=/dev/null R_LIBS=$SP/lib Rscript -e 'roxygen2::roxygenise(".", load_code = roxygen2::load_installed)'

# test FILTER
(cd tests && R_PROFILE_USER=/dev/null R_LIBS=$SP/lib Rscript -e 'testthat::test_dir("testthat", package = "occJSDM", load_package = "installed", filter = "FILTER", stop_on_failure = TRUE)')
```

Before Task 1, commit the spec and this plan:

```bash
git add dev/superpowers/specs/2026-11-20-embedding-gradients-design.md dev/superpowers/plans/2026-11-20-embedding-gradients.md
git commit -m "Design and plan for gradients of the embedding response"
```

## Review Focus

Each line is pinned by a test in the task named.

- Site names that are numbers in the fit and in the covariate table: they must match as text (Task 3).
- One gradient, where a matrix can silently become a vector (Tasks 2 and 4).
- One chain: no agreement, and only "pooled" and "1" in the comparison (Tasks 2 and 3).
- Duplicated or collinear covariates: finite results, and `unique` of 0 for each of a duplicated pair (Task 4).
- An occupancy fit, where sites are sorted and have several rows of data each: one site score per site, named by `infos$siteNames` (Task 2).

Not testable in the suite, for the reviewer to read: that no code path subsets `B_output` by chain or by embedding rows as a whole.

---

### Task 1: Share the effect array and the agreement measure

**Files:**
- Modify: `R/grouped-convergence.R` (`grouped_effect_checks()` lines 76 to 90, `grouped_convergence()` lines 217 to 222)
- Create: `tests/testthat/test-embedding-gradients.R`

**Interfaces:**
- Produces: `embedding_effect(chain_mean, X, is_emb, S)` returning `list(B_mean, E)`, where `chain_mean` is `blockRhatESS()$chain_mean` (rows are coefficients, coefficient i of species j in row `i + (j - 1) * p`; columns are chains), `B_mean` is K by S by chains and `E` is `nrow(X)` by S by chains.
- Produces: `gradient_agreement(E, r)` returning a chains by chains matrix with `NA` on the diagonal and entry `[a, b]` equal to `sum(crossprod(U_b, E[, , a])^2) / sum(d_a[1:r]^2)`, where `U_b` holds chain b's top `r` left singular vectors and `d_a` chain a's singular values.

- [ ] **Step 1: Write the failing tests** in `tests/testthat/test-embedding-gradients.R`:

```r
# Gradients of the embedding response and their comparison with covariates
# (R/embedding-gradients.R; design in
# dev/superpowers/specs/2026-11-20-embedding-gradients-design.md).

test_that("gradient_agreement() gives the pairs grouped_effect_checks() reports", {
  set.seed(1); n <- 40L; S <- 12L
  U <- qr.Q(qr(matrix(rnorm(n * 3), n, 3))); V <- matrix(rnorm(S * 3), S, 3)
  E1 <- U[, 1:2] %*% diag(c(9, 5)) %*% t(V[, 1:2])
  E2 <- U[, 2:3] %*% diag(c(9, 5)) %*% t(V[, 1:2])
  E <- array(c(E1, E2), c(n, S, 2L))
  got <- occJSDM:::gradient_agreement(E, 2L)
  expect_identical(dim(got), c(2L, 2L)); expect_true(all(is.na(diag(got))))
  expect_equal(got, occJSDM:::grouped_effect_checks(E)$agreement_pairs)
  expect_equal(occJSDM:::gradient_agreement(array(c(E1, E1), c(n, S, 2L)), 2L)[1, 2], 1)
})

test_that("embedding_effect() rebuilds each chain's effect from the chain means", {
  set.seed(2); n <- 10L; p <- 5L; S <- 3L; is_emb <- c(FALSE, TRUE, TRUE, FALSE, TRUE)
  X <- matrix(rnorm(n * p), n, p); B1 <- matrix(rnorm(p * S), p, S); B2 <- matrix(rnorm(p * S), p, S)
  got <- occJSDM:::embedding_effect(cbind(as.vector(B1), as.vector(B2)), X, is_emb, S)
  expect_identical(dim(got$B_mean), c(3L, S, 2L)); expect_identical(dim(got$E), c(n, S, 2L))
  expect_equal(got$B_mean[, , 2], B2[is_emb, ]); expect_equal(got$E[, , 1], X[, is_emb] %*% B1[is_emb, ])
})
```

- [ ] **Step 2: Run to see them fail.** Run: install, then test `embedding-gradients`. Expected: both fail with "object 'gradient_agreement' not found" and "object 'embedding_effect' not found".

- [ ] **Step 3: Extract the two helpers** in `R/grouped-convergence.R`, each with a comment in the file's style, and call them from `grouped_effect_checks()` (for `agreement_pairs`) and `grouped_convergence()` (for `B_mean` and `E`). The code moves; nothing else changes. `grouped_effect_checks()` still computes `r`, `disagreement_pairs` and everything else itself.

- [ ] **Step 4: Run.** Run: install, then test `embedding-gradients|grouped-convergence|grouped-prior-fits|diagnostics`. Expected: all pass, with the same number of expectations in the three existing files as before the change (run them once on the untouched branch first and note the counts).

- [ ] **Step 5: Commit.**

```bash
git add R/grouped-convergence.R tests/testthat/test-embedding-gradients.R
git commit -m "Share the embedding effect and chain agreement of the convergence check"
```

### Task 2: `returnEmbeddingGradients()`

**Files:**
- Create: `R/embedding-gradients.R`, `tests/testthat/helper-embedding-gradients.R`
- Modify: `tests/testthat/test-embedding-gradients.R`, `tests/testthat/test-grouped-prior-fits.R` (append only), `NAMESPACE` and `man/` (by document)

**Interfaces:**
- Consumes: `embedding_effect()`, `gradient_agreement()`, `blockRhatESS()`, `is_grouped_fit()`, `GROUPED_LIMITS`.
- Produces: `returnEmbeddingGradients(fitModel, n_gradients = NULL)` returning the list the spec's Design section gives: `n_gradients` (integer), `rule` (`"hold 0.90"` or `"set"`), `site_scores` (n by r, row names the site names as text, columns `gradient1` on), `species_loadings` (S by r), `weights` (K by r), `share` (data frame `gradient`, `share`, `cumulative`, one row per direction up to the numerical rank), `size`, `agreement`, `agreement_pairs`, `by_chain` (a list with one entry per chain holding `site_scores`, `species_loadings`, `share`; `NULL` for a chain whose effect has rank below r).

- [ ] **Step 1: Write the test helper** `tests/testthat/helper-embedding-gradients.R`:

```r
# An n x K embedding block with centred, uncorrelated columns of standard
# deviation 1, and K x S coefficients that load on the columns `cols` only, so
# the effect X_emb %*% B has singular values `d` and its gradients are those
# columns.
low_rank_parts <- function(n = 40L, K = 6L, S = 12L, cols = 1:2, d = c(30, 15), seed = 1L) {
  set.seed(seed)
  Q <- qr.Q(qr(scale(matrix(rnorm(n * K), n, K), scale = FALSE)))
  V <- qr.Q(qr(matrix(rnorm(S * length(cols)), S)))
  B <- matrix(0, K, S); B[cols, ] <- diag(d, length(cols)) %*% t(V) / sqrt(n - 1)
  list(X_emb = Q * sqrt(n - 1), B = B)
}

# A fit object holding only what the gradient functions read. `B_emb` is a
# list with one K x S matrix per chain: that chain's posterior mean of the
# embedding coefficients. A measured column comes first.
fake_gradient_fit <- function(X_emb, B_emb, niter = 20L, seed = 1L,
                              site_names = paste0("s", seq_len(nrow(X_emb)))) {
  set.seed(seed)
  n <- nrow(X_emb); K <- ncol(X_emb); S <- ncol(B_emb[[1]]); nchain <- length(B_emb)
  colnames(X_emb) <- paste0("emb", seq_len(K))
  B <- array(0, c(K + 1L, S, niter, nchain))
  for (ch in seq_len(nchain)) for (j in seq_len(S)) for (i in seq_len(K + 1L)) {
    noise <- rnorm(niter, sd = 0.01)
    B[i, j, , ch] <- (if (i == 1L) 0.3 else B_emb[[ch]][i - 1L, j]) + noise - mean(noise)
  }
  list(X_psi = cbind(meas1 = as.numeric(scale(rnorm(n))), X_emb),
       infos = list(coef_prior = list(type = "grouped", is_emb = c(FALSE, rep(TRUE, K))),
                    siteNames = site_names, speciesNames = paste0("sp", seq_len(S))),
       results_output = list(jsdm_output = list(B_output = B)))
}
```

- [ ] **Step 2: Write the failing tests**, appended to `tests/testthat/test-embedding-gradients.R`:

```r
parts <- low_rank_parts()
fit2 <- fake_gradient_fit(parts$X_emb, list(parts$B, parts$B))

test_that("the gradients rebuild the pooled effect and are weighted sums of the embedding", {
  g <- returnEmbeddingGradients(fit2)
  expect_identical(g$n_gradients, 2L); expect_identical(g$rule, "hold 0.90")
  expect_equal(unname(g$site_scores %*% t(g$species_loadings)), parts$X_emb %*% parts$B, tolerance = 1e-6)
  expect_equal(unname(parts$X_emb %*% g$weights), unname(g$site_scores), tolerance = 1e-6)
  expect_equal(unname(colMeans(g$site_scores)), c(0, 0)); expect_equal(unname(colMeans(g$site_scores^2)), c(1, 1))
  expect_identical(rownames(g$site_scores), fit2$infos$siteNames)
  expect_identical(colnames(g$site_scores), c("gradient1", "gradient2"))
  expect_identical(rownames(g$species_loadings), fit2$infos$speciesNames)
  expect_identical(rownames(g$weights), paste0("emb", 1:6))
  expect_identical(names(g$share), c("gradient", "share", "cumulative"))
  expect_equal(g$share$share[1:2], c(900, 225) / 1125, tolerance = 1e-6)
  expect_equal(g$agreement, 1); expect_length(g$by_chain, 2L)
  for (k in 1:2) expect_gt(g$weights[which.max(abs(g$weights[, k])), k], 0)
})

test_that("n_gradients overrides the rule and is checked", {
  g1 <- returnEmbeddingGradients(fit2, n_gradients = 1L)
  expect_identical(g1$n_gradients, 1L); expect_identical(g1$rule, "set")
  expect_identical(dim(g1$site_scores), c(40L, 1L)); expect_identical(dim(g1$species_loadings), c(12L, 1L))
  expect_identical(dim(g1$weights), c(6L, 1L)); expect_identical(dim(g1$by_chain[[1]]$site_scores), c(40L, 1L))
  expect_error(returnEmbeddingGradients(fit2, n_gradients = 0), "single positive whole number")
  expect_error(returnEmbeddingGradients(fit2, n_gradients = 1.5), "single positive whole number")
  expect_error(returnEmbeddingGradients(fit2, n_gradients = 5L), "rank")
})

test_that("chains with different gradients are reported and warned about", {
  other <- low_rank_parts(cols = 3:4)
  fit <- fake_gradient_fit(parts$X_emb, list(parts$B, other$B))
  expect_warning(g <- returnEmbeddingGradients(fit, n_gradients = 2L), "different gradients")
  expect_lt(g$agreement, 0.90); expect_identical(dim(g$agreement_pairs), c(2L, 2L))
  one <- g$by_chain[[1]]
  expect_equal(unname(one$site_scores %*% t(one$species_loadings)), parts$X_emb %*% parts$B, tolerance = 1e-6)
  # Three gradients are more than either chain's own effect holds (each has rank 2).
  expect_warning(expect_warning(g3 <- returnEmbeddingGradients(fit, n_gradients = 3L), "different gradients"), "rank")
  expect_length(g3$by_chain, 2L); expect_null(g3$by_chain[[1]]); expect_null(g3$by_chain[[2]])
})

test_that("a small effect still returns gradients, with a warning", {
  small <- low_rank_parts(d = c(0.9, 0.5))
  fit <- fake_gradient_fit(small$X_emb, list(small$B, small$B))
  expect_warning(g <- returnEmbeddingGradients(fit), "small")
  expect_identical(g$n_gradients, 2L); expect_lt(g$size, 0.5)
})

test_that("one chain gives no agreement, and other fits are refused", {
  g <- returnEmbeddingGradients(fake_gradient_fit(parts$X_emb, list(parts$B)))
  expect_true(is.na(g$agreement)); expect_length(g$by_chain, 1L)
  bad <- fit2; bad$infos$coef_prior$type <- "fixed"
  expect_error(returnEmbeddingGradients(bad), "needs a fit with embeddingVars")
  nonfinite <- fit2; nonfinite$results_output$jsdm_output$B_output[2, 1, 1, 1] <- NA
  expect_error(returnEmbeddingGradients(nonfinite), "missing or non-finite")
})

test_that("the default matches the convergence check's own numbers", {
  other <- low_rank_parts(cols = 2:3)
  fit <- fake_gradient_fit(parts$X_emb, list(parts$B, other$B))
  B <- fit$results_output$jsdm_output$B_output; dB <- dim(B)
  cm <- occJSDM:::blockRhatESS(B, dB[1] * dB[2], dB[3], dB[4])$chain_mean
  chk <- occJSDM:::grouped_effect_checks(
    occJSDM:::embedding_effect(cm, fit$X_psi, fit$infos$coef_prior$is_emb, dB[2])$E)
  g <- suppressWarnings(returnEmbeddingGradients(fit))
  expect_true(chk$assessed); expect_identical(g$n_gradients, chk$n_gradients)
  expect_equal(g$agreement, chk$agreement); expect_equal(g$agreement_pairs, chk$agreement_pairs)
  expect_equal(g$size, chk$size)
})
```

And appended to the end of `tests/testthat/test-grouped-prior-fits.R`:

```r
test_that("returnEmbeddingGradients() gives one named score per site on real fits", {
  for (fit in list(grouped_fit(),
                   fit_grouped(simulate_grouped_fixture(model = "occupancy", p_emb = 6L, S = 8L)))) {
    g <- suppressWarnings(returnEmbeddingGradients(fit))
    expect_identical(nrow(g$site_scores), nrow(fit$X_psi))
    expect_identical(rownames(g$site_scores), as.character(fit$infos$siteNames))
    expect_identical(nrow(g$weights), sum(fit$infos$coef_prior$is_emb))
    conv <- suppressWarnings(returnGroupedConvergence(fit))
    if (conv$effect$assessed) {
      expect_identical(g$n_gradients, conv$effect$n_gradients); expect_equal(g$agreement, conv$effect$agreement)
    }
  }
})
```

- [ ] **Step 3: Run to see them fail.** Run: install, then test `embedding-gradients`. Expected: the new tests fail with "could not find function \"returnEmbeddingGradients\"".

- [ ] **Step 4: Implement `returnEmbeddingGradients(fitModel, n_gradients = NULL)`** in `R/embedding-gradients.R`, with roxygen (`@export`, `@md`, `@seealso` the `grouped-prior` page and `compareEmbeddingGradients()`), following the spec's "Gradients" mathematics and its Design section. Points the spec and tests fix:
    - Acceptance and message as `returnEmbeddingShare()`: stop with "returnEmbeddingGradients() needs a fit with embeddingVars".
    - Chain means from `blockRhatESS(B_output, p * S, niter, nchain)$chain_mean`, then `embedding_effect()`. Stop with a message containing "missing or non-finite" when the effect is not finite, and when it is all zero.
    - Numerical rank of the pooled effect: the count of singular values above `max(n, S) * d[1] * .Machine$double.eps`. `n_gradients` must be `NULL` or a single positive whole number (message "n_gradients must be a single positive whole number") no larger than that rank (message containing "rank" and its value).
    - Signs by the largest absolute weight, applied to scores, loadings and weights together; per chain by that chain's own weights.
    - Warnings: one containing "small" when `size` is below `GROUPED_LIMITS$small_effect`; one containing "different gradients", naming `r` and `returnGroupedConvergence()`, when `agreement` is below `GROUPED_LIMITS$gradients`; and a single one containing "rank", naming the chains, when any chain's own effect has rank below `r` and its entry in `by_chain` is therefore `NULL`. `by_chain` always has one entry per chain.
    - `drop = FALSE` wherever a matrix is subset, and `diag(x, nrow = r)` wherever a diagonal matrix is built, so one gradient works.

- [ ] **Step 5: Run.** Run: install, document, install, then test `embedding-gradients|grouped-prior-fits|grouped-convergence`. Expected: all pass; `NAMESPACE` gains `export(returnEmbeddingGradients)` and `man/returnEmbeddingGradients.Rd` exists.

- [ ] **Step 6: Commit.**

```bash
git add R/embedding-gradients.R NAMESPACE man/returnEmbeddingGradients.Rd tests/testthat/helper-embedding-gradients.R tests/testthat/test-embedding-gradients.R tests/testthat/test-grouped-prior-fits.R
git commit -m "Add returnEmbeddingGradients(): gradients of the whole embedding response"
```

### Task 3: `compareEmbeddingGradients()`: sites, covariates and check 1

**Files:**
- Modify: `R/embedding-gradients.R`, `NAMESPACE` and `man/` (by document)
- Create: `tests/testthat/test-compare-embedding-gradients.R`

**Interfaces:**
- Consumes: the list from `returnEmbeddingGradients()`. Only `n_gradients`, `site_scores`, `species_loadings` and `by_chain` are read, so a hand-built list with those four works; `by_chain` may be an empty list.
- Produces: `compareEmbeddingGradients(gradients, covariates, site_col = "Site", squared = FALSE)` returning a list with `covariate_from_gradients`, `gradients_from_covariates`, `unmeasured` and `sites`. This task builds the first and last; the middle two are `NULL` until Task 4. `sites` is `list(n_used, not_in_fit, not_in_covariates, missing_values)`, the last three being character vectors of site names. Sites used keep the fit's order.

- [ ] **Step 1: Write the failing tests** in `tests/testthat/test-compare-embedding-gradients.R`:

```r
# compareEmbeddingGradients() (R/embedding-gradients.R; design in
# dev/superpowers/specs/2026-11-20-embedding-gradients-design.md).

parts <- low_rank_parts()
g <- returnEmbeddingGradients(fake_gradient_fit(parts$X_emb, list(parts$B, parts$B)))
Z <- g$site_scores
set.seed(3)
covs <- data.frame(Site = rownames(Z), inplane = as.numeric(Z %*% c(2, -1)) + 5, noise = rnorm(nrow(Z)))

test_that("check 1 finds a covariate that lies in the gradients, and not one that does not", {
  out <- compareEmbeddingGradients(g, covs)
  expect_identical(names(out), c("covariate_from_gradients", "gradients_from_covariates", "unmeasured", "sites"))
  c1 <- out$covariate_from_gradients
  expect_identical(names(c1), c("covariate", "chain", "n_gradients", "r_squared", "adj_r_squared"))
  expect_setequal(c1$chain, c("pooled", "1", "2")); expect_setequal(c1$n_gradients, 1:2)
  pick <- function(cv, k) c1[c1$covariate == cv & c1$chain == "pooled" & c1$n_gradients == k, ]
  expect_equal(pick("inplane", 2)$r_squared, 1); expect_lt(pick("inplane", 1)$r_squared, 1)
  expect_lt(abs(pick("noise", 2)$adj_r_squared), 0.2)
  expect_identical(out$sites$n_used, 40L)
})

test_that("sites are matched by name, whatever the row order or the type of the names", {
  ref <- compareEmbeddingGradients(g, covs)
  set.seed(9); expect_equal(compareEmbeddingGradients(g, covs[sample(nrow(covs)), ]), ref)
  f <- covs; f$Site <- factor(f$Site); expect_equal(compareEmbeddingGradients(g, f), ref)
  named <- covs; names(named)[1] <- "plot"
  expect_equal(compareEmbeddingGradients(g, named, site_col = "plot"), ref)
  gn <- returnEmbeddingGradients(fake_gradient_fit(parts$X_emb, list(parts$B, parts$B), site_names = 1:40))
  cn <- covs; cn$Site <- 1:40
  expect_equal(compareEmbeddingGradients(gn, cn)$covariate_from_gradients, ref$covariate_from_gradients)
})

test_that("one chain gives pooled and chain 1 only", {
  g1 <- returnEmbeddingGradients(fake_gradient_fit(parts$X_emb, list(parts$B)))
  expect_setequal(compareEmbeddingGradients(g1, covs)$covariate_from_gradients$chain, c("pooled", "1"))
})

test_that("rows and sites that cannot be used are left out and counted", {
  extra <- rbind(covs, data.frame(Site = "elsewhere", inplane = 1, noise = 0))
  expect_message(o1 <- compareEmbeddingGradients(g, extra), "1 row")
  expect_identical(o1$sites$not_in_fit, "elsewhere"); expect_identical(o1$sites$n_used, 40L)
  partial <- covs[-(1:3), ]; partial$noise[1] <- NA
  expect_message(o2 <- compareEmbeddingGradients(g, partial), "3 fitted site")
  expect_identical(o2$sites$n_used, 36L)
  expect_identical(o2$sites$not_in_covariates, rownames(Z)[1:3])
  expect_identical(o2$sites$missing_values, rownames(Z)[4])
})

test_that("unusable input stops", {
  expect_error(compareEmbeddingGradients(g, covs, site_col = "nope"), "site_col")
  expect_error(compareEmbeddingGradients(g, rbind(covs, covs[1, ])), "more than once")
  const <- covs; const$flat <- 1; expect_error(compareEmbeddingGradients(g, const), "flat")
  expect_error(suppressMessages(compareEmbeddingGradients(g, covs[1:4, ])), "too few sites")
  stranger <- covs; stranger$Site <- paste0("x", 1:40)
  expect_error(suppressMessages(compareEmbeddingGradients(g, stranger)), "No site")
  expect_error(compareEmbeddingGradients(list(1), covs), "returnEmbeddingGradients")
})

test_that("categorical covariates are assessed level by level in check 1", {
  cat3 <- covs; cat3$habitat <- rep(c("bog", "fen", "pond", "fen"), 10); cat3$two <- rep(c("a", "b"), 20)
  c1 <- compareEmbeddingGradients(g, cat3)$covariate_from_gradients
  expect_setequal(unique(c1$covariate),
                  c("inplane", "noise", "habitat: bog", "habitat: fen", "habitat: pond", "two"))
})
```

- [ ] **Step 2: Run to see them fail.** Run: install, then test `compare-embedding`. Expected: fail with "could not find function \"compareEmbeddingGradients\"".

- [ ] **Step 3: Implement the function's inputs and check 1**, with roxygen (`@export`, `@md`), following the spec. Points the spec and tests fix:
    - `gradients` must be a list with `n_gradients`, `site_scores` and `species_loadings`; otherwise stop naming `returnEmbeddingGradients()`. Only the row names of `site_scores` are relied on; gradient names in the output are made here (`gradient1` on), never read from column names. `sites$n_used` is an integer.
    - Read the site column with `covariates[[site_col]]` and compare `as.character()` values with `rownames(gradients$site_scores)`. Stop when the column is absent (message containing "site_col"), when a site name is repeated (containing "more than once"), when no site matches (starting "No site"), or when `m` is under `q + 3` or under `r + 2` (containing "too few sites").
    - Messages, each only when its count is over 0: "N row(s) of covariates are for sites not in the fit", "N fitted site(s) are not in covariates", "N site(s) have missing covariate values".
    - Logical columns become 0 and 1; character columns become factors; a covariate constant over the sites used stops with its name in the message.
    - The design matrix comes from `stats::model.matrix()` on the covariates; every projection is by `qr.resid()` or `qr.fitted()`, never by inverting `F'F`; `q` is the QR rank minus 1.
    - Check 1 rows: each numeric or logical covariate under its own name; a two-level categorical under its own name, as the indicator of its second level; a categorical with more levels as one row set per level named `"<covariate>: <level>"`. For every chain entry of `by_chain` that is not `NULL`, labelled by its position, and for `"pooled"`. Adjusted R-squared uses `m - k - 1`.

- [ ] **Step 4: Run.** Run: install, document, install, then test `compare-embedding|embedding-gradients`. Expected: all pass.

- [ ] **Step 5: Commit.**

```bash
git add R/embedding-gradients.R NAMESPACE man/compareEmbeddingGradients.Rd tests/testthat/test-compare-embedding-gradients.R
git commit -m "Add compareEmbeddingGradients(): match sites and predict covariates from gradients"
```

### Task 4: Checks 2 and 3: what the covariates explain, and the unmeasured gradients

**Files:**
- Modify: `R/embedding-gradients.R`, `tests/testthat/test-compare-embedding-gradients.R`, `man/compareEmbeddingGradients.Rd` (by document)

**Interfaces:**
- Consumes: Task 3's site matching and design matrix.
- Produces: `gradients_from_covariates`, a list of `overall` (data frame `chain`, `explained`, `adj_explained`), `by_covariate` (data frame `covariate`, `alone`, `unique`; pooled; one row per column of `covariates`, in their order) and `by_gradient` (data frame `gradient`, `share`, `r_squared`, `adj_r_squared`; pooled). And `unmeasured`, a list of `site_scores` (sites used by components kept, columns `unmeasured1` on), `species_loadings` and `share` (data frame `gradient`, `share`, `cumulative`).

- [ ] **Step 1: Write the failing tests**, appended to `tests/testthat/test-compare-embedding-gradients.R`:

```r
pooled <- function(out) out$gradients_from_covariates$overall$explained[
  out$gradients_from_covariates$overall$chain == "pooled"]
set.seed(4)
mix <- data.frame(Site = rownames(Z), a = Z[, 1] + rnorm(40, sd = 0.5), b = Z[, 1] + Z[, 2] + rnorm(40, sd = 0.5))
mix$dup <- mix$a

test_that("covariates that span the gradients explain everything and leave nothing unmeasured", {
  out <- compareEmbeddingGradients(g, data.frame(Site = rownames(Z), a = Z[, 1] + Z[, 2], b = Z[, 1] - 2 * Z[, 2]))
  expect_identical(names(out$gradients_from_covariates), c("overall", "by_covariate", "by_gradient"))
  expect_setequal(out$gradients_from_covariates$overall$chain, c("pooled", "1", "2"))
  expect_equal(pooled(out), 1)
  expect_identical(ncol(out$unmeasured$site_scores), 0L); expect_identical(ncol(out$unmeasured$species_loadings), 0L)
  expect_identical(nrow(out$unmeasured$share), 0L)
})

test_that("a covariate that is one gradient leaves the other as the unmeasured gradient", {
  out <- compareEmbeddingGradients(g, data.frame(Site = rownames(Z), a = Z[, 1]))
  bg <- out$gradients_from_covariates$by_gradient; u <- out$unmeasured
  expect_equal(pooled(out), 0.8, tolerance = 1e-6)
  expect_equal(bg$r_squared, c(1, 0)); expect_equal(bg$share, c(0.8, 0.2), tolerance = 1e-6)
  expect_identical(dim(u$site_scores), c(40L, 1L)); expect_identical(colnames(u$site_scores), "unmeasured1")
  expect_identical(rownames(u$site_scores), rownames(Z))
  expect_equal(abs(cor(u$site_scores[, 1], Z[, 2])), 1); expect_equal(u$share$share, 0.2, tolerance = 1e-6)
  expect_equal(unname(u$site_scores %*% t(u$species_loadings)),
               unname(Z[, 2, drop = FALSE] %*% t(g$species_loadings[, 2, drop = FALSE])), tolerance = 1e-6)
})

test_that("alone, unique and the unmeasured shares add up, with a duplicated covariate", {
  out <- compareEmbeddingGradients(g, mix); c2 <- out$gradients_from_covariates
  expect_equal(sum(out$unmeasured$share$share), 1 - pooled(out))
  expect_identical(c2$by_covariate$covariate, c("a", "b", "dup"))
  expect_true(all(c2$by_covariate$alone <= pooled(out) + 1e-12))
  expect_equal(c2$by_covariate$unique[c(1, 3)], c(0, 0)); expect_gt(c2$by_covariate$unique[2], 0)
  expect_equal(c2$by_covariate$alone[1], c2$by_covariate$alone[3])
  expect_equal(sum(c2$by_gradient$share * c2$by_gradient$r_squared), pooled(out))
  expect_lt(max(abs(cor(out$unmeasured$site_scores, cbind(mix$a, mix$b)))), 1e-8)
})

test_that("findings do not depend on how the gradients are rotated", {
  th <- 0.7; Q <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2)
  plain <- g; plain$by_chain <- list()
  rot <- plain; rot$site_scores <- g$site_scores %*% Q; rot$species_loadings <- g$species_loadings %*% Q
  rownames(rot$site_scores) <- rownames(Z)
  a <- compareEmbeddingGradients(plain, mix); b <- compareEmbeddingGradients(rot, mix)
  expect_equal(b$gradients_from_covariates$overall, a$gradients_from_covariates$overall)
  expect_equal(b$gradients_from_covariates$by_covariate, a$gradients_from_covariates$by_covariate)
  expect_equal(b$unmeasured$share, a$unmeasured$share)
  top <- function(x) { t <- x$covariate_from_gradients; t[t$n_gradients == 2, ] }
  expect_equal(top(b), top(a))
  expect_equal(abs(diag(cor(b$unmeasured$site_scores, a$unmeasured$site_scores))),
               rep(1, ncol(a$unmeasured$site_scores)))
})

test_that("one gradient works throughout", {
  g1 <- g; g1$n_gradients <- 1L; g1$by_chain <- list()
  g1$site_scores <- g$site_scores[, 1, drop = FALSE]; g1$species_loadings <- g$species_loadings[, 1, drop = FALSE]
  out <- compareEmbeddingGradients(g1, mix)
  expect_identical(nrow(out$gradients_from_covariates$by_gradient), 1L)
  expect_equal(out$gradients_from_covariates$by_gradient$r_squared, pooled(out))
  expect_lte(ncol(out$unmeasured$site_scores), 1L)
})

test_that("squared terms catch a humped response", {
  set.seed(5); x <- rnorm(40); rms1 <- function(v) v / sqrt(mean(v^2))
  hand <- list(n_gradients = 2L, by_chain = list(),
               site_scores = cbind(rms1(x - mean(x)), rms1(resid(lm(I(x^2) ~ x)))),
               species_loadings = qr.Q(qr(matrix(rnorm(24), 12))) %*% diag(c(3, 2)))
  rownames(hand$site_scores) <- paste0("s", 1:40)
  hump <- data.frame(Site = paste0("s", 1:40), x = x)
  expect_equal(pooled(compareEmbeddingGradients(hand, hump)), 9 / 13)
  sq <- compareEmbeddingGradients(hand, hump, squared = TRUE)
  expect_equal(pooled(sq), 1)
  expect_identical(sq$gradients_from_covariates$by_covariate$covariate, "x")
  expect_identical(unique(sq$covariate_from_gradients$covariate), "x")
})
```

- [ ] **Step 2: Run to see them fail.** Run: install, then test `compare-embedding`. Expected: the six new tests fail (the two elements are `NULL`); Task 3's still pass.

- [ ] **Step 3: Implement checks 2 and 3**, following the spec's mathematics with `G = Z_I %*% (L'L)^{1/2}` (the symmetric square root, from `eigen()` of `crossprod(L)`), so that a rotated gradients list gives the same answers. Points the spec and tests fix:
    - `explained` is 1 minus the residual sum of squares of `G` on `F` over the sum of squares of centred `G`; `adj_explained` uses `m - q - 1`.
    - `alone` uses the intercept and that covariate's columns of the design matrix (its indicator columns, and its square when `squared = TRUE`); `unique` is `explained` minus `explained` without those columns, set to 0 when it is below 0 by rounding.
    - `by_gradient` uses the columns of `Z_I` scaled by `sqrt(colSums(L^2))`, and `share` is each column's share of the centred sum of squares.
    - With `squared = TRUE`, add the square of the centred values of every numeric covariate that has more than two distinct values. Squares enter checks 2 and 3 only.
    - Unmeasured components are kept when their share is above `1e-10`. Site scores are `sqrt(m)` times the left singular vectors; loadings follow the spec's `L_u`; each component is signed so that its species loading of largest absolute value is positive. With none kept, the two matrices have zero columns and `share` zero rows.
    - Per chain, `overall` is computed from each non-`NULL` entry of `by_chain`; `by_covariate`, `by_gradient` and `unmeasured` are pooled only.

- [ ] **Step 4: Run.** Run: install, document, install, then test `compare-embedding|embedding-gradients`. Expected: all pass.

- [ ] **Step 5: Commit.**

```bash
git add R/embedding-gradients.R man/compareEmbeddingGradients.Rd tests/testthat/test-compare-embedding-gradients.R
git commit -m "compareEmbeddingGradients(): what covariates explain, and unmeasured gradients"
```

### Task 5: Help page, site index and changelog

**Files:**
- Modify: `R/coef-prior.R` (the `grouped-prior` roxygen block, lines 444 to 670), `R/output.R` (the roxygen of `returnEmbeddingBlends()`, line 335), `_pkgdown.yml` (line 70), `NEWS.md`, `man/` (by document)

**Interfaces:**
- Consumes: the two exported functions and their help pages.

- [ ] **Step 1: Add the help section.** In `R/coef-prior.R`, after the section "Reading a management covariate", add `@section Interpreting the gradients:` with the text of the spec's appendix, as roxygen markdown, keeping every caution and the recipe for passing `returnOrdinationScores()` output. Check that recipe against `returnOrdinationScores()`'s actual return value before writing it, and write the two or three lines of R that do it.

- [ ] **Step 2: Update the pointers.** In the same block, change "Five functions summarise these" to "Seven" and add the two new functions to the list, one line each. In the roxygen of `returnEmbeddingBlends()`, after the paragraph on rotation, add: "This function works from the shared part of the responses only. \[returnEmbeddingGradients()\] works from the whole response, and is the one to use for interpretation." (written as a roxygen link). In `_pkgdown.yml` add `returnEmbeddingGradients` and `compareEmbeddingGradients` after `returnEmbeddingBlends`.

- [ ] **Step 3: Add the changelog entry.** Re-read `NEWS.md`, then add under "New features" of the development version one bullet in the style of its neighbours naming both functions and what each returns, and under "Documentation" one bullet for the new help section. Re-read the file afterwards.

- [ ] **Step 4: Regenerate and check.** Run: install, document, install. Then:

```bash
grep -c "Interpreting the gradients" man/grouped-prior.Rd
grep -c "returnEmbeddingGradients" man/returnEmbeddingBlends.Rd NAMESPACE _pkgdown.yml
git status --short
```

Expected: 1 for the first; at least 1 for each file in the second; and `git status` showing only the files this task names plus the three `.Rd` files.

- [ ] **Step 5: Run the whole suite.** Run: test with the filter removed (drop `filter = "FILTER", ` from the command), output redirected to `$SP/full-test.log`, and read its tail. Expected: no failures and no errors. Skips and warnings that were there before this branch are not this plan's.

- [ ] **Step 6: Commit.**

```bash
git add R/coef-prior.R R/output.R _pkgdown.yml NEWS.md man/grouped-prior.Rd man/returnEmbeddingBlends.Rd man/returnEmbeddingGradients.Rd man/compareEmbeddingGradients.Rd
git commit -m "Document how to interpret the embedding gradients"
```
