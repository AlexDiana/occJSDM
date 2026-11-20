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
  # The second chain loads on other embedding columns and, with the species
  # reversed, on other species, so the pooled effect has four directions.
  other <- low_rank_parts(cols = 3:4)
  fit <- fake_gradient_fit(parts$X_emb, list(parts$B, other$B[, 12:1]))
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

test_that("a fit without usable site names is refused, not numbered", {
  nameless <- fit2; nameless$infos$siteNames <- NULL
  expect_error(returnEmbeddingGradients(nameless), "siteNames")
  short <- fit2; short$infos$siteNames <- fit2$infos$siteNames[1:5]
  expect_error(returnEmbeddingGradients(short), "siteNames")
})

test_that("an embedding with no effect at all is refused", {
  zero <- fit2; zero$results_output$jsdm_output$B_output[-1, , , ] <- 0
  expect_error(returnEmbeddingGradients(zero), "zero at every fitted site")
})
