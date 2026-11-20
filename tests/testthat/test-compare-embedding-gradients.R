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
  expect_message(expect_message(o2 <- compareEmbeddingGradients(g, partial), "3 fitted site"),
                 "1 site has missing covariate values")
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
  expect_equal(unname(abs(diag(cor(b$unmeasured$site_scores, a$unmeasured$site_scores)))),
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

test_that("a covariate with a large offset is still used", {
  out <- compareEmbeddingGradients(g, data.frame(Site = rownames(Z), when = 1.7e9 + Z[, 1]))
  expect_equal(pooled(out), 0.8, tolerance = 1e-6)
  expect_equal(out$gradients_from_covariates$by_covariate$alone, 0.8, tolerance = 1e-6)
  expect_equal(out$gradients_from_covariates$by_covariate$unique, 0.8, tolerance = 1e-6)
})

test_that("a label column is named when it leaves too few sites", {
  labelled <- covs; labelled$notes <- paste("plot", 1:40)
  expect_error(compareEmbeddingGradients(g, labelled), "too few sites.*notes")
})
