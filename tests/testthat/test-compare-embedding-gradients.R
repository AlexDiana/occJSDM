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
