# Run from the repository root with the study library first in R_LIBS:
# Rscript dev/simstudy/waic-cross-validation/test-cv-lib.R
# Uses short fits of the real Lesson 4 survey; about a minute.
library(testthat)
source("dev/simstudy/waic-cross-validation/cv-lib.R")
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
data <- lesson$input$sim$data_list
sites <- sort(unique(data$info$Site))
short <- list(nchain = 2L, nburn = 20L, niter = 30L, nthin = 1L)

test_that("jobs, seeds and folds are fixed and complete", {
  jobs <- cv_job_names()
  expect_length(jobs, 22L)
  expect_false(anyDuplicated(jobs) > 0)
  parsed <- lapply(jobs, cv_parse_job)
  seeds <- vapply(parsed, `[[`, integer(1), "seed")
  expect_false(anyDuplicated(seeds) > 0)
  expect_identical(cv_parse_job("full-two")$seed, 20260921L)
  expect_identical(cv_parse_job("full-one")$seed, 20260926L)
  expect_identical(cv_parse_job("fold07-one")$fold, 7L)
  expect_error(cv_parse_job("fold7-one"))
  folds <- cv_make_folds(sites)
  expect_identical(folds, cv_make_folds(sites))
  expect_setequal(folds$Site, sites)
  expect_true(all(table(folds$fold) == 10L))
})

folds <- cv_make_folds(sites)
heldout_sites <- folds$Site[folds$fold == 1L]
train <- cv_subset(data, setdiff(sites, heldout_sites))
heldout <- cv_subset(data, heldout_sites)
fit <- suppressWarnings(cv_fit(train, "two", seed = 1L, mcmc = short))$fit

test_that("subsetting keeps every PCR row of the chosen sites", {
  expect_identical(nrow(train$info) + nrow(heldout$info), nrow(data$info))
  expect_identical(nrow(heldout$info), 10L * 3L * 2L * 6L)
  expect_setequal(unique(fit$infos$siteNames), setdiff(sites, heldout_sites))
})

test_that("the scoring copy reproduces the fitted design for training sites", {
  copy <- cv_heldout_fit(fit, train)
  expect_identical(unname(copy$X_psi), unname(fit$X_psi))
  expect_identical(colnames(copy$X_psi), colnames(fit$X_psi))
  expect_identical(unname(copy$X_theta), unname(fit$X_theta))
  expect_identical(colnames(copy$X_theta), colnames(fit$X_theta))
  for (k in c("idx_z_w", "idx_w_k", "idx_p_k")) {
    expect_identical(as.integer(copy$infos$list_idx[[k]]), as.integer(fit$infos$list_idx[[k]]))
  }
  expect_identical(unname(as.matrix(copy$infos$OTU)), unname(as.matrix(fit$infos$OTU)))
  expect_identical(as.character(copy$infos$siteNames), as.character(fit$infos$siteNames))
  own <- suppressWarnings(occJSDM::computeSiteWAIC(fit))
  via_copy <- cv_score(fit, train)
  expect_identical(via_copy$log_lik, own$log_lik)
})

copy <- cv_heldout_fit(fit, heldout)

test_that("withheld covariates use the training constants", {
  first <- heldout$info[!duplicated(heldout$info$Site), ]
  first <- first[order(first$Site), ]
  m <- fit$infos$list_X_psi_mat
  manual <- sweep(sweep(as.matrix(first[, m$names_df]), 2, m$mean_df, "-"), 2, m$sd_df, "/")
  expect_equal(unname(copy$X_psi), unname(manual), tolerance = 1e-14)
  samples <- heldout$info[!duplicated(heldout$info[c("Site", "Sample")]), ]
  samples <- samples[order(samples$Site, samples$Sample), ]
  t <- fit$infos$list_X_theta_mat
  expect_equal(unname(copy$X_theta[, 1]), rep(1, 30))
  expect_equal(unname(copy$X_theta[, 2]), (samples$X_theta - t$mean_df) / t$sd_df,
               tolerance = 1e-14)
  expect_false(isTRUE(all.equal(m$mean_df, colMeans(as.matrix(first[, m$names_df])))))
})

test_that("a site's score does not depend on which other sites are scored", {
  all_ten <- cv_score(fit, heldout)
  alone <- cv_score(fit, cv_subset(data, heldout_sites[3]))
  site <- as.character(heldout_sites[3])
  expect_lt(max(abs(alone$log_lik[, site] - all_ten$log_lik[, site])), 1e-3)
})

# Independent reference: the observed-data likelihood of one withheld site at
# one posterior draw, written directly from the formula in site-waic-plan.md
# from the raw PCR rows, integrating the two factors on a fine product grid.
brute_force <- function(fit, data, site, it, ch) {
  o <- fit$results_output; j <- o$jsdm_output
  rows <- data$info$Site == site
  info <- data$info[rows, ]; y <- 1 * (data$OTU[rows, , drop = FALSE] >= 1)
  S <- ncol(y); d <- fit$infos$n_factors
  m <- fit$infos$list_X_psi_mat
  x <- (unlist(info[1, m$names_df]) - m$mean_df) / m$sd_df
  eta <- j$B0_output[, it, ch] + drop(x %*% matrix(j$B_output[, , it, ch], length(x), S))
  t <- fit$infos$list_X_theta_mat
  primer <- match(info$Primer, fit$infos$primerNames)
  p <- matrix(o$p_output[, , it, ch], ncol = S); q <- matrix(o$q_output[, , it, ch], ncol = S)
  bt <- matrix(o$beta_theta_output[, , it, ch], ncol = S); theta0 <- o$theta0_output[, it, ch]
  g0 <- g1 <- numeric(S)
  for (s in seq_len(S)) for (sample in unique(info$Sample)) {
    r <- info$Sample == sample
    a <- prod(ifelse(y[r, s] == 1, p[primer[r], s], 1 - p[primer[r], s]))
    b <- prod(ifelse(y[r, s] == 1, q[primer[r], s], 1 - q[primer[r], s]))
    theta <- stats::plogis(bt[1, s] + bt[2, s] * (info$X_theta[r][1] - t$mean_df) / t$sd_df)
    g1[s] <- g1[s] + log(theta * a + (1 - theta) * b)
    g0[s] <- g0[s] + log(theta0[s] * a + (1 - theta0[s]) * b)
  }
  L <- matrix(j$L_output[, , it, ch], d, S) * j$sigmah_output[it, ch]
  step <- 0.02; u <- seq(-8, 8, by = step)
  grid <- as.matrix(expand.grid(rep(list(u), d)))
  log_w <- rowSums(stats::dnorm(grid, log = TRUE)) + d * log(step)
  shifted <- sweep(grid %*% L, 2, eta, "+")
  conditional <- rowSums(log(stats::plogis(-shifted) * rep(exp(g0), each = nrow(grid)) +
                               stats::plogis(shifted) * rep(exp(g1), each = nrow(grid))))
  total <- conditional + log_w
  max(total) + log(sum(exp(total - max(total))))
}

test_that("withheld scores match a brute-force likelihood from the raw data", {
  score <- cv_score(fit, heldout)
  iterations <- dim(fit$results_output$theta0_output)[2]
  for (check in list(c(it = 5L, ch = 1L), c(it = 30L, ch = 2L))) {
    draw <- check[["it"]] + (check[["ch"]] - 1L) * iterations
    for (site in heldout_sites[c(1, 6, 10)]) {
      expect_equal(unname(score$log_lik[draw, as.character(site)]),
                   brute_force(fit, data, site, check[["it"]], check[["ch"]]), tolerance = 1e-6)
    }
  }
})

test_that("fold summaries agree with the scorer and with each other", {
  score <- cv_score(fit, heldout)
  iterations <- dim(fit$results_output$theta0_output)[2]
  summary <- cv_fold_summary(score, iterations, 2L)
  expect_equal(summary$sites$lppd, unname(score$pointwise$lppd), tolerance = 1e-12)
  expect_equal(summary$sites$chain1,
               unname(apply(score$log_lik[1:iterations, ], 2, cv_log_mean_exp)), tolerance = 1e-12)
  expect_true(is.finite(summary$sum_mcse) && summary$sum_mcse > 0)
})

test_that("the screen covers every scored quantity and fitted cell", {
  screen <- cv_screen(fit)
  counts <- table(screen$param)
  expect_equal(as.integer(counts[["psi"]]), 90L * 10L)
  expect_equal(as.integer(counts[["beta_psi"]]), 2L * 10L)
  expect_equal(as.integer(counts[["p"]]), 2L * 10L)
  expect_type(screen$flagged, "logical")
})

test_that("one-factor fits score through the same route", {
  one <- suppressWarnings(cv_fit(train, "one", seed = 2L, mcmc = short))$fit
  expect_identical(cv_score(one, train)$log_lik,
                   suppressWarnings(occJSDM::computeSiteWAIC(one))$log_lik)
  score <- cv_score(one, heldout)
  expect_equal(unname(score$log_lik[3, as.character(heldout_sites[2])]),
               brute_force(one, data, heldout_sites[2], 3L, 1L), tolerance = 1e-6)
})
