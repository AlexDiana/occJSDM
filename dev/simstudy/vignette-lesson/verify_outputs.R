# Independent numerical checks for the Lesson 3 summaries; no MCMC is run.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
x <- readRDS("vignettes/teaching-data/output-lesson.rds")
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
b <- lesson$input
t <- b$sim$true_params$jsdmParams_true
sp <- colnames(b$sim$data_list$OTU)
env <- grep("^X_psi", names(b$sim$data_list$info), value = TRUE)
traits <- colnames(b$sim$data_list$traits)
close <- function(a, b, tolerance = 1e-10) {
  stopifnot(identical(dim(a), dim(b)), length(a) == length(b),
            all(is.finite(a)), all(is.finite(b)), max(abs(a - b)) < tolerance)
}
stopifnot(identical(x$source_hashes, lesson_source_hashes()),
          identical(x$lesson_md5, unname(tools::md5sum("vignettes/teaching-data/nonspatial-lesson.rds"))),
          identical(x$legacy_fit_md5, unname(tools::md5sum("data/sampleresults.rda"))),
          identical(x$exporter_md5, unname(tools::md5sum("dev/simstudy/vignette-lesson/summarise_outputs.R"))))

for (arm in c("perfect", "default")) {
  path <- file.path(archive, x$fit_manifests[[arm]]$file)
  stopifnot(identical(unname(tools::md5sum(path)), x$fit_manifests[[arm]]$md5))
  r <- readRDS(path)
  f <- r$fit
  validate_lesson_fit_identity(f, b)
  j <- f$results_output$jsdm_output
  stopifnot(identical(r$source_hashes, b$source_hashes), identical(r$input_md5, lesson$input_md5))
  close(unname(f$Tr), unname(scale(b$sim$data_list$traits)))
  close(unname(sweep(f$X_psi %*% t$B + t$U %*% t$L, 2, t$B0, "+")),
        unname(t$eta))

  coefficients <- x$coefficients[x$coefficients$arm == arm, ]
  stopifnot(nrow(coefficients) == 24L)
  for (i in seq_len(nrow(coefficients))) {
    a <- coefficients[i, ]
    e <- match(a$covariate, env)
    if (a$block == "Environment") {
      s <- match(a$term, sp)
      draws <- j$B_output[e, s, , ]
      expected_truth <- t$B[e, s]
    } else {
      g <- match(a$term, traits)
      draws <- j$G_output[g, e, , ]
      expected_truth <- t$G[g, e] * sd(b$sim$data_list$traits[, g])
    }
    close(as.numeric(a[c("estimate", "lower", "upper", "truth")]),
          c(mean(draws), quantile(draws, c(.025, .975)), expected_truth))
    close(a$rhat, posterior::rhat(draws))
    close(a$ess_bulk, posterior::ess_bulk(draws))
    close(a$ess_tail, posterior::ess_tail(draws))
    stopifnot(a$excludes_zero == (min(draws |> quantile(c(.025, .975))) > 0 ||
                                   max(draws |> quantile(c(.025, .975))) < 0))
  }

  # Check all truth curves and selected fitted endpoints by directly summing
  # components, independently of returnOccupancyGradient().
  curves <- x$gradients[x$gradients$arm == arm, ]
  stopifnot(nrow(curves) == 800L, all(is.finite(as.matrix(curves[c("truth", "med", "low", "high")]))))
  for (i in seq_len(nrow(curves))) {
    a <- curves[i, ]
    s <- match(a$species, sp)
    e <- match(a$covariate, env)
    values <- apply(f$X_psi, 2, median)
    values[e] <- a$x
    close(a$truth, plogis(t$B0[s] + sum(values * t$B[, s])))
    if (i %% 40L %in% c(0L, 1L)) {
      eta <- j$B0_output[s, , ]
      for (k in seq_along(env)) eta <- eta + values[k] * j$B_output[k, s, , ]
      close(as.numeric(a[c("low", "med", "high")]),
            unname(quantile(plogis(eta), c(.025, .5, .975))))
    }
  }

  # Independently reconstruct combined factor contributions for three sites
  # and every species. This also detects a site/species table-order swap.
  for (s in seq_along(sp)) for (i in c(1L, 19L, 100L)) {
    a <- x$residual[x$residual$arm == arm & x$residual$species == sp[s] &
                      x$residual$Site == as.character(f$infos$siteNames[i]), ]
    draws <- j$U_output[i, 1, , ] * j$L_output[1, s, , ] +
      j$U_output[i, 2, , ] * j$L_output[2, s, , ]
    close(a$estimate, mean(draws))
    close(a$truth, sum(t$U[i, ] * t$L[, s]))
  }

  # Correlation intervals from the raw loading arrays, not the public getter.
  for (s1 in c(1L, 4L, 8L)) for (s2 in c(2L, 5L, 10L)) {
    l1 <- matrix(j$L_output[, s1, , ], nrow = 2)
    l2 <- matrix(j$L_output[, s2, , ], nrow = 2)
    draws <- colSums(l1 * l2) / sqrt(colSums(l1^2) * colSums(l2^2))
    a <- x$correlations[x$correlations$arm == arm & x$correlations$species1 == sp[s1] &
                         x$correlations$species2 == sp[s2], ]
    close(as.numeric(a[c("lower", "estimate", "upper")]),
          unname(quantile(draws, c(.025, .5, .975))))
    if (s1 == 4L) stopifnot(is.na(a$truth))
  }

  # In this non-spatial model the package partition reduces to two component
  # shares. Recompute it using ordinary sd(), independently of its helper.
  for (s in seq_along(sp)) {
    E <- t$B0[s] + drop(f$X_psi %*% t$B[, s])
    F <- drop(t$U %*% t$L[, s])
    ve <- sd(plogis(E)); vf <- sd(plogis(F)); vef <- sd(plogis(E + F))
    ce <- ve + max(vef - vf, 0)
    cf <- vf + max(vef - ve, 0)
    a <- x$variation[x$variation$arm == arm & x$variation$species == sp[s], ]
    close(a$truth[a$component == "Environmental"], ce / (ce + cf))
    close(sum(a$truth), 1)
    close(sum(a$estimate), 1)
  }
  cat(arm, ": identities, coefficients, profiles, residual products, correlations and partition truth verified\n")
  rm(f, r, j)
  invisible(gc())
}

# Verify the trait-cancellation algebra from a different regression calculation.
X <- cbind(1, scale(b$sim$data_list$traits))
net <- qr.solve(X, t(t$B))[-1, , drop = FALSE]
close(unname(net), unname(x$realized_trait_slopes))
for (g in seq_along(traits)) for (e in seq_along(env)) {
  a <- x$trait_components[x$trait_components$trait == traits[g] &
                           x$trait_components$covariate == env[e], ]
  close(sum(a$value), net[g, e])
}

# Independent expectation: sum over the possible number of collected samples,
# using a Binomial(M, theta) probability for each count.
f <- readRDS(file.path(archive, x$fit_manifests$default$file))$fit
collection_truth <- b$sim$true_params$beta_theta_true
collection_truth[1, ] <- collection_truth[1, ] +
  f$infos$list_X_theta_mat$mean_df * collection_truth[2, ]
collection_truth[2, ] <- f$infos$list_X_theta_mat$sd_df * collection_truth[2, ]
samples <- f$infos$data_info[!duplicated(f$infos$data_info$Sample), ]
close(unname(f$X_theta %*% collection_truth),
      unname(cbind(1, samples$X_theta) %*% b$sim$true_params$beta_theta_true))
for (i in seq_len(nrow(x$collection))) {
  a <- x$collection[i, ]
  cov <- match(a$covariate, colnames(f$X_theta))
  s <- match(a$term, sp)
  draws <- f$results_output$beta_theta_output[cov, s, , ]
  close(as.numeric(a[c("estimate", "lower", "upper", "truth")]),
        c(mean(draws), quantile(draws, c(.025, .975)), collection_truth[cov, s]))
}
theta <- plogis(b$sim$true_params$beta_theta_true[1, ] +
                  f$infos$list_X_theta_mat$mean_df * b$sim$true_params$beta_theta_true[2, ])
p <- b$params$p * pnorm(log(1.5), b$params$mu1, b$params$sigma1, lower.tail = FALSE)
for (i in seq_len(nrow(x$detection_effort))) {
  a <- x$detection_effort[i, ]
  expected <- sum(vapply(seq_along(sp), function(s) {
    sum(vapply(0:a$M, function(collected) {
      missed <- prod((1 - p[, s])^(a$K * collected))
      dbinom(collected, a$M, theta[s]) * (1 - missed)
    }, numeric(1)))
  }, numeric(1)))
  close(a$truth, expected)
}
# Independently reconstruct fitted expected-detection draws at both extremes.
for (i in c(1L, nrow(x$detection_effort))) {
  a <- x$detection_effort[i, ]
  expected_draws <- 0
  for (s in seq_along(sp)) {
    theta_draws <- plogis(c(f$results_output$beta_theta_output[1, s, , ]))
    missed_draws <- (1 - c(f$results_output$p_output[1, s, , ]))^a$K *
      (1 - c(f$results_output$p_output[2, s, , ]))^a$K
    expected_draws <- expected_draws + 1 - (1 - theta_draws * (1 - missed_draws))^a$M
  }
  close(as.numeric(a[c("estimate", "lower", "upper")]),
        c(mean(expected_draws), quantile(expected_draws, c(.025, .975))))
}
stopifnot(all(x$detection_effort$K <= 6L),
          all(x$detection_effort$lower <= x$detection_effort$estimate),
          all(x$detection_effort$upper >= x$detection_effort$estimate))
cat("Trait decomposition, sampling-effort truth, source hashes and fit hashes verified.\n")
