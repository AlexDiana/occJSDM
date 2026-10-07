# Independent checks of the fitted Bernoulli estimand and known-w Beta update.
# Does not fit a full occupancy model or change the package.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1]) else normalizePath("..")
out <- file.path(root, "q-audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(root, "library"), .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
ns <- asNamespace("occJSDM")
source(file.path(root, "source-main/tests/testthat/helper-simstudy.R"))
writeLines(c(paste("Source revision:", readLines(file.path(root, "source-revision.txt"))),
             paste("Installed package:", find.package("occJSDM")),
             capture.output(sessionInfo())), file.path(out, "session-info.txt"))

# For positive threshold t and rounded nonnegative read counts, a read is
# detected iff X >= log(ceiling(t) + .5), except a probability-zero tie.
retention <- function(mu, sigma, threshold = 1) {
  pnorm(log(ceiling(threshold) + .5), mu, sigma, lower.tail = FALSE)
}
ret <- data.frame(stage = c("p", "q"), mu = c(5, 1.5), sigma = 1)
ret$retention <- mapply(retention, ret$mu, ret$sigma)
write.csv(ret, file.path(out, "threshold-retention.csv"), row.names = FALSE)

# The former high-K report used nominal read-event probabilities as its
# truth. Only this aggregate mean-error correction can be reconstructed
# from the rounded report; element-wise coverage needs the original rows.
historical <- data.frame(cell = c("qnear_K3", "qnear_K30", "qfar_K3", "qfar_K30"),
                         nominal_mean = c(.0298, .0298, .2244, .2244),
                         reported_nominal_bias = c(-.0001, -.004, -.0365, -.0314))
historical$effective_mean <- historical$nominal_mean * ret$retention[2]
historical$threshold_difference <- historical$effective_mean - historical$nominal_mean
historical$approx_effective_bias <- historical$reported_nominal_bias - historical$threshold_difference
write.csv(historical, file.path(out, "historical-rounded-bias-correction.csv"), row.names = FALSE)

# Independently repeat the public generator's Bernoulli -> lognormal ->
# rounded read count mechanism. This tests the probability being compared,
# not a posterior fit, and does not call the package simulator.
set.seed(914071L)
Nmc <- 1000000L
mc_grid <- expand.grid(nominal = c(.03, .2244, .30, .60), threshold = c(1, 3))
mc_grid$mu <- ifelse(mc_grid$nominal < .25, 1.5, 5)
mc_grid$sigma <- 1
mc_rows <- lapply(seq_len(nrow(mc_grid)), function(i) {
  a <- mc_grid[i, ]
  event <- rbinom(Nmc, 1, a$nominal)
  x <- numeric(Nmc)
  x[event == 1] <- rnorm(sum(event), a$mu, a$sigma)
  y <- pmax(0, round(exp(x) - 1))
  actual <- mean(y >= a$threshold)
  expected <- a$nominal * retention(a$mu, a$sigma, a$threshold)
  se <- sqrt(expected * (1 - expected) / Nmc)
  cbind(a, observations = Nmc, empirical = actual, expected = expected,
        standard_error = se, discrepancy_se = (actual - expected) / se)
})
mc <- do.call(rbind, mc_rows)
write.csv(mc, file.path(out, "empirical-threshold-check.csv"), row.names = FALSE)
stopifnot(max(abs(mc$discrepancy_se)) < 5)

# Count using actual row labels, independently of createDataIdx(). Return
# posterior shapes for known sample state w and thresholded observations.
known_w_shapes <- function(y, missing, w, sample, primer, maxP, priors) {
  rows <- list()
  for (s in seq_len(ncol(w))) for (l in seq_len(maxP)) {
    use <- !missing[, s] & primer == l
    wrep <- w[sample, s]
    for (stage in c("p", "q")) {
      pick <- use & wrep == as.integer(stage == "p")
      n <- sum(pick)
      success <- sum(y[pick, s])
      alpha <- priors[paste0("a_", stage)] + success
      beta <- priors[paste0("b_", stage)] + n - success
      rows[[length(rows) + 1L]] <- data.frame(species = s, primer = l,
        stage = stage, n = n, successes = success,
        alpha = unname(alpha), beta = unname(beta),
        posterior_mean = unname(alpha / (alpha + beta)),
        posterior_variance = unname(alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))))
    }
  }
  do.call(rbind, rows)
}

# Construct irregular global primer identities, missing observations, and
# interleaved sample rows. The p/q updater must use the supplied mappings.
set.seed(914072L)
S <- 3L
sample <- as.integer(c(4,1,6,2,1,3,6,2,4,5,5,6,3,1,4,2,6,5,1,4,2,5,6,3))
primer <- as.integer(c(3,1,3,1,2,2,1,3,1,2,3,2,3,1,3,1,1,3,2,2,3,1,3,2))
w <- matrix(as.double(c(1,0,1,0,0,1, 0,1,1,0,1,0, 1,1,0,0,1,0)), 6, S)
y <- matrix(rbinom(length(sample)*S, 1, .4), length(sample), S)
missing <- matrix(FALSE, nrow(y), ncol(y))
missing[cbind(c(2,7,11,20), c(1,2,1,3))] <- TRUE
# An entirely missing species/primer cell must return its prior.
missing[primer == 2, 3] <- TRUE
cimk <- y * (2 - w[sample, , drop = FALSE])
cimk[missing] <- NA_real_
mask <- missing * 1L
priors <- c(a_p = 5, b_p = 1, a_q = 1, b_q = 20)
shape <- known_w_shapes(y, missing, w, sample, primer, 3L, priors)
sampler <- get("sample_pq_cpp_parallel", ns)
draw_cpp <- function() sampler(cimk, mask, w, primer, sample, 3L,
                              priors[1], priors[2], priors[3], priors[4])
set.seed(914073L)
cpp <- draw_cpp()
seed_cpp <- .Random.seed
set.seed(914073L)
r_draw <- rbeta(nrow(shape), shape$alpha, shape$beta)
seed_r <- .Random.seed
cpp_order <- unlist(lapply(seq_len(S), function(s) {
  unlist(lapply(seq_len(3L), function(l) c(cpp$p[l, s], cpp$q[l, s])))
}))
stopifnot(identical(unname(cpp_order), unname(r_draw)), identical(seed_cpp, seed_r))
shape$cpp_single_draw <- cpp_order
shape$independent_r_draw <- r_draw

# Distributional check is additional to the exact same-seed comparison.
# Each call is an independent conditional draw with w held fixed.
set.seed(914074L)
B <- 20000L
draws <- replicate(B, {
  a <- draw_cpp()
  unlist(lapply(seq_len(S), function(s) {
    unlist(lapply(seq_len(3L), function(l) c(a$p[l, s], a$q[l, s])))
  }))
})
shape$empirical_mean <- rowMeans(draws)
shape$mean_standard_error <- sqrt(shape$posterior_variance / B)
shape$mean_discrepancy_se <- (shape$empirical_mean - shape$posterior_mean) / shape$mean_standard_error
shape$empirical_variance <- apply(draws, 1L, var)
shape$variance_ratio <- shape$empirical_variance / shape$posterior_variance
stopifnot(max(abs(shape$mean_discrepancy_se)) < 5)
write.csv(shape, file.path(out, "known-w-conditional-sampling.csv"), row.names = FALSE)

# On each original scenario with spatial effects disabled, compute the
# exact conditional known-w reference using the public simulator. These
# references do not use fitted w and are diagnostic information ceilings.
specs <- simstudy_scenarios()
specs <- specs[vapply(specs, function(s) s$label %in% historical$cell, logical(1))]
public_rows <- list()
for (spec in specs) {
  spec$useSpatField <- FALSE
  spec$n_supportpoints <- NULL
  truth <- draw_truth(spec, seed = 914075L)
  sim <- simstudy_simulate(truth)
  y <- (sim$data_list$OTU >= 1) * 1
  inf <- sim$data_list$info
  w <- sim$true_params$w_true
  shapes <- known_w_shapes(y, is.na(y), w, as.integer(inf$Sample),
                          as.integer(inf$Primer), spec$P, priors)
  shapes$cell <- spec$label
  shapes$nominal_truth <- mapply(function(stage, l, s) truth$params[[stage]][l, s],
                                 shapes$stage, shapes$primer, shapes$species)
  shapes$effective_truth <- shapes$nominal_truth * ret$retention[match(shapes$stage, ret$stage)]
  shapes$posterior_error_effective <- shapes$posterior_mean - shapes$effective_truth
  shapes$posterior_error_nominal <- shapes$posterior_mean - shapes$nominal_truth
  shapes$known_w_expected_prior_bias <- mapply(function(stage, n, prob) {
    a <- priors[paste0("a_",stage)]; b <- priors[paste0("b_",stage)]
    (a - (a+b)*prob)/(n+a+b)
  }, shapes$stage, shapes$n, shapes$effective_truth)
  public_rows[[spec$label]] <- shapes
}
public <- do.call(rbind, public_rows)
write.csv(public, file.path(out, "public-simulator-known-w-reference.csv"), row.names = FALSE)
public_summary <- aggregate(cbind(posterior_error_effective, posterior_error_nominal,
                                  known_w_expected_prior_bias) ~ cell + stage,
                             data = public, FUN = mean)
write.csv(public_summary, file.path(out, "public-simulator-known-w-summary.csv"), row.names = FALSE)

# These are illustrative sensitivity priors, not recommended new defaults.
# Independent p ~ Beta(5,1) gives P(p<=q) = E(q^5), in closed form.
prior_grid <- data.frame(label = c("default", "q_mean_0.10", "q_mean_0.20"),
                         a_p = 5, b_p = 1, a_q = c(1,1,2), b_q = c(20,9,8))
prior_grid$q_mean <- with(prior_grid, a_q / (a_q+b_q))
prior_grid$prior_probability_p_le_q <- with(prior_grid, exp(lbeta(a_q+5,b_q)-lbeta(a_q,b_q)))
write.csv(prior_grid, file.path(out, "illustrative-prior-separation.csv"), row.names = FALSE)
cat("All independent diagnostics completed. No full-model fits were run.\n")
print(ret)
print(historical)
cat("Maximum empirical threshold discrepancy (SE):", max(abs(mc$discrepancy_se)), "\n")
cat("Same-seed conditional C++/independent-R draws and final seeds identical.\n")
cat("Maximum conditional draw mean discrepancy (SE):", max(abs(shape$mean_discrepancy_se)), "\n")
cat("Conditional draw variance ratio range:", range(shape$variance_ratio), "\n")
print(public_summary)
print(prior_grid)
