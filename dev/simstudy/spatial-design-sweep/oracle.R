# Oracle for the Lesson 2 sweep: the diagnosis elliptical-slice sampler,
# sourced unmodified, estimating only the spatial field with every other
# parameter supplied at its generating value. See PLAN.md.

oracle_source_files <- function(repo) {
  file.path(repo, c("dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp",
                    "dev/simstudy/spatial-amplitude-prior/robust.R"))
}

oracle_source_hashes <- function(repo) tools::md5sum(oracle_source_files(repo))

compile_oracle <- function(repo, cache = NULL) {
  files <- oracle_source_files(repo)
  source(files[2])
  if (is.null(cache)) cache <- file.path(tempdir(), "sweep-oracle-cpp")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(files[1], cacheDir = cache, env = globalenv())
  invisible(oracle_source_hashes(repo))
}

oracle_covariance <- function(xy, range, jitter = 1e-8) {
  K <- sq_exp(xy, xy, range) + diag(jitter, nrow(xy))
  (K + t(K)) / 2
}

run_oracle <- function(input, arrangement, species, seed, nburn = 1000L, niter = 2000L, nchain = 4L) {
  tr <- input$surveys[[arrangement]]$truth; land <- input$landscape
  lower <- t(chol(oracle_covariance(tr$xy, land$range)))
  offset <- land$B0[species] + land$B[species] * tr$environment
  n <- nrow(lower)
  set.seed(seed)
  initial <- cbind(rep(0, n), lower %*% matrix(rnorm(n * (nchain - 1L)), n, nchain - 1L))
  fit <- ellipse_draws_cpp(lower, tr$z[, species], 1L, offset, FALSE, initial, nburn, niter)
  list(draws = fit$draws, mean_proposals = fit$mean_proposals, offset = offset, lower = lower)
}

score_oracle <- function(draws, truth_field, psi, offset) {
  n <- dim(draws)[1]; ni <- dim(draws)[2]; nc <- dim(draws)[3]
  field <- matrix(draws, n, ni * nc)
  probabilities <- plogis(field + offset)
  med <- apply(field, 1, median)
  tc <- truth_field - mean(truth_field); fc <- med - mean(med)
  diagnostics <- list()
  add <- function(values, label) {
    values <- matrix(values, ni, nc)
    diagnostics[[length(diagnostics) + 1L]] <<- data.frame(quantity = label,
      as.list(robust_trace_diagnostics(values)),
      ess_mean = as.numeric(posterior::ess_mean(values)))
  }
  for (i in seq_len(n)) add(field[i, ], paste0("field_", i))
  add(colMeans(field), "field_average")
  add(colSums(field * tc) / sum(tc^2), "field_truth_projection")
  add(colMeans(probabilities), "probability_average")
  diagnostics <- do.call(rbind, diagnostics)
  flags <- robust_flag_rows(diagnostics) | !is.finite(diagnostics$ess_mean) | diagnostics$ess_mean < 100
  pmean <- rowMeans(probabilities)
  metrics <- data.frame(
    centred_rmse = sqrt(mean((fc - tc)^2)), centred_correlation = cor(fc, tc),
    centred_slope = sum(fc * tc) / sum(tc^2), raw_rmse = sqrt(mean((med - truth_field)^2)),
    zero_field_rmse = sqrt(mean(tc^2)), median_field_rms = sqrt(mean(med^2)),
    occupancy_bias = mean(pmean - psi), occupancy_mae = mean(abs(pmean - psi)),
    flag_count = sum(flags), max_rhat = max(diagnostics$rhat, na.rm = TRUE),
    min_ess = min(as.matrix(diagnostics[c("ess_bulk", "ess_median", "ess_q025", "ess_q975", "ess_mean")]), na.rm = TRUE))
  list(metrics = metrics, diagnostics = diagnostics, field_median = med, probability_mean = pmean)
}

oracle_lattice <- function(draws, xy_sites, xy_lattice, range, offset_lattice, thin = 4L, seed) {
  n <- dim(draws)[1]; ni <- dim(draws)[2]; nc <- dim(draws)[3]
  Kss <- oracle_covariance(xy_sites, range)
  Kls <- sq_exp(xy_lattice, xy_sites, range)
  A <- Kls %*% solve(Kss)
  cond <- sq_exp(xy_lattice, xy_lattice, range) - A %*% t(Kls)
  cond <- (cond + t(cond)) / 2 + diag(1e-8, nrow(xy_lattice))
  L <- t(chol(cond))
  keep <- seq(1L, ni, by = thin)
  set.seed(seed)
  prob <- numeric(nrow(xy_lattice)); fmean <- numeric(nrow(xy_lattice)); count <- 0L
  for (ch in seq_len(nc)) for (it in keep) {
    f <- as.vector(A %*% draws[, it, ch]) + as.vector(L %*% rnorm(nrow(xy_lattice)))
    prob <- prob + plogis(offset_lattice + f); fmean <- fmean + f; count <- count + 1L
  }
  list(probability_mean = prob / count, field_mean = fmean / count, draws_used = count)
}
