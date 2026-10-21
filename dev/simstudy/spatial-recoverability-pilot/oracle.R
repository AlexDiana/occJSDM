# Oracle wrappers for the recoverability pilot. Truth 2 uses the sweep's ellipse.cpp
# unmodified; truth 1 uses block-ellipse.cpp (one ellipse angle per block; see the spec's
# revision record for 10 October, evening). joint-ellipse.cpp, which the first truth 1 run
# used, is kept and compiled for the record. All know every parameter but the fields.
# Requires kernel.R, seeds.R and generator.R to be sourced first.
ELLIPSE_MD5 <- "8154a611d1b7ab926eb32bf441e9af3e"

pilot_sampler_files <- function(repo = ".") c(
  ellipse = file.path(repo, "dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp"),
  joint = file.path(repo, "dev/simstudy/spatial-recoverability-pilot/joint-ellipse.cpp"),
  block = file.path(repo, "dev/simstudy/spatial-recoverability-pilot/block-ellipse.cpp"))

compile_pilot_samplers <- function(repo = ".", cache = NULL) {
  files <- pilot_sampler_files(repo)
  if (unname(tools::md5sum(files[["ellipse"]])) != ELLIPSE_MD5) stop("ellipse.cpp has changed; the sweep sampler must be sourced unmodified")
  if (is.null(cache)) cache <- file.path(tempdir(), "recoverability-pilot-cpp")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  for (f in files) Rcpp::sourceCpp(f, cacheDir = cache, env = globalenv())
  invisible(tools::md5sum(files))
}

training_offset <- function(cm, idx) sweep(cm$environment[idx] %o% cm$B, 2, cm$B0, "+")

factor_record <- function(f, range, what) data.frame(what = what, range = range, route = f$route,
                                                     condition = f$condition, rank = f$rank, error = f$error,
                                                     stringsAsFactors = FALSE)

# draws_full: n x niter x nchain for one species. Returns niter x nchain x 3 traces.
trace_matrix <- function(draws_full, offset, truth_centred) {
  ni <- dim(draws_full)[2]; nc <- dim(draws_full)[3]
  out <- array(NA_real_, c(ni, nc, 3))
  tss <- sum(truth_centred^2)
  for (ch in seq_len(nc)) {
    m <- draws_full[, , ch]
    out[, ch, 1] <- colMeans(m)
    out[, ch, 2] <- if (tss > 0) colSums(m * truth_centred) / tss else 0
    out[, ch, 3] <- colMeans(plogis(offset + m))
  }
  out
}

run_species_oracle <- function(cm, arrangement, species, nburn = 3000L, niter = 4000L, thin = 10L, seed, nchain = 4L) {
  idx <- cm$index[[arrangement]]; xy <- cm$points[idx, , drop = FALSE]; n <- length(idx)
  f <- stable_factor(sq_exp(xy, xy, cm$ranges[species]))
  offset <- training_offset(cm, idx)[, species]
  set.seed(seed)
  initial <- cbind(rep(0, n), f$lower %*% matrix(rnorm(n * (nchain - 1L)), n, nchain - 1L))
  fit <- ellipse_draws_cpp(f$lower, cm$z[idx, species], 1L, offset, FALSE, initial, nburn, niter)
  truth <- cm$spatial[idx, species]; tc <- truth - mean(truth)
  traces <- trace_matrix(fit$draws, offset, tc)
  keep <- seq(1L, niter, by = thin)
  site_draws <- list(matrix(fit$draws[, keep, , drop = FALSE], n, length(keep) * nchain))
  list(site_draws = site_draws, traces = traces, offset = offset,
       factor_record = factor_record(f, cm$ranges[species], paste0("prior-", arrangement, "-", cm$species$species[species])),
       mean_proposals = fit$mean_proposals, nburn = nburn, niter = niter, thin = thin, seed = seed)
}

joint_weights <- function(cm) sqrt(cm$shared_fraction) * cm$loadings

run_joint_oracle <- function(cm, arrangement, nburn = 3000L, niter = 4000L, thin = 10L, seed, nchain = 4L) {
  stopifnot(cm$truth == "t1")
  idx <- cm$index[[arrangement]]; xy <- cm$points[idx, , drop = FALSE]; n <- length(idx); S <- nrow(cm$species)
  f <- stable_factor(sq_exp(xy, xy, cm$range))
  offset <- training_offset(cm, idx)
  w <- joint_weights(cm); c_indep <- sqrt(1 - cm$shared_fraction)
  truth <- cm$spatial[idx, , drop = FALSE]; tc <- sweep(truth, 2, colMeans(truth), "-")
  blocks <- 2L + S; d <- blocks * n
  set.seed(seed)
  initial <- cbind(rep(0, d), do.call(rbind, replicate(blocks, f$lower %*% matrix(rnorm(n * (nchain - 1L)), n, nchain - 1L), simplify = FALSE)))
  fit <- block_ellipse_draws_cpp(f$lower, cm$z[idx, , drop = FALSE], offset, w, 2L, c_indep, tc, initial, nburn, niter, thin)
  nkeep <- dim(fit$draws)[2]
  site_draws <- lapply(seq_len(blocks), function(b) matrix(fit$draws[(b - 1L) * n + seq_len(n), , , drop = FALSE], n, nkeep * nchain))
  list(site_draws = site_draws, traces = fit$traces, offset = offset,
       factor_record = factor_record(f, cm$range, paste0("prior-", arrangement, "-joint")),
       mean_proposals = fit$mean_proposals, block_proposals = fit$block_proposals,
       block_mean_abs_angle = fit$block_mean_abs_angle, nburn = nburn, niter = niter, thin = thin, seed = seed)
}

screen_traces <- function(traces, rhat_max = 1.01, ess_min = 400) {
  q <- dim(traces)[3]
  labels <- c("field_average", "field_truth_projection", "probability_average")
  rows <- lapply(seq_len(q), function(j) {
    m <- traces[, , j]
    constant <- any(apply(m, 2, function(z) all(z == z[1])))
    rhat <- if (constant) NA_real_ else suppressWarnings(posterior::rhat(m))
    bulk <- if (constant) NA_real_ else suppressWarnings(posterior::ess_bulk(m))
    tail <- if (constant) NA_real_ else suppressWarnings(posterior::ess_tail(m))
    data.frame(quantity = labels[(j - 1L) %% 3L + 1L], species = (j - 1L) %/% 3L + 1L,
               rhat = rhat, ess_bulk = bulk, ess_tail = tail,
               pass = is.finite(rhat) && is.finite(bulk) && is.finite(tail) && rhat <= rhat_max && bulk >= ess_min && tail >= ess_min,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
