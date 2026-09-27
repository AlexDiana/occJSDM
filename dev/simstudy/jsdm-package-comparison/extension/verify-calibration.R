# Independent arithmetic and archive checks; does not source calibration helpers.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
bundle <- readRDS(args[1]); root <- normalizePath(args[2])
stopifnot(bundle$complete, nrow(bundle$prediction) == 626L,
  nrow(bundle$manifest) == 640L, all(bundle$manifest$completed),
  sum(!bundle$manifest$fit_ok) == 14L,
  !anyDuplicated(bundle$prediction$job),
  all(bundle$prediction$job %in% bundle$manifest$job[bundle$manifest$scored]))
for (h in bundle$provenance$input_hashes) stopifnot(identical(tools::md5sum(names(h)), h))
stopifnot(identical(tools::md5sum(names(bundle$provenance$sources)), bundle$provenance$sources))

# A newly recorded hash only detects subsequent changes. Reconstruct the
# prespecified experiment as well, before trusting the archived truth files.
simulator <- new.env(parent = globalenv())
sys.source(file.path(root, "worker-normal-integral-v1/extension/study.R"), envir = simulator)
for (dataset in unique(bundle$prediction$dataset)) {
  j <- bundle$prediction[match(dataset, bundle$prediction$dataset), ]
  declared <- simulator$simulate_community(j$replicate, j$scenario)
  archived <- readRDS(file.path(root, "truth", paste0(dataset, ".rds")))
  stopifnot(isTRUE(all.equal(declared, archived, tolerance = 0)))
}

for (k in seq_len(nrow(bundle$prediction))) {
  j <- bundle$prediction[k, ]
  s <- readRDS(file.path(root, "jobs", j$job, "score.rds"))
  e <- as.vector(s$estimate) - as.vector(s$truth)
  stopifnot(abs(sum(e)/length(e) - j$bias) < 1e-12,
    abs(sqrt(sum(e*e)/length(e)) - j$rmse) < 1e-12,
    abs(sum(abs(e))/length(e) - j$mae) < 1e-12)
  grid_rows <- subset(bundle$rows, job == j$job & target == "Marginal probability")
  t <- readRDS(file.path(root, "truth", paste0(j$dataset, ".rds")))
  input <- readRDS(file.path(root, "inputs", paste0(j$input, ".rds")))
  declared_input <- simulator$make_input(t, j$n_sites, j$n_species, j$response, j$use_traits)
  stopifnot(isTRUE(all.equal(input, declared_input, tolerance = 0)),
    identical(s$observations, t$y[301:600, seq_len(j$n_species), drop = FALSE]))
  # Recompute test truth once per dataset; both point and grid analyses must
  # refer to the same community. No fitted predictions enter this calculation.
  truth_key <- paste0("test_", j$dataset)
  if (!exists(truth_key, envir = simulator, inherits = FALSE)) {
    assign(truth_key, simulator$marginal_truth(t, t$raw_x[301:600, ]), envir = simulator)
  }
  declared_test_truth <- get(truth_key, envir = simulator)[, seq_len(j$n_species), drop = FALSE]
  stopifnot(identical(dim(s$truth), dim(declared_test_truth)),
            max(abs(s$truth - declared_test_truth)) < 1e-9)
  # Independent truth calculation, including the curved response before integration.
  for (species in seq_len(j$n_species)) for (site in seq_len(5L)) {
    env <- unlist(bundle$grid[site, ])
    mu <- t$beta[1,species] + sum(env*t$beta[-1,species]) + env[1]^2*t$curvature[species]
    variance <- sum(t$loading[,species]^2)
    truth <- integrate(function(z) plogis(mu + sqrt(variance)*z)*dnorm(z),
                       -Inf, Inf, rel.tol=1e-9)$value
    stopifnot(abs(grid_rows$truth[(species-1L)*5L+site] - truth) < 1e-9)
  }
  if (j$package %in% c("gllvm","sjSDM")) {
    stopifnot(all(is.na(grid_rows$lower)), all(is.na(grid_rows$upper)))
  }
  coefficient_rows <- subset(bundle$rows, job == j$job & target == "Environmental coefficient")
  if (j$package == "Hmsc" || (j$scenario == "curved" && j$response == "linear")) {
    stopifnot(nrow(coefficient_rows) == 0L)
  } else {
    r <- readRDS(file.path(root,"jobs",j$job,"result.rds"))
    input <- readRDS(file.path(root,"inputs",paste0(j$input,".rds")))
    shape <- dim(r$parameters$beta)
    # Matrix multiplication applies the affine change of predictor basis to all draws.
    transform <- diag(shape[1]); transform[-1,-1] <- diag(1/input$spread)
    transform[1,-1] <- -input$centre/input$spread
    raw <- transform %*% matrix(r$parameters$beta, nrow=shape[1])
    draws <- matrix(raw, nrow=shape[1]*shape[2])
    stopifnot(max(abs(rowMeans(draws) - coefficient_rows$estimate)) < 1e-10)
    if (j$package == "occJSDM") {
      for (idx in seq_len(nrow(draws))) {
        q <- quantile(draws[idx,], c(.025,.975))
        stopifnot(abs(q[1]-coefficient_rows$lower[idx]) < 1e-10,
                  abs(q[2]-coefficient_rows$upper[idx]) < 1e-10)
      }
    }
  }
}

keys <- c("package","scenario","n_sites","n_species","response","use_traits","target","term")
for (k in seq_len(nrow(bundle$summary))) {
  s <- bundle$summary[k, ]
  keep <- rep(TRUE, nrow(bundle$rows))
  for (key in keys) keep <- keep & bundle$rows[[key]] == s[[key]]
  if (s$population == "Passed diagnostics only") keep <- keep & bundle$rows$diagnostic_pass
  d <- bundle$rows[keep, ]
  e <- d$estimate-d$truth
  mean_errors <- tapply(e, d$replicate, mean)
  mean_squares <- tapply(e^2, d$replicate, mean)
  stopifnot(abs(mean(mean_errors)-s$bias) < 1e-12,
            abs(sqrt(mean(mean_squares))-s$rmse) < 1e-12,
            length(mean_errors) == s$communities,
            s$elements == nrow(d))
  available <- is.finite(d$lower) & is.finite(d$upper) & d$lower <= d$upper
  stopifnot(sum(available) == s$intervals)
  if (!any(available)) {
    stopifnot(is.na(s$coverage), is.na(s$width), s$interval_communities == 0L)
  } else {
    a <- d[available, ]
    cover <- tapply(as.numeric(a$truth >= a$lower & a$truth <= a$upper), a$replicate, mean)
    width <- tapply(a$upper-a$lower, a$replicate, mean)
    stopifnot(abs(mean(cover)-s$coverage) < 1e-12, abs(mean(width)-s$width) < 1e-12,
              length(cover) == s$interval_communities)
    if (length(cover)>1L) stopifnot(abs(sd(cover)/sqrt(length(cover))-s$coverage_mcse)<1e-12)
  }
}
stopifnot(all(bundle$integration$endpoint_change[bundle$integration$package == "occJSDM"] <= 1e-4))
cat("Verified all 626 point scores, archive/source hashes, grid truths, coefficient transformations and", nrow(bundle$summary), "community-weighted summaries. No fits run.\n")
