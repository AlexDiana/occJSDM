#!/usr/bin/env Rscript
# Read only completed result/fit pairs. Cache small derived records so
# repeated runs never need to reread unchanged large fits. No MCMC.

postprocessing_version <- 1L

audit_trace <- function(x, iterations, chains) {
  trace <- matrix(x, iterations, chains)
  if (all(trace == trace[1])) return(c(rhat = NA_real_, ess_mean = NA_real_,
                                      mcse = 0, chain_gap = 0))
  c(rhat = posterior::rhat(trace), ess_mean = posterior::ess_mean(trace),
    mcse = posterior::mcse_mean(trace),
    chain_gap = if (chains > 1L) diff(range(colMeans(trace))) else NA_real_)
}

derive_completed_fit_audit <- function(result, fit) {
  prior <- result$job$prior_name
  key <- result$job$key
  scenario <- sub(paste0("-", prior, "-[0-9]+$"), "", key)
  replicate <- as.integer(sub(".*-", "", key))
  identity <- data.frame(key = key, scenario = scenario, prior = prior,
                         replicate = replicate, model = fit$infos$model)
  truth_for <- function(metric, n) {
    rows <- result$elements[result$elements$metric == metric, ]
    rows <- rows[order(rows$element), ]
    stopifnot(nrow(rows) == n, identical(rows$element, seq_len(n)))
    rows$truth
  }
  ro <- fit$results_output
  beta <- ro$beta_theta_output
  dims <- dim(beta)
  stopifnot(length(dims) == 4L, dims[1] == 2L)
  S <- dims[2]; iterations <- dims[3]; chains <- dims[4]
  species <- as.character(fit$infos$speciesNames)
  stopifnot(length(species) == S)
  intercept <- matrix(beta[1, , , , drop = FALSE], nrow = S)
  slope <- matrix(beta[2, , , , drop = FALSE], nrow = S)
  intercept_truth <- truth_for("collection_intercept", S)
  slope_truth <- truth_for("collection_slope", S)
  sign_name <- c("negative", "zero", "positive")[sign(slope_truth) + 2L]

  # Apply the probability transform to every joint posterior draw. The
  # contrast uses paired probabilities from the same draw, retaining their
  # posterior dependence, rather than plugging in posterior mean slopes.
  point_draws <- list(theta_at_minus1SD = plogis(intercept - slope),
                      theta_at_mean = plogis(intercept),
                      theta_at_plus1SD = plogis(intercept + slope))
  point_truth <- list(theta_at_minus1SD = plogis(intercept_truth - slope_truth),
                      theta_at_mean = plogis(intercept_truth),
                      theta_at_plus1SD = plogis(intercept_truth + slope_truth))
  point_draws$theta_change_mean_to_plus1SD <-
    point_draws$theta_at_plus1SD - point_draws$theta_at_mean
  point_truth$theta_change_mean_to_plus1SD <-
    point_truth$theta_at_plus1SD - point_truth$theta_at_mean
  per_species <- list(); groups <- list()
  for (metric in names(point_draws)) {
    draws <- point_draws[[metric]]
    truth <- point_truth[[metric]]
    estimates <- rowMeans(draws)
    quantiles <- t(apply(draws, 1L, quantile, c(.025, .975), names = FALSE))
    per_species[[metric]] <- cbind(identity,
      data.frame(metric = metric, species_index = seq_len(S), species = species,
                 slope_sign = sign_name, standardized_slope_truth = slope_truth,
                 truth = truth, estimate = estimates, bias = estimates - truth,
                 bias_percentage_points = 100 * (estimates - truth),
                 posterior_lower = quantiles[, 1], posterior_upper = quantiles[, 2]))
    for (group in c("all", "negative", "zero", "positive")) {
      use <- if (group == "all") seq_len(S) else which(sign_name == group)
      if (!length(use)) next
      trace <- colMeans(draws[use, , drop = FALSE])
      errors <- estimates[use] - truth[use]
      bounds <- quantile(trace, c(.025, .975), names = FALSE)
      groups[[paste(metric, group)]] <- cbind(identity,
        data.frame(metric = metric, slope_sign = group, n_species = length(use),
                   truth = mean(truth[use]), estimate = mean(estimates[use]),
                   bias = mean(errors), bias_percentage_points = 100 * mean(errors),
                   mae = mean(abs(errors)), rmse = sqrt(mean(errors^2)),
                   posterior_lower = bounds[1], posterior_upper = bounds[2],
                   as.list(audit_trace(trace, iterations, chains))))
    }
  }

  ordering_pairs <- NULL; ordering_summary <- NULL
  if (fit$infos$model == "two_stage") {
    p <- ro$p_output; q <- ro$q_output
    stopifnot(identical(dim(p), dim(q)), length(dim(p)) == 4L,
              dim(p)[2] == S, dim(p)[3] == iterations, dim(p)[4] == chains)
    P <- dim(p)[1]; n_pairs <- P * S
    pmat <- matrix(p, nrow = n_pairs); qmat <- matrix(q, nrow = n_pairs)
    reversed <- pmat <= qmat
    fraction <- rowMeans(reversed)
    chain_fraction <- vapply(seq_len(chains), function(ch) {
      cols <- ((ch - 1L) * iterations + 1L):(ch * iterations)
      rowMeans(reversed[, cols, drop = FALSE])
    }, numeric(n_pairs))
    chain_fraction <- matrix(chain_fraction, n_pairs, chains)
    p_true <- truth_for("p", n_pairs); q_true <- truth_for("q", n_pairs)
    primers <- as.character(fit$infos$primerNames)
    stopifnot(length(primers) == P)
    pair_labels <- expand.grid(primer_index = seq_len(P), species_index = seq_len(S))
    ordering_pairs <- cbind(identity, pair_labels,
      data.frame(primer = primers[pair_labels$primer_index],
                 species = species[pair_labels$species_index],
                 p_effective_truth = p_true, q_effective_truth = q_true,
                 true_p_le_q = p_true <= q_true, true_p_minus_q = p_true - q_true,
                 p_posterior_mean = rowMeans(pmat), q_posterior_mean = rowMeans(qmat),
                 probability_p_le_q = fraction,
                 maximum_chain_probability_p_le_q = apply(chain_fraction, 1L, max),
                 chain_probability_gap = apply(chain_fraction, 1L, function(x) diff(range(x))),
                 majority_reversed = fraction > .5,
                 any_chain_majority_reversed = apply(chain_fraction > .5, 1L, any)))
    chain_records <- as.data.frame(chain_fraction)
    names(chain_records) <- paste0("chain", seq_len(chains), "_probability_p_le_q")
    ordering_pairs <- cbind(ordering_pairs, chain_records)
    flag <- if (any(fraction > .5)) "possible_reversed_detection_mode" else
      if (any(chain_fraction > .5)) "possible_chain_specific_reversed_mode" else
        "no_majority_reversed_pair_or_chain"
    ordering_summary <- cbind(identity,
      data.frame(n_pairs = n_pairs, n_draws_per_pair = ncol(pmat),
                 overall_probability_p_le_q = mean(reversed),
                 maximum_pair_probability_p_le_q = max(fraction),
                 maximum_chain_pair_probability_p_le_q = max(chain_fraction),
                 pairs_majority_reversed = sum(fraction > .5),
                 pairs_any_chain_majority_reversed = sum(apply(chain_fraction > .5, 1L, any)),
                 true_pairs_p_le_q = sum(p_true <= q_true),
                 minimum_true_p_minus_q = min(p_true - q_true),
                 ordering_flag = flag))
  }
  list(identity = identity, collection_species = do.call(rbind, per_species),
       collection_groups = do.call(rbind, groups),
       ordering_pairs = ordering_pairs, ordering_summary = ordering_summary)
}

# Hand-checked small arrays verify draw ordering and the probability contrast
# without any data generation or calls to a model-fitting function.
self_check_completed_fit_audit <- function() {
  beta <- array(0, c(2, 3, 2, 2))
  beta[2, , , ] <- array(rep(c(1, -1, 0), 4), c(3, 2, 2))
  fit <- list(infos = list(model = "two_stage", speciesNames = c("a", "b", "c"),
                           primerNames = "one"),
              results_output = list(beta_theta_output = beta,
                p_output = array(rep(c(.8, .2, .6), 4), c(1, 3, 2, 2)),
                q_output = array(rep(c(.1, .3, .4), 4), c(1, 3, 2, 2))))
  result <- list(job = list(key = "synthetic-default-01", prior_name = "default"),
                 elements = data.frame(metric = rep(c("collection_intercept", "collection_slope", "p", "q"), each = 3),
                                       element = rep(1:3, 4), truth = c(rep(0, 3), 1, -1, 0, rep(.8, 3), rep(.1, 3))))
  audit <- derive_completed_fit_audit(result, fit)
  stopifnot(audit$ordering_summary$overall_probability_p_le_q == 1/3,
            audit$ordering_summary$maximum_pair_probability_p_le_q == 1,
            audit$ordering_summary$pairs_majority_reversed == 1L,
            audit$ordering_summary$true_pairs_p_le_q == 0L)
  contrasts <- audit$collection_species[audit$collection_species$metric == "theta_change_mean_to_plus1SD", ]
  stopifnot(max(abs(contrasts$truth - c(.2310585786300049, -.2310585786300049, 0))) < 1e-14,
            all(contrasts$bias == 0),
            identical(contrasts$slope_sign, c("positive", "negative", "zero")))
  # Mixing intercept draws at -2 and +2 changes the mean probability
  # contrast. Plugging the mean intercept into plogis would incorrectly
  # reproduce .23106 rather than this independently calculated .11076.
  fit$results_output$beta_theta_output[1, 1, , ] <- matrix(c(2, -2, 2, -2), 2)
  mixed <- derive_completed_fit_audit(result, fit)$collection_species
  estimate <- mixed$estimate[mixed$metric == "theta_change_mean_to_plus1SD" &
                               mixed$species_index == 1L]
  stopifnot(abs(estimate - .110757774096214) < 1e-14)
}

run_completed_fit_audit <- function(root, results_dir = file.path(root, c("main-results", "operational-k6")),
                                    output_dir = file.path(root, "q-audit", "completed-fits")) {
  results_dir <- normalizePath(results_dir, mustWork = TRUE)
  dir.create(file.path(output_dir, "cache"), recursive = TRUE, showWarnings = FALSE)
  self_check_completed_fit_audit()
  # Inventory once: files completing after this snapshot are picked up on
  # the next run. A fit without a completed result marker is never opened.
  results <- list.files(results_dir, pattern = "-result\\.rds$", full.names = TRUE)
  output <- list(); new <- 0L; reused <- 0L
  for (result_file in results) {
    key <- sub("-result\\.rds$", "", basename(result_file))
    fit_file <- sub("-result\\.rds$", "-fit.rds", result_file)
    stopifnot(file.exists(fit_file))
    info <- file.info(c(result_file, fit_file))
    signature <- list(paths = normalizePath(c(result_file, fit_file)),
                       sizes = info$size, mtimes = as.numeric(info$mtime))
    cache_file <- file.path(output_dir, "cache", paste0(key, "-audit.rds"))
    cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
    if (!is.null(cached) && identical(cached$version, postprocessing_version) &&
        identical(cached$signature, signature)) {
      audit <- cached$audit
      reused <- reused + 1L
    } else {
      result <- readRDS(result_file)
      fit <- readRDS(fit_file)$fit
      stopifnot(identical(result$job$key, key))
      audit <- derive_completed_fit_audit(result, fit)
      rm(fit, result)
      temporary <- paste0(cache_file, ".tmp")
      saveRDS(list(version = postprocessing_version, signature = signature, audit = audit), temporary)
      stopifnot(file.rename(temporary, cache_file))
      new <- new + 1L
    }
    output[[key]] <- audit
  }
  for (name in c("ordering_summary", "ordering_pairs", "collection_groups", "collection_species")) {
    combined <- do.call(rbind, lapply(output, `[[`, name))
    if (is.null(combined)) combined <- data.frame()
    rownames(combined) <- NULL
    write.csv(combined, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  }
  status <- data.frame(time = format(Sys.time(), tz = "UTC", usetz = TRUE),
                        result_directory = paste(results_dir, collapse = ";"), completed_snapshot = length(results),
                        newly_processed = new, reused_cache = reused)
  write.csv(status, file.path(output_dir, "status.csv"), row.names = FALSE)
  print(status)
  ordering <- do.call(rbind, lapply(output, `[[`, "ordering_summary"))
  if (!is.null(ordering)) {
    cat("Two-stage fits:", nrow(ordering), "\n")
    cat("Maximum pair Pr(p <= q):", max(ordering$maximum_pair_probability_p_le_q), "\n")
    cat("Pairs majority reversed:", sum(ordering$pairs_majority_reversed), "\n")
    flagged <- ordering[ordering$pairs_majority_reversed > 0 | ordering$pairs_any_chain_majority_reversed > 0, ]
    if (nrow(flagged)) print(flagged)
  }
  invisible(status)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
  .libPaths(c(file.path(root, "library"), .libPaths()))
  run_completed_fit_audit(root,
    results_dir = if (length(args) >= 2L) strsplit(args[2], ",", fixed = TRUE)[[1]] else file.path(root, c("main-results", "operational-k6")),
    output_dir = if (length(args) >= 3L) args[3] else file.path(root, "q-audit", "completed-fits"))
}
