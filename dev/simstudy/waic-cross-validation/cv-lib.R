# Ten-fold site cross-validation of observed-survey predictions for the
# Lesson 4 factor-count example. Research code, not package API.
#
# A withheld site is scored with occJSDM's own computeSiteWAIC() likelihood:
# its complete observed survey, with occupancy and collection states summed
# out and the shared site factors integrated out. The scorer is applied to a
# copy of the training fit that carries the withheld sites' data, prepared
# the way runOccJSDM() prepares data, with covariates standardized by the
# training fit's constants. Only the posterior arrays come from training.

cv_mcmc <- list(nchain = 4L, nburn = 3000L, niter = 6000L, nthin = 1L)
cv_models <- c(one = 1L, two = 2L)
# The archived Lesson 4 fits used these seeds (prediction-build.R and
# helpers.R); reusing them lets the full-data refits be compared with them.
cv_full_seeds <- c(one = 20260926L, two = 20260921L)
cv_fold_assignment_seed <- 20261013L
cv_fold_seed <- function(fold, model) 2026101300L + 10L * as.integer(fold) + cv_models[[model]]

cv_job_names <- function(K = 10L) {
  c(paste0("full-", names(cv_models)),
    as.vector(outer(sprintf("fold%02d", seq_len(K)), names(cv_models), paste, sep = "-")))
}

cv_parse_job <- function(job) {
  parts <- strsplit(job, "-", fixed = TRUE)[[1]]
  stopifnot(length(parts) == 2L, parts[2] %in% names(cv_models))
  fold <- if (parts[1] == "full") NA_integer_ else as.integer(sub("^fold", "", parts[1]))
  stopifnot(parts[1] == "full" || grepl("^fold[0-9]{2}$", parts[1]))
  list(job = job, fold = fold, model = parts[2],
       seed = if (is.na(fold)) cv_full_seeds[[parts[2]]] else cv_fold_seed(fold, parts[2]))
}

cv_make_folds <- function(sites, K = 10L, seed = cv_fold_assignment_seed) {
  stopifnot(length(sites) %% K == 0L, !anyDuplicated(sites))
  set.seed(seed)
  shuffled <- sample(sites)
  data.frame(Site = shuffled, fold = rep(seq_len(K), each = length(sites) / K))
}

cv_subset <- function(data, sites) {
  keep <- data$info$Site %in% sites
  stopifnot(setequal(unique(data$info$Site[keep]), sites))
  list(info = data$info[keep, , drop = FALSE], OTU = data$OTU[keep, , drop = FALSE],
       traits = data$traits)
}

# Same arguments as the archived fits (helpers.R lesson_fit() and
# prediction-build.R), apart from the factor count.
cv_fit_args <- function(data, model, mcmc = cv_mcmc) {
  list(data = data, listParams = list(n_factors = cv_models[[model]], n_lattrait = 1L),
       threshold = 1, occCovariates = grep("^X_psi", names(data$info), value = TRUE),
       collCovariates = "X_theta", spatCovariates = NULL, MCMCparams = mcmc,
       listPriors = list(), summarisedLatentPresences = TRUE)
}

cv_fit <- function(data, model, seed, mcmc = cv_mcmc) {
  args <- cv_fit_args(data, model, mcmc)
  warnings <- character()
  set.seed(seed)
  started <- Sys.time()
  fit <- withCallingHandlers(do.call(occJSDM::runOccJSDM, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
  })
  list(fit = fit, model = model, seed = seed, mcmc = mcmc,
       args = args[setdiff(names(args), "data")], warnings = warnings,
       started = started, finished = Sys.time())
}

# Prepare `data` (any sites) for scoring under `fit`'s posterior. Mirrors the
# two-stage data preparation in runOccJSDM(); the training reproduction test
# checks that it returns exactly the fitted design for the training sites.
cv_heldout_fit <- function(fit, data) {
  info <- fit$infos
  stopifnot(info$model == "two_stage", length(info$ps) == 1L, info$ps == 0,
            identical(colnames(data$OTU), as.character(info$speciesNames)))
  di <- as.data.frame(data$info) |>
    dplyr::mutate(rownum = dplyr::row_number()) |>
    dplyr::arrange(Site, Sample, Primer)
  OTU <- data$OTU[di$rownum, , drop = FALSE]
  di <- dplyr::select(di, -rownum)
  samples <- di |> dplyr::group_by(Site, Sample) |> dplyr::slice(1) |> dplyr::ungroup()
  M_df <- samples |> dplyr::group_by(Site) |> dplyr::summarise(M = dplyr::n(), .groups = "keep")
  P_df <- di |> dplyr::group_by(Site, Sample, Primer) |> dplyr::slice(1) |>
    dplyr::group_by(Site, Sample) |> dplyr::summarise(P = dplyr::n(), .groups = "keep")
  data_K <- di |> dplyr::group_by(Site, Sample, Primer) |>
    dplyr::summarise(K = dplyr::n(), .groups = "keep") |> dplyr::ungroup()
  # Primer identity indexes the training primers, as p and q do.
  primer_id <- match(data_K$Primer, info$primerNames)
  if (anyNA(primer_id)) stop("Withheld data use a primer the training fit has not seen.")
  n <- nrow(M_df)
  list_idx <- occJSDM:::createDataIdx(n, M_df$M, P_df$P, data_K$K, TRUE, primerId_p = primer_id)
  psi_names <- info$list_X_psi_mat$names_df
  theta_names <- info$list_X_theta_mat$names_df
  site_covariates <- di |> dplyr::group_by(Site) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(psi_names), ~ dplyr::first(.x)), .groups = "drop")
  di$SiteSample <- paste(di$Site, di$Sample, sep = "-")
  pseudo <- fit
  pseudo$X_psi <- occJSDM:::transform_new_covariates(
    site_covariates[, psi_names, drop = FALSE], info$list_X_psi_mat, remove_intercept = TRUE)
  pseudo$X_theta <- occJSDM:::transform_new_covariates(
    samples[, theta_names, drop = FALSE], info$list_X_theta_mat, remove_intercept = FALSE)
  pseudo$infos$OTU <- OTU
  pseudo$infos$n <- n
  pseudo$infos$siteNames <- M_df$Site
  pseudo$infos$list_idx <- list_idx
  pseudo$infos$data_info <- di
  pseudo$cv_scoring_copy <- TRUE
  pseudo
}

# computeSiteWAIC() warns about p_waic for every set of sites it scores; for a
# held-out score that penalty is not used, so only that warning is muffled.
cv_score <- function(fit, data, draws = NULL, ...) {
  withCallingHandlers(occJSDM::computeSiteWAIC(cv_heldout_fit(fit, data), draws = draws, ...),
    warning = function(w) {
      if (grepl("p_waic > 0.4", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
    })
}

cv_log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Monte Carlo SE of a sum over sites of log posterior-mean likelihoods, by the
# delta method on the per-draw influence sum_i lik_i / mean(lik_i). Rows of
# log_lik must be all draws, ordered by iteration within chain.
cv_sum_lppd_mcse <- function(log_lik, iterations, chains) {
  stopifnot(nrow(log_lik) == iterations * chains)
  scaled <- exp(sweep(log_lik, 2, apply(log_lik, 2, max), "-"))
  influence <- rowSums(sweep(scaled, 2, colMeans(scaled), "/"))
  posterior::mcse_mean(matrix(influence, iterations, chains))
}

# Held-out summaries for one scored fold: per-site lppd overall and per chain,
# and the Monte Carlo SE of the fold's summed lppd.
cv_fold_summary <- function(score, iterations, chains) {
  ll <- score$log_lik
  stopifnot(nrow(ll) == iterations * chains,
            identical(score$draw_ids$chain, rep(seq_len(chains), each = iterations)))
  per_chain <- sapply(seq_len(chains), function(ch) {
    apply(ll[score$draw_ids$chain == ch, , drop = FALSE], 2, cv_log_mean_exp)
  })
  per_chain <- matrix(per_chain, ncol(ll), chains,
                      dimnames = list(colnames(ll), paste0("chain", seq_len(chains))))
  list(sites = data.frame(site = colnames(ll), lppd = apply(ll, 2, cv_log_mean_exp),
                          per_chain, row.names = NULL),
       sum_mcse = cv_sum_lppd_mcse(ll, iterations, chains))
}

# Posterior occupancy probabilities of the fitted sites, draws in the last
# two dimensions. These are identified even though factors are not.
cv_psi_draws <- function(fit) {
  j <- fit$results_output$jsdm_output
  n <- nrow(fit$X_psi); S <- dim(j$B0_output)[1]
  iterations <- dim(j$B0_output)[2]; chains <- dim(j$B0_output)[3]
  d <- fit$infos$n_factors
  out <- array(NA_real_, c(n * S, iterations, chains))
  for (ch in seq_len(chains)) for (it in seq_len(iterations)) {
    B <- matrix(j$B_output[, , it, ch], ncol(fit$X_psi), S)
    U <- matrix(j$U_output[, , it, ch], n, d)
    L <- matrix(j$L_output[, , it, ch], d, S)
    eta <- sweep(fit$X_psi %*% B + U %*% L, 2, j$B0_output[, it, ch], "+")
    out[, it, ch] <- stats::plogis(eta)
  }
  out
}

# The project's strict screen (interval-calibration REPORT.md): nonfinite
# Rhat, Rhat above 1.05, or bulk or tail ESS below 400, on every scored
# parameter plus the fitted sites' occupancy probabilities.
cv_screen <- function(fit, rhat_max = 1.05, ess_min = 400) {
  o <- fit$results_output
  arrays <- list(beta0_psi = o$jsdm_output$B0_output, beta_psi = o$jsdm_output$B_output,
                 beta_theta = o$beta_theta_output, p = o$p_output, q = o$q_output,
                 theta0 = o$theta0_output, psi = cv_psi_draws(fit))
  rows <- lapply(names(arrays), function(name) {
    a <- arrays[[name]]
    dims <- dim(a); k <- length(dims)
    flat <- array(a, c(prod(dims[seq_len(k - 2L)]), dims[k - 1L], dims[k]))
    stats <- t(apply(flat, 1, function(x) {
      x <- matrix(x, dims[k - 1L], dims[k])
      withCallingHandlers(
        c(rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
          ess_tail = posterior::ess_tail(x)),
        warning = function(w) {
          if (grepl("ESS has been capped", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
        })
    }))
    data.frame(param = name, element = seq_len(nrow(stats)), stats)
  })
  table <- do.call(rbind, rows)
  table$flagged <- !is.finite(table$rhat) | table$rhat > rhat_max |
    !(table$ess_bulk >= ess_min) | !(table$ess_tail >= ess_min)
  table
}
