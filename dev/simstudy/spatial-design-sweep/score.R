# Scoring for the Lesson 2 sweep full fits. Reuses the September helpers.
score_source_files <- function(repo) file.path(repo, c(
  "dev/simstudy/spatial-targeted-recheck/score.R", "dev/simstudy/spatial-design-sweep/score.R"))
score_source_hashes <- function(repo) tools::md5sum(score_source_files(repo))
local({
  here <- if (exists("repo", inherits = TRUE)) get("repo", inherits = TRUE) else "."
  source(file.path(here, "dev/simstudy/spatial-targeted-recheck/score.R"))
})

reconstruct_field_draws <- function(fit, thin = 4L) {
  js <- fit$results_output$jsdm_output
  n <- nrow(fit$Xs); S <- dim(js$B0_output)[1]; ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  H <- independent_bases(fit); keep <- seq(1L, ni, by = thin)
  out <- array(NA_real_, c(n, S, length(keep) * nc)); k <- 0L
  for (ch in seq_len(nc)) for (it in keep) {
    k <- k + 1L
    out[, , k] <- H[[js$idx_ls_output[it, ch]]] %*% matrix(js$Bs_output[, , it, ch], fit$infos$ps, S)
  }
  out
}

reconstruct_lattice_draws <- function(fit, lattice_environment, lattice_xy, thin = 4L) {
  js <- fit$results_output$jsdm_output
  S <- dim(js$B0_output)[1]; ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  X_new <- occJSDM:::transform_new_covariates(data.frame(environment = lattice_environment),
                                               fit$infos$list_X_psi_mat, remove_intercept = TRUE)
  Ks <- occJSDM:::createSpatialPredMatrix(data.frame(longitude = lattice_xy[, 1], latitude = lattice_xy[, 2]),
                                          fit$infos$l_s_grid, fit$infos$list_Xs$X_tilde, fit$infos$list_Xs_mat)
  keep <- seq(1L, ni, by = thin); m <- nrow(lattice_xy)
  with <- matrix(0, m, S); without <- matrix(0, m, S); count <- 0L
  for (ch in seq_len(nc)) for (it in keep) {
    count <- count + 1L
    base <- sweep(as.matrix(X_new) %*% matrix(js$B_output[, , it, ch], ncol(X_new), S), 2, js$B0_output[, it, ch], "+")
    field <- Ks[, , js$idx_ls_output[it, ch]] %*% matrix(js$Bs_output[, , it, ch], fit$infos$ps, S)
    with <- with + plogis(base + field); without <- without + plogis(base)
  }
  list(with = with / count, without = without / count, draws_used = count, Ks = Ks, X_new = X_new)
}

sweep_groups <- function(psi, prevalence) {
  g <- list(all = matrix(TRUE, nrow(psi), ncol(psi)), low = psi < .2, medium = psi >= .2 & psi <= .8, high = psi > .8)
  for (p in c(.05, .25, .75)) g[[paste0("prevalence_", p * 100, "pct")]] <-
    matrix(rep(abs(prevalence - p) < 1e-10, each = nrow(psi)), nrow(psi), ncol(psi))
  g
}

distance_bins <- function(lattice_xy, site_xy) {
  d <- sqrt(outer(lattice_xy[, 1], site_xy[, 1], "-")^2 + outer(lattice_xy[, 2], site_xy[, 2], "-")^2)
  nn <- apply(d, 1, min)
  list(distance = nn, bin = cut(nn, c(-Inf, .02, .05, .1, Inf),
       labels = c("up to 0.02", "0.02 to 0.05", "0.05 to 0.1", "above 0.1")))
}

score_sweep_fit <- function(fit, input, arrangement, arm, lattice_subset = 50L) {
  land <- input$landscape; tr <- input$surveys[[arrangement]]$truth
  n <- nrow(tr$xy); S <- ncol(tr$psi); knots <- 100L
  js <- fit$results_output$jsdm_output; ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  stopifnot(fit$infos$ps == knots, fit$infos$n_factors == 0L,
            identical(fit$infos$speciesNames, land$species),
            fit$infos$model == if (arm == "binary") "binary" else "two_stage",
            max(abs(fit$Xs - scale(tr$xy))) < 1e-10, max(abs(fit$X_psi - scale(tr$environment))) < 1e-10)
  env_sd <- sd(tr$environment); xy_sd <- apply(tr$xy, 2, sd)
  reconstructed <- reconstruct_spatial_draws(fit)
  x <- reconstructed$probability
  estimate <- matrix(rowMeans(matrix(x, n * S)), n, S)
  saved_difference <- verify_saved_probability(estimate, fit$results_output$psi_output)
  H <- independent_bases(fit)
  native <- suppressMessages(occJSDM:::precomputeSORmatrices(fit$infos$l_s_grid, fit$infos$list_Xs))
  basis_difference <- max(vapply(seq_along(H), function(g) max(abs(H[[g]] -
    occJSDM:::KsBproduct(native$Ks_all[, , g], diag(knots), fit$infos$list_Xs$Xs_centers))), numeric(1)))
  stopifnot(basis_difference < 1e-10)
  main <- score_draw_block(x, tr$psi, "occupancy")
  elements <- main$elements; groups <- list(main$summary)
  masks <- sweep_groups(tr$psi, land$prevalence)
  for (g in setdiff(names(masks), "all")) {
    idx <- which(masks[[g]]); if (!length(idx)) next
    el <- elements[idx, ]; trace <- apply(x[idx, , , drop = FALSE], c(2, 3), mean)
    groups[[length(groups) + 1L]] <- data.frame(metric = "occupancy", group = g, n = length(idx),
      truth = mean(el$truth), estimate = mean(el$estimate), bias = mean(el$bias), mae = mean(abs(el$bias)),
      rmse = sqrt(mean(el$bias^2)), coverage = mean(el$covered), interval_width = mean(el$upper - el$lower),
      as.list(trace_diagnostics(trace)))
  }
  standardised_range <- land$range / mean(xy_sd)
  blocks <- list(intercept = list(js$B0_output, land$B0 + land$B * mean(tr$environment)),
                 environment_slope = list(js$B_output, land$B * env_sd),
                 range = list(matrix(fit$infos$l_s_grid[js$idx_ls_output], ni, nc), standardised_range),
                 spatial_sd = list(js$sigmabs_output, land$field_sd))
  if (arm != "binary") blocks <- c(blocks, list(
    collection_coefficient = list(fit$results_output$beta_theta_output, input$detection$beta_theta),
    theta0 = list(fit$results_output$theta0_output, input$detection$theta0),
    p = list(fit$results_output$p_output, input$detection$p),
    q = list(fit$results_output$q_output, input$detection$q)))
  for (nm in names(blocks)) {
    b <- score_draw_block(blocks[[nm]][[1]], blocks[[nm]][[2]], nm)
    elements <- rbind(elements, b$elements); groups[[length(groups) + 1L]] <- b$summary
  }
  species <- do.call(rbind, lapply(seq_len(S), function(s) {
    idx <- ((s - 1L) * n + 1L):(s * n); el <- elements[elements$metric == "occupancy", ][idx, ]
    data.frame(species = land$species[s], target = land$prevalence[s], occupied = sum(tr$z[, s]),
               detections = sum(input$surveys[[arrangement]][[arm]]$OTU[, s]),
               truth = mean(el$truth), estimate = mean(el$estimate), bias = mean(el$bias),
               mae = mean(abs(el$bias)), rmse = sqrt(mean(el$bias^2)), coverage = mean(el$covered),
               as.list(trace_diagnostics(apply(x[idx, , , drop = FALSE], c(2, 3), mean))))
  }))
  field_draws <- reconstruct_field_draws(fit)
  field_median <- apply(field_draws, c(1, 2), median)
  tc <- sweep(tr$field, 2, colMeans(tr$field), "-"); fc <- sweep(field_median, 2, colMeans(field_median), "-")
  field <- do.call(rbind, lapply(seq_len(S), function(s) data.frame(species = land$species[s],
    target = land$prevalence[s], centred_rmse = sqrt(mean((fc[, s] - tc[, s])^2)),
    centred_correlation = suppressWarnings(cor(fc[, s], tc[, s])), centred_slope = sum(fc[, s] * tc[, s]) / sum(tc[, s]^2),
    zero_field_rmse = sqrt(mean(tc[, s]^2)), median_field_rms = sqrt(mean(field_median[, s]^2)))))
  grid <- fit$infos$l_s_grid; nearest <- which.min(abs(grid - standardised_range))
  range_table <- do.call(rbind, lapply(seq_len(nc), function(ch) data.frame(chain = ch, range = grid,
    frequency = tabulate(js$idx_ls_output[, ch], nbins = length(grid)) / ni)))
  range_summary <- data.frame(standardised_truth = standardised_range, nearest_grid = grid[nearest],
    mass_within_one_step = mean(abs(js$idx_ls_output - nearest) <= 1L),
    mass_at_boundary = mean(js$idx_ls_output %in% c(1L, length(grid))),
    posterior_mean = mean(grid[js$idx_ls_output]))
  amp <- quantile(js$sigmabs_output, c(.025, .5, .975), names = FALSE)
  amplitude <- data.frame(truth = land$field_sd, lower = amp[1], median = amp[2], upper = amp[3])
  lattice_index <- land$index$lattice
  lat <- reconstruct_lattice_draws(fit, land$environment[lattice_index], land$points[lattice_index, ])
  lattice_truth <- land$psi[lattice_index, ]
  bins <- distance_bins(land$points[lattice_index, ], tr$xy)
  lattice <- do.call(rbind, lapply(c("with", "without"), function(term) {
    est <- lat[[term]]; err <- est - lattice_truth
    rbind(data.frame(spatial_term = term, bin = "all", n = length(err), bias = mean(err), mae = mean(abs(err)),
                     rmse = sqrt(mean(err^2))),
          do.call(rbind, lapply(levels(bins$bin), function(b) { rows <- bins$bin == b
            data.frame(spatial_term = term, bin = b, n = sum(rows) * S, bias = mean(err[rows, ]),
                       mae = mean(abs(err[rows, ])), rmse = sqrt(mean(err[rows, ]^2))) })))
  }))
  set.seed(1L); subset <- sort(sample(length(lattice_index), lattice_subset))
  native <- suppressMessages(occJSDM::predictNewSites(fit,
    X_psi = data.frame(environment = land$environment[lattice_index][subset]),
    X_s = data.frame(longitude = land$points[lattice_index[subset], 1], latitude = land$points[lattice_index[subset], 2]),
    useEnvCov = TRUE, useSpatial = TRUE, useBiotic = FALSE, confidence = .95, verbose = FALSE))
  full <- reconstruct_lattice_draws(fit, land$environment[lattice_index][subset], land$points[lattice_index[subset], ], thin = 1L)
  # Native output is quantile x site x species; compare its median with the median of the unthinned draws.
  draws <- array(NA_real_, c(length(subset), S, ni * nc)); k <- 0L
  for (ch in seq_len(nc)) for (it in seq_len(ni)) { k <- k + 1L
    draws[, , k] <- plogis(sweep(as.matrix(full$X_new) %*% matrix(js$B_output[, , it, ch], ncol(full$X_new), S), 2, js$B0_output[, it, ch], "+") +
      full$Ks[, , js$idx_ls_output[it, ch]] %*% matrix(js$Bs_output[, , it, ch], fit$infos$ps, S)) }
  lattice_difference <- max(abs(apply(draws, c(1, 2), median) - native[2, , ]))
  list(groups = do.call(rbind, groups), elements = elements, species = species, field = field,
       range = range_table, range_summary = range_summary, amplitude = amplitude, lattice = lattice,
       probability = estimate, field_median = field_median, lattice_mean = lat$with,
       lattice_mean_nospatial = lat$without, lattice_distance = bins$distance,
       verification = list(basis = basis_difference, saved_probability = saved_difference,
                           lattice = lattice_difference, thin = 4L))
}
