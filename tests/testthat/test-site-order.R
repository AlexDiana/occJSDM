# Site order. Every per-site matrix of a fit must follow one order, the order
# of the sites in fit$infos$data_info:
#
#   - data with one row per site (binary and continuous fits) keep the rows of
#     data$info as given, and their sites are named by the Site column;
#   - occupancy and two-stage data are sorted by Site, then Sample, then
#     Primer, and their sites are named in fit$infos$siteNames.
#
# Until 0.3.1 the covariates of a one-row-per-site fit were sorted by Site
# while the response stayed in input order, so any data whose Site column was
# not already sorted (FHT1, FHT2, ..., FHT10 is file order, not sorted order)
# was fitted with covariates attached to the wrong sites, silently. These
# tests plant site labels that are not in sorted order, as character, factor
# and integer, and check each matrix against the input rows by identity.

SITE_ORDER_MCMC <- list(nchain = 1, nburn = 10, niter = 10, nthin = 1)

# Labels for n sites, in input order, that sorting would rearrange.
unsorted_sites <- function(kind, n) {
  switch(kind,
    character = paste0("FHT", seq_len(n)),
    factor    = factor(paste0("FHT", seq_len(n))),
    integer   = c(seq(2L, n, by = 2L), seq(1L, n, by = 2L)),
    stop("unknown kind"))
}
site_kinds <- c("character", "factor", "integer")

# One row per site: a numeric and a categorical occupancy covariate, two
# coordinates, and either a binary or a continuous response driven by x.
single_visit_data <- function(n = 40L, S = 4L, response = "binary", seed = 11L) {
  set.seed(seed)
  x <- stats::rnorm(n)
  habitat <- sample(c("forest", "open"), n, replace = TRUE)
  eta <- -0.5 + 1.5 * x
  Y <- if (response == "binary") {
    sapply(seq_len(S), function(s) stats::rbinom(n, 1, stats::plogis(eta)))
  } else {
    sapply(seq_len(S), function(s) eta + stats::rnorm(n, sd = 0.5))
  }
  colnames(Y) <- paste0("sp", seq_len(S))
  info <- data.frame(x = x, habitat = habitat,
                     east = stats::runif(n), north = stats::runif(n))
  list(info = info, OTU = Y)
}

fit_site_order <- function(data, occCovariates = c("x", "habitat"),
                           spatCovariates = NULL, collCovariates = NULL,
                           listParams = list(n_factors = 0, n_lattrait = 0),
                           MCMCparams = SITE_ORDER_MCMC, seed = 5L) {
  set.seed(seed)
  suppressMessages(suppressWarnings(runOccJSDM(
    data, occCovariates = occCovariates, spatCovariates = spatCovariates,
    collCovariates = collCovariates, listParams = listParams,
    MCMCparams = MCMCparams)))
}

standardise <- function(v) (v - mean(v)) / stats::sd(v)

# --- one row per site ------------------------------------------------------

test_that("binary covariates stay on their input rows when Site is unsorted", {
  data <- single_visit_data()
  n <- nrow(data$info)
  for (kind in site_kinds) {
    data$info$Site <- unsorted_sites(kind, n)
    fit <- fit_site_order(data)
    expect_equal(fit$infos$model, "binary", info = kind)
    expect_equal(unname(fit$X_psi[, "x"]), standardise(data$info$x), info = kind)
    expect_equal(unname(fit$X_psi[, "habitatopen"]),
                 as.numeric(data$info$habitat == "open"), info = kind)
    expect_equal(unname(fit$infos$OTU), unname(data$OTU), info = kind)
    expect_equal(as.character(fit$infos$data_info$Site),
                 as.character(data$info$Site), info = kind)
  }
})

test_that("continuous covariates stay on their input rows when Site is unsorted", {
  data <- single_visit_data(response = "continuous")
  n <- nrow(data$info)
  for (kind in site_kinds) {
    data$info$Site <- unsorted_sites(kind, n)
    fit <- fit_site_order(data)
    expect_equal(fit$infos$model, "continuous", info = kind)
    expect_equal(unname(fit$X_psi[, "x"]), standardise(data$info$x), info = kind)
    expect_equal(unname(fit$X_psi[, "habitatopen"]),
                 as.numeric(data$info$habitat == "open"), info = kind)
  }
})

test_that("spatial coordinates stay on their input rows when Site is unsorted", {
  for (response in c("binary", "continuous")) {
    data <- single_visit_data(response = response)
    n <- nrow(data$info)
    for (kind in site_kinds) {
      data$info$Site <- unsorted_sites(kind, n)
      fit <- fit_site_order(data, spatCovariates = c("east", "north"),
                            listParams = list(n_factors = 0, n_lattrait = 0,
                                              n_supportpoints = 8))
      label <- paste(response, kind)
      expect_equal(unname(fit$Xs[, 1]), standardise(data$info$east), info = label)
      expect_equal(unname(fit$Xs[, 2]), standardise(data$info$north), info = label)
      expect_equal(unname(fit$X_psi[, "x"]), standardise(data$info$x), info = label)
    }
  }
})

test_that("site labels alone do not change a one-row-per-site fit", {
  # The same rows under sorted and unsorted labels are the same data, so a
  # seeded fit must give the same draws, and site-level summaries the same
  # numbers in the same rows.
  for (response in c("binary", "continuous")) {
    data <- single_visit_data(response = response)
    n <- nrow(data$info)
    data$info$Site <- seq_len(n)
    reference <- fit_site_order(data, listParams = list(n_factors = 1, n_lattrait = 0))
    for (kind in site_kinds) {
      data$info$Site <- unsorted_sites(kind, n)
      fit <- fit_site_order(data, listParams = list(n_factors = 1, n_lattrait = 0))
      label <- paste(response, kind)
      expect_identical(fit$results_output, reference$results_output, info = label)
      expect_identical(fit$X_psi, reference$X_psi, info = label)
      # The scores are the same numbers; only their site labels differ.
      expect_identical(unname(returnOrdinationScores(fit)),
                       unname(returnOrdinationScores(reference)), info = label)
    }
  }
})

test_that("a one-row-per-site fit names its sites by Site", {
  # Until 0.3.1 these fits stored 1:n in infos$siteNames whatever Site held,
  # so tables and plots labelled their sites by row number.
  data <- single_visit_data()
  n <- nrow(data$info)
  for (kind in site_kinds) {
    sites <- unsorted_sites(kind, n)
    data$info$Site <- sites
    fit <- fit_site_order(data, listParams = list(n_factors = 1, n_lattrait = 0))
    expect_identical(fit$infos$siteNames, sites, info = kind)
    expect_identical(dimnames(returnOrdinationScores(fit))[[2]], as.character(sites),
                     info = kind)
  }
  # With no Site column there are no names to store: the row number stands in.
  data$info$Site <- NULL
  fit <- fit_site_order(data)
  expect_identical(fit$infos$siteNames, seq_len(n))
})

test_that("a planted slope is recovered whatever the site labels", {
  skip_on_cran()
  set.seed(1)
  n <- 150L
  x <- stats::rnorm(n)
  truth <- stats::plogis(-0.5 + 1.5 * x)
  Y <- sapply(1:6, function(s) stats::rbinom(n, 1, truth))
  colnames(Y) <- paste0("sp", 1:6)
  for (kind in site_kinds) {
    data <- list(info = data.frame(Site = unsorted_sites(kind, n), x = x), OTU = Y)
    fit <- fit_site_order(data, occCovariates = "x",
                          MCMCparams = list(nchain = 1, nburn = 300, niter = 300, nthin = 1))
    # The fitted slope is on the standardised covariate; sd(x) is close to 1.
    # Misaligned covariates give a slope near 0 (measured: -0.04 to 0.15).
    slope <- mean(fit$results_output$jsdm_output$B_output)
    expect_gt(slope, 0.8)
    expect_lt(slope, 2)
    # Prediction at the fitted sites' own covariates is close to the truth
    # (measured mean absolute error: 0.03 aligned, 0.22 misaligned).
    predicted <- suppressMessages(predictNewSites(fit, X_psi = data.frame(x = x),
                                                  verbose = FALSE))
    expect_lt(mean(abs(predicted[2, , ] - truth)), 0.1)
  }
})

test_that("site WAIC scores each site against its own covariates", {
  data <- single_visit_data()
  n <- nrow(data$info)
  data$info$Site <- seq_len(n)
  reference <- suppressWarnings(computeSiteWAIC(fit_site_order(data)))
  for (kind in site_kinds) {
    data$info$Site <- unsorted_sites(kind, n)
    score <- suppressWarnings(computeSiteWAIC(fit_site_order(data)))
    expect_equal(score$pointwise$site, as.character(data$info$Site), info = kind)
    expect_equal(score$pointwise$waic, reference$pointwise$waic, info = kind)
  }
})

# --- repeated visits -------------------------------------------------------

# Occupancy (one row per sample) or two-stage (two PCR rows per sample) data
# with site-level and sample-level covariates and coordinates that identify
# each site and sample, so that every row of the fit can be traced to its
# input by name. Rows are shuffled so no ordering survives by accident.
repeated_visit_data <- function(kind, model, n = 40L, M = 2L, K = 2L, S = 3L,
                                seed = 21L) {
  set.seed(seed)
  sites <- unsorted_sites(kind, n)
  site_table <- data.frame(Site = sites, x = stats::rnorm(n),
                           east = stats::runif(n), north = stats::runif(n))
  samples <- data.frame(site_row = rep(seq_len(n), each = M),
                        Sample = rep(seq_len(M), n))
  samples$effort <- stats::rnorm(nrow(samples))
  rows <- if (model == "two_stage") rep(seq_len(nrow(samples)), each = K) else
    seq_len(nrow(samples))
  info <- cbind(site_table[samples$site_row[rows], ],
                samples[rows, c("Sample", "effort")])
  if (model == "two_stage") info$Primer <- "A"
  otu <- matrix(stats::rbinom(nrow(info) * S, 1, 0.4), nrow(info), S,
                dimnames = list(NULL, paste0("sp", seq_len(S))))
  # A row identifier carried as an unused column and as the first species'
  # reads would both be changed by a fit; keep the key outside the data.
  key <- paste(info$Site, info$Sample, stats::ave(seq_len(nrow(info)),
               info$Site, info$Sample, FUN = seq_along))
  perm <- sample(nrow(info))
  list(data = list(info = info[perm, ], OTU = otu[perm, , drop = FALSE]),
       site_table = site_table, key = key[perm])
}

test_that("occupancy and two-stage matrices follow the sorted site order", {
  for (model in c("occupancy", "two_stage")) {
    for (kind in site_kinds) {
      built <- repeated_visit_data(kind, model)
      data <- built$data
      truth <- built$site_table
      fit <- fit_site_order(data, occCovariates = "x", collCovariates = "effort",
                            spatCovariates = c("east", "north"),
                            listParams = list(n_factors = 0, n_lattrait = 0,
                                              n_supportpoints = 8))
      label <- paste(model, kind)
      expect_equal(fit$infos$model, model, info = label)

      # Sites are in sorted order and named.
      site_names <- as.character(fit$infos$siteNames)
      at <- match(site_names, as.character(truth$Site))
      expect_false(anyNA(at), info = label)
      expect_equal(site_names, as.character(sort(truth$Site)), info = label)
      expect_equal(names(fit$infos$M), site_names, info = label)

      # Site-level matrices: row i belongs to site_names[i].
      expect_equal(unname(fit$X_psi[, "x"]), standardise(truth$x)[at], info = label)
      expect_equal(unname(fit$Xs[, 1]), standardise(truth$east)[at], info = label)
      expect_equal(unname(fit$Xs[, 2]), standardise(truth$north)[at], info = label)
      expect_equal(rownames(computePredictiveOccupancyProbs(fit)), site_names, info = label)
      expect_equal(rownames(computeConditionalOccupancyProbs(fit)), site_names, info = label)

      # Observation rows: the stored response is the input response, row for
      # row, under the stored data_info.
      stored <- fit$infos$data_info
      stored_key <- paste(stored$Site, stored$Sample,
                          stats::ave(seq_len(nrow(stored)), stored$Site, stored$Sample,
                                     FUN = seq_along))
      # Rows of one (Site, Sample) are exchangeable replicates, so compare
      # their detections as sorted sets.
      input_by_sample <- split(data$OTU[, 1], paste(data$info$Site, data$info$Sample))
      stored_by_sample <- split(fit$infos$OTU[, 1], paste(stored$Site, stored$Sample))
      expect_equal(lapply(stored_by_sample, sort),
                   lapply(input_by_sample[names(stored_by_sample)], sort), info = label)
      expect_equal(length(stored_key), nrow(data$info), info = label)

      # The site index of every sample points at that sample's own site.
      sample_rows <- stored[!duplicated(stored[, c("Site", "Sample")]), ]
      expect_equal(site_names[fit$infos$list_idx$idx_z_w],
                   as.character(sample_rows$Site), info = label)
      effort_by_sample <- data$info$effort[match(
        paste(sample_rows$Site, sample_rows$Sample),
        paste(data$info$Site, data$info$Sample))]
      expect_equal(unname(fit$X_theta[, "effort"]), standardise(effort_by_sample),
                   info = label)
    }
  }
})

test_that("returnLatentPresences() rows carry their own site's probabilities", {
  built <- repeated_visit_data("character", "two_stage")
  fit <- fit_site_order(built$data, occCovariates = "x")
  table <- returnLatentPresences(fit, idx_species = 1)
  occupancy <- computeConditionalOccupancyProbs(fit)
  expect_equal(table$CondOccProb, unname(occupancy[as.character(table$Site), 1]))
  expect_equal(table$x, built$site_table$x[match(table$Site, built$site_table$Site)])
})

test_that("an occupancy slope is recovered under unsorted site names", {
  skip_on_cran()
  set.seed(31)
  n <- 150L; M <- 3L
  x <- stats::rnorm(n)
  z <- sapply(1:6, function(s) stats::rbinom(n, 1, stats::plogis(-0.5 + 1.5 * x)))
  site <- rep(paste0("FHT", seq_len(n)), each = M)
  Y <- z[rep(seq_len(n), each = M), ] * stats::rbinom(n * M * 6, 1, 0.7)
  colnames(Y) <- paste0("sp", 1:6)
  info <- data.frame(Site = site, Sample = rep(seq_len(M), n), x = rep(x, each = M))
  perm <- sample(nrow(info))
  fit <- fit_site_order(list(info = info[perm, ], OTU = Y[perm, ]), occCovariates = "x",
                        MCMCparams = list(nchain = 1, nburn = 300, niter = 300, nthin = 1))
  expect_equal(fit$infos$model, "occupancy")
  expect_gt(mean(fit$results_output$jsdm_output$B_output), 0.8)
})

# --- traits ----------------------------------------------------------------

test_that("traits are matched to species by name, not by row", {
  data <- single_visit_data()
  species <- colnames(data$OTU)
  traits <- data.frame(size = c(4, 1, 3, 2), row.names = species)
  data$traits <- traits[c(3, 1, 4, 2), , drop = FALSE]
  fit <- fit_site_order(data)
  expect_equal(rownames(fit$Tr), species)
  expect_equal(unname(fit$Tr[, "size"]), standardise(traits$size))
})
