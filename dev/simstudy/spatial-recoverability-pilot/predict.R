# Lattice prediction by Gaussian-process conditioning and the spatially uninformed reference.
# Requires kernel.R, generator.R and oracle.R to be sourced first.
uninformed_reference <- function(eta, sd = 1, gh = gauss_hermite(40L)) {
  vapply(eta, function(e) sum(gh$weights * plogis(e + sqrt(2) * sd * gh$nodes)) / sqrt(pi), numeric(1))
}

# `oracle` is list(species = integer vector, oracles = list of run_species_oracle results)
# for truth 2, or list(joint = run_joint_oracle result) for truth 1.
predict_lattice <- function(cm, arrangement, oracle, seed, tolerance = 1e-5) {
  idx <- cm$index[[arrangement]]; xy <- cm$points[idx, , drop = FALSE]; lat <- cm$points[cm$index$lattice, , drop = FALSE]
  offset_lat <- training_offset(cm, cm$index$lattice); n <- length(idx); S <- nrow(cm$species)
  set.seed(seed)
  records <- list(); cache <- list()
  conditional_for <- function(range, what) {
    key <- as.character(range)
    if (!is.null(cache[[key]])) return(cache[[key]])
    cf <- lattice_conditional(xy, lat, range)
    if (cf$error > tolerance) stop("Lattice conditional reconstruction error ", cf$error, " exceeds ", tolerance, " for ", what)
    records[[length(records) + 1L]] <<- data.frame(what = what, range = range, route = cf$route, rank = cf$rank,
                                                   min_eigen = cf$min_eigen, error = cf$error, stringsAsFactors = FALSE)
    cache[[key]] <<- cf
    cf
  }
  draw_lattice <- function(cf, site_block) {
    nd <- ncol(site_block)
    cf$A %*% site_block + cf$lower %*% matrix(rnorm(ncol(cf$lower) * nd), ncol(cf$lower), nd)
  }
  if (cm$truth == "t1") {
    o <- oracle$joint; cf <- conditional_for(cm$range, paste0("conditional-", arrangement, "-joint"))
    blocks <- lapply(o$site_draws, function(b) draw_lattice(cf, b))
    w <- joint_weights(cm); c_indep <- sqrt(1 - cm$shared_fraction)
    species_index <- seq_len(S)
    species_lattice <- function(j) w[j, 1] * blocks[[1]] + w[j, 2] * blocks[[2]] + c_indep * blocks[[2L + j]]
    species_site <- function(j) w[j, 1] * o$site_draws[[1]] + w[j, 2] * o$site_draws[[2]] + c_indep * o$site_draws[[2L + j]]
    offsets_site <- o$offset
  } else {
    species_index <- oracle$species
    species_lattice <- function(j) {
      r <- cm$ranges[species_index[j]]
      draw_lattice(conditional_for(r, paste0("conditional-", arrangement, "-r", r)), oracle$oracles[[j]]$site_draws[[1]])
    }
    species_site <- function(j) oracle$oracles[[j]]$site_draws[[1]]
    offsets_site <- do.call(cbind, lapply(oracle$oracles, function(o) o$offset))
  }
  m <- length(species_index)
  probability_mean <- matrix(NA_real_, nrow(lat), m); field_mean <- probability_mean
  site_field_mean <- matrix(NA_real_, n, m); site_probability_mean <- site_field_mean
  for (j in seq_len(m)) {
    s <- species_index[j]
    H <- species_lattice(j)
    probability_mean[, j] <- rowMeans(plogis(offset_lat[, s] + H)); field_mean[, j] <- rowMeans(H)
    Hs <- species_site(j)
    site_field_mean[, j] <- rowMeans(Hs); site_probability_mean[, j] <- rowMeans(plogis(offsets_site[, j] + Hs))
  }
  nm <- cm$species$species[species_index]
  colnames(probability_mean) <- colnames(field_mean) <- colnames(site_field_mean) <- colnames(site_probability_mean) <- nm
  list(species = species_index, probability_mean = probability_mean, field_mean = field_mean,
       site_field_mean = site_field_mean, site_probability_mean = site_probability_mean,
       conditional_records = do.call(rbind, records))
}

# Geometric descriptors of a training layout relative to the lattice (spec section 4).
layout_geometry <- function(cm, arrangement) {
  idx <- cm$index[[arrangement]]; xy <- cm$points[idx, , drop = FALSE]; lat <- cm$points[cm$index$lattice, , drop = FALSE]
  d2 <- outer(lat[, 1], xy[, 1], "-")^2 + outer(lat[, 2], xy[, 2], "-")^2
  nearest <- sqrt(apply(d2, 1, min))
  hc <- half_correlation(c(0.03, 0.08, 0.15))
  data.frame(arrangement = arrangement, median_nearest_training = median(nearest),
             within_half_short = mean(nearest <= hc[1]), within_half_medium = mean(nearest <= hc[2]),
             within_half_long = mean(nearest <= hc[3]),
             centre_separation = if (arrangement == "clustered") cm$arrangements$separation else NA_real_,
             stringsAsFactors = FALSE)
}
