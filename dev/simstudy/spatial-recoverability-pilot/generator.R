# Generator for the spatial recoverability pilot. Standalone; no occJSDM code.
# Requires kernel.R and seeds.R to be sourced first.
pilot_species <- function() {
  data.frame(species = sprintf("species%02d", 1:12),
             occupancy = rep(c(0.05, 0.25, 0.50, 0.75), each = 3),
             class = rep(c("short", "medium", "long"), 4),
             log_range = rep(log(c(0.03, 0.08, 0.15)), 4),
             slope = rep(c(0.4, -0.4), 6), stringsAsFactors = FALSE)
}

pilot_lattice <- function(m = 40L) {
  g <- (seq_len(m) - 0.5) / m
  xy <- as.matrix(expand.grid(x = g, y = g)); dimnames(xy) <- list(NULL, c("x", "y")); xy
}

pilot_cells <- function() {
  t1 <- expand.grid(range = c(0.03, 0.08, 0.15), n = c(100L, 300L, 900L))
  cells <- rbind(data.frame(truth = "t1", range = t1$range, n = t1$n),
                 data.frame(truth = "t2", range = NA_real_, n = c(100L, 300L, 900L)))
  cells$cell <- ifelse(cells$truth == "t1", sprintf("t1-r%.2f-n%d", cells$range, cells$n), sprintf("t2-n%d", cells$n))
  cells$truth <- as.character(cells$truth)
  cells
}

cluster_separation <- function(k) 0.2 * sqrt(10 / k)

axis_ratio_ok <- function(xy, max_ratio = 0.1) abs(sd(xy[, 1]) / sd(xy[, 2]) - 1) <= max_ratio

draw_separated_centres <- function(k, min_sep, max_tries = 200000L) {
  centres <- matrix(numeric(0), 0, 2)
  for (i in seq_len(max_tries)) {
    cand <- runif(2)
    ok <- nrow(centres) == 0L || min(sqrt(colSums((t(centres) - cand)^2))) >= min_sep
    if (ok) centres <- rbind(centres, cand)
    if (nrow(centres) == k) return(unname(centres))
  }
  stop("Could not place ", k, " cluster centres at separation ", min_sep)
}

make_arrangements <- function(n, seed, cluster_radius = 0.02, max_ratio = 0.1) {
  stopifnot(n %% 10L == 0L)
  set.seed(seed); rejections <- 0L
  repeat { uniform <- matrix(runif(2 * n), n, 2); if (axis_ratio_ok(uniform, max_ratio)) break; rejections <- rejections + 1L }
  k <- n %/% 10L; separation <- cluster_separation(k)
  repeat {
    centres <- draw_separated_centres(k, separation)
    # A point falling outside the unit square is redrawn within its cluster rather than
    # clipped: clipping stacks points on a corner and produces duplicate sites.
    clustered <- do.call(rbind, lapply(seq_len(k), function(j) {
      pts <- matrix(NA_real_, 10, 2)
      for (i in 1:10) repeat {
        r <- cluster_radius * sqrt(runif(1)); a <- runif(1, 0, 2 * pi)
        p <- centres[j, ] + r * c(cos(a), sin(a))
        if (all(p >= 0 & p <= 1)) { pts[i, ] <- p; break }
      }
      pts
    }))
    if (axis_ratio_ok(clustered, max_ratio)) break
    rejections <- rejections + 1L
  }
  name <- function(m) { m <- unname(m); dimnames(m) <- list(NULL, c("x", "y")); m }
  list(uniform = name(uniform), clustered = name(clustered), centres = centres,
       separation = separation, rejections = rejections)
}

solve_intercepts <- function(offset_lattice, targets) {
  vapply(seq_along(targets), function(s) uniroot(function(b) mean(plogis(b + offset_lattice[, s])) - targets[s],
                                                 c(-30, 30), tol = 1e-12)$root, numeric(1))
}

make_community <- function(truth, n, replicate, range = NA_real_, env_range = 0.5, env_noise_sd = 0.3,
                           shared_fraction = 0.9) {
  stopifnot(truth %in% c("t1", "t2"), truth == "t2" || is.finite(range))
  cell <- if (truth == "t1") sprintf("t1-r%.2f-n%d", range, n) else sprintf("t2-n%d", n)
  sp <- pilot_species(); S <- nrow(sp)
  seeds <- c(sites = pilot_seed(pilot_label("sites", cell, paste0("rep", replicate))),
             data = pilot_seed(pilot_label("data", cell, paste0("rep", replicate))))
  arr <- make_arrangements(n, seeds[["sites"]])
  lattice <- pilot_lattice()
  points <- rbind(arr$uniform, arr$clustered, lattice)
  index <- list(uniform = seq_len(n), clustered = n + seq_len(n), lattice = 2L * n + seq_len(nrow(lattice)))
  stopifnot(!anyDuplicated(points))
  N <- nrow(points)
  set.seed(seeds[["data"]])
  conditioning <- list()
  factor_for <- function(r, label) {
    f <- stable_factor(sq_exp(points, points, r))
    conditioning[[length(conditioning) + 1L]] <<- data.frame(what = label, range = r, route = f$route,
                                                             condition = f$condition, rank = f$rank, error = f$error)
    f$lower
  }
  environment <- as.vector(factor_for(env_range, "environment") %*% rnorm(N)) + rnorm(N, 0, env_noise_sd)
  if (truth == "t1") {
    Lf <- factor_for(range, "shared-and-independent")
    shared_fields <- Lf %*% matrix(rnorm(N * 2), N, 2)
    independent <- Lf %*% matrix(rnorm(N * S), N, S)
    loadings <- matrix(rnorm(2 * S), S, 2)
    loadings <- loadings / sqrt(rowSums(loadings^2))
    shared_part <- sqrt(shared_fraction) * shared_fields %*% t(loadings)
    spatial <- shared_part + sqrt(1 - shared_fraction) * independent
    ranges <- rep(range, S)
  } else {
    ranges <- exp(sp$log_range)
    distinct <- sort(unique(ranges))
    factors <- lapply(distinct, function(r) factor_for(r, paste0("class-", r)))
    names(factors) <- as.character(distinct)
    spatial <- vapply(seq_len(S), function(s) as.vector(factors[[as.character(ranges[s])]] %*% rnorm(N)), numeric(N))
    shared_fields <- NULL; independent <- NULL; loadings <- NULL; shared_part <- NULL
  }
  offset <- environment %o% sp$slope + spatial
  B0 <- solve_intercepts(offset[index$lattice, , drop = FALSE], sp$occupancy)
  psi <- plogis(sweep(offset, 2, B0, "+"))
  z <- matrix(as.integer(runif(N * S) < psi), N, S)
  colnames(spatial) <- colnames(psi) <- colnames(z) <- sp$species
  list(truth = truth, cell = cell, n = n, replicate = replicate, range = range, species = sp,
       points = points, index = index, arrangements = arr, environment = environment,
       spatial = spatial, shared_fields = shared_fields, independent = independent, loadings = loadings,
       shared_part = shared_part, shared_fraction = shared_fraction, ranges = ranges, spatial_sd = rep(1, S),
       psi = psi, z = z, B0 = B0, B = sp$slope, env_range = env_range, env_noise_sd = env_noise_sd,
       conditioning = do.call(rbind, conditioning), seeds = seeds)
}
