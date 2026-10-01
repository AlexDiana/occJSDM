# Generator for the Lesson 2 site-arrangement sweep. Standalone: no occJSDM
# simulation or fitting internals. See PLAN.md for the frozen design.

sweep_seed <- function(label, replicate) {
  stopifnot(is.character(label), length(label) == 1L, nzchar(label),
            length(replicate) == 1L, replicate >= 1L, replicate == as.integer(replicate))
  as.integer(sum(as.integer(charToRaw(label))) * 1000L + replicate)
}

sq_exp <- function(a, b, range) {
  exp(-(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2) / (2 * range^2))
}

make_lattice <- function(m = 40L) {
  g <- (seq_len(m) - 0.5) / m
  xy <- as.matrix(expand.grid(x = g, y = g))
  dimnames(xy) <- list(NULL, c("x", "y"))
  xy
}

axis_sd_ratio <- function(xy) sd(xy[, 1]) / sd(xy[, 2])

draw_separated_centres <- function(k, min_sep, max_tries = 100000L) {
  centres <- matrix(numeric(0), 0, 2)
  for (i in seq_len(max_tries)) {
    cand <- runif(2)
    ok <- nrow(centres) == 0L || min(sqrt(colSums((t(centres) - cand)^2))) >= min_sep
    if (ok) centres <- rbind(centres, cand)
    if (nrow(centres) == k) return(unname(centres))
  }
  stop("Could not place ", k, " cluster centres at separation ", min_sep)
}

make_arrangements <- function(seed, n = 100L, pair_distance = 0.01, cluster_radius = 0.02,
                              centre_separation = 0.2, max_ratio = 0.1) {
  set.seed(seed)
  rejections <- 0L
  repeat {
    spread <- matrix(runif(2 * n), n, 2)
    if (abs(axis_sd_ratio(spread) - 1) <= max_ratio) break
    rejections <- rejections + 1L
  }
  base <- spread[1:80, , drop = FALSE]
  partner_of <- sample(80L, 20L)
  partner <- matrix(NA_real_, 20L, 2L)
  for (i in 1:20) {
    repeat {
      angle <- runif(1, 0, 2 * pi)
      p <- base[partner_of[i], ] + pair_distance * c(cos(angle), sin(angle))
      if (all(p >= 0 & p <= 1)) break
    }
    partner[i, ] <- p
  }
  pairs <- rbind(base, partner)
  repeat {
    centres <- draw_separated_centres(10L, centre_separation)
    clustered <- do.call(rbind, lapply(1:10, function(k) {
      r <- cluster_radius * sqrt(runif(10)); a <- runif(10, 0, 2 * pi)
      pts <- cbind(centres[k, 1] + r * cos(a), centres[k, 2] + r * sin(a))
      pmin(pmax(pts, 0), 1)
    }))
    if (abs(axis_sd_ratio(clustered) - 1) <= max_ratio) break
    rejections <- rejections + 1L
  }
  g <- (1:10 - 0.5) / 10
  grid <- as.matrix(expand.grid(x = g, y = g))
  out <- lapply(list(spread = spread, pairs = pairs, clustered = clustered, grid = grid),
                function(m) { m <- unname(m); dimnames(m) <- list(NULL, c("x", "y")); m })
  attr(out, "cluster_centres") <- centres
  attr(out, "partner_of") <- partner_of
  attr(out, "rejections") <- rejections
  out
}

make_landscape <- function(replicate, range = 0.03, env_range = 0.5, env_noise_sd = 0.3,
                           prevalence = c(.05, .05, .25, .25, .25, .75, .75, .75),
                           slopes = rep(c(.4, -.4), 4L), jitter = 1e-8) {
  sites <- make_arrangements(sweep_seed("spatial-design-sites", replicate))
  lattice <- make_lattice()
  points <- rbind(sites$spread, sites$pairs[81:100, ], sites$clustered, sites$grid, lattice)
  index <- list(spread = 1:100, pairs = c(1:80, 101:120), clustered = 121:220,
                grid = 221:320, lattice = 321:1920)
  stopifnot(nrow(points) == 1920L, !anyDuplicated(points))
  set.seed(sweep_seed("spatial-design-data", replicate))
  S <- length(prevalence); N <- nrow(points)
  env_lower <- t(chol(sq_exp(points, points, env_range) + diag(jitter, N)))
  environment <- as.vector(env_lower %*% rnorm(N)) + rnorm(N, 0, env_noise_sd)
  field_lower <- t(chol(sq_exp(points, points, range) + diag(jitter, N)))
  field <- field_lower %*% matrix(rnorm(N * S), N, S)
  offset <- environment %o% slopes + field
  B0 <- vapply(seq_len(S), function(s) uniroot(function(b)
    mean(plogis(b + offset[index$lattice, s])) - prevalence[s], c(-30, 30), tol = 1e-12)$root, numeric(1))
  psi <- plogis(sweep(offset, 2, B0, "+"))
  z <- matrix(as.integer(runif(N * S) < psi), N, S)
  species <- sprintf("species%02d", seq_len(S))
  colnames(field) <- colnames(psi) <- colnames(z) <- species
  list(points = points, index = index, sites = sites, environment = environment,
       field = field, psi = psi, z = z, B0 = B0, B = slopes, prevalence = prevalence,
       range = range, field_sd = 1, env_range = env_range, env_noise_sd = env_noise_sd,
       species = species, cluster_centres = attr(sites, "cluster_centres"),
       rejections = attr(sites, "rejections"))
}

make_detection_parameters <- function(replicate, S = 8L, P = 2L) {
  set.seed(sweep_seed("spatial-design-detection", replicate))
  list(beta_theta = rbind(qlogis(runif(S, .2, .5)), rep(c(1, -1), length.out = S)),
       theta0 = runif(S, .02, .1), p = matrix(runif(P * S, .3, .6), P, S),
       q = matrix(.01 + .04 * runif(P * S), P, S))
}

make_survey <- function(land, arrangement, replicate, det, M = 2L, P = 2L, K = 6L) {
  idx <- land$index[[arrangement]]; n <- length(idx); S <- ncol(land$psi)
  xy <- land$points[idx, , drop = FALSE]; env <- land$environment[idx]
  z <- land$z[idx, , drop = FALSE]
  site <- data.frame(Site = seq_len(n), environment = env, longitude = xy[, 1], latitude = xy[, 2])
  binary <- list(info = site, OTU = z, traits = NULL)
  set.seed(sweep_seed(paste0("spatial-design-survey-", arrangement), replicate))
  sample_site <- rep(seq_len(n), each = M)
  collection <- rnorm(n * M)
  Xt <- cbind(1, as.numeric(scale(collection)))
  theta <- plogis(Xt %*% det$beta_theta)
  wprob <- ifelse(z[sample_site, ] == 1L, theta, matrix(det$theta0, n * M, S, byrow = TRUE))
  w <- matrix(as.integer(runif(n * M * S) < wprob), n * M, S)
  sample <- rep(seq_len(n * M), each = P * K)
  primer <- rep(rep(seq_len(P), each = K), times = n * M)
  row_site <- sample_site[sample]
  info <- data.frame(Site = row_site, Sample = sample, Primer = primer,
                     environment = env[row_site], longitude = xy[row_site, 1],
                     latitude = xy[row_site, 2], collection = collection[sample])
  prob <- ifelse(w[sample, ] == 1L, det$p[primer, ], det$q[primer, ])
  y <- matrix(as.integer(runif(length(prob)) < prob), nrow(prob), S,
              dimnames = list(NULL, land$species))
  list(binary = binary, two_stage = list(info = info, OTU = y, traits = NULL),
       truth = list(sites = idx, xy = xy, environment = env, z = z,
                    psi = land$psi[idx, , drop = FALSE], field = land$field[idx, , drop = FALSE],
                    w = w, Xt = Xt, theta = theta, collection = collection))
}

arrangement_statistics <- function(xy, range, environment_sites, environment_lattice) {
  d <- as.matrix(dist(xy)); diag(d) <- Inf
  nn <- apply(d, 1, min)
  threshold <- range * sqrt(2 * log(2))
  K <- sq_exp(xy, xy, range)
  sx <- range / sd(xy[, 1]); sy <- range / sd(xy[, 2])
  data.frame(mean_nearest_neighbour = mean(nn),
             fraction_with_half_neighbour = mean(nn < threshold),
             mean_neighbours_above_half = mean(rowSums(d < threshold)),
             effective_rank = sum(diag(K))^2 / sum(K^2),
             standardised_range_x = sx, standardised_range_y = sy,
             axis_sd_ratio = axis_sd_ratio(xy),
             environment_span = diff(range(environment_sites)) / diff(range(environment_lattice)),
             in_grid = sx >= 0.01 & sx <= 0.30 & sy >= 0.01 & sy <= 0.30)
}

make_sweep_input <- function(replicate) {
  land <- make_landscape(replicate)
  det <- make_detection_parameters(replicate)
  arrangements <- c("spread", "pairs", "clustered", "grid")
  surveys <- lapply(arrangements, function(a) make_survey(land, a, replicate, det))
  names(surveys) <- arrangements
  statistics <- do.call(rbind, lapply(arrangements, function(a) data.frame(
    replicate = replicate, arrangement = a,
    arrangement_statistics(land$points[land$index[[a]], ], land$range,
                           land$environment[land$index[[a]]],
                           land$environment[land$index$lattice]))))
  list(replicate = replicate, landscape = land, detection = det, surveys = surveys,
       statistics = statistics,
       settings = list(n = 100L, S = 8L, M = 2L, P = 2L, K = 6L, range = land$range,
                       seeds = c(sites = sweep_seed("spatial-design-sites", replicate),
                                 data = sweep_seed("spatial-design-data", replicate),
                                 detection = sweep_seed("spatial-design-detection", replicate)),
                       observation_model = "Direct Bernoulli detections; p and q are positive-observation probabilities"))
}
