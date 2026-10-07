# Lesson 2 Site-Arrangement Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the prespecified site-arrangement sweep in `dev/simstudy/spatial-design-sweep/PLAN.md` and turn its results into the worked part of Lesson 2.

**Architecture:** A standalone R generator builds one landscape per community and four 100-site arrangements over it; an oracle elliptical-slice sampler (sourced from the September diagnosis, unmodified) gives the perfect-observation ceiling per cell; a frozen copy of the package is fitted on every cell in binary and two-stage arms by detached workers; scoring, selection, audit and plotting scripts follow the September study conventions; an exporter writes a compact bundle that Lesson 2 renders from.

**Tech Stack:** R 4.5, occJSDM (frozen main revision), Rcpp/RcppArmadillo for the oracle, `posterior` for diagnostics, testthat for research tests, dplyr/tidyr/ggplot2/knitr/rmarkdown for the lesson.

**Spec:** `dev/simstudy/spatial-design-sweep/PLAN.md` (committed on branch `codex/lesson-2-spatial-sweep`, worktree `~/src/occJSDM-worktrees/lesson-2-spatial-sweep`). This plan lives beside it because `docs/` is gitignored in this repository.

## Global Constraints

- Work only in the worktree `~/src/occJSDM-worktrees/lesson-2-spatial-sweep` on branch `codex/lesson-2-spatial-sweep`. Never touch `~/src/occJSDM` (Doug's checkout) or anything under `~/Documents` (unreadable).
- Markdown written in this repo: one line per paragraph, no hard wrapping, no pipe tables in files edited in RStudio Visual mode (lesson Rmd and plan files), no em-dashes anywhere, escape `\>` `\<` `\~` `\|` in prose, never inside code. Commit messages: no em-dashes; end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.
- Every script refuses to overwrite a completed output. Raw outputs go under `dev/simstudy/results/spatial-design-20261001/` (gitignored). Compact outputs are committed under `dev/simstudy/spatial-design-sweep/results/`.
- Long jobs are launched with `nohup ... &` from Bash so they survive the session; the Bash tool's own background mode is not used for fits.
- The oracle sampler `dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp` is sourced, never copied or edited. The September helpers `dev/simstudy/spatial-targeted-recheck/score.R` and `dev/simstudy/spatial-amplitude-prior/robust.R` are sourced for `trace_diagnostics`, `independent_bases`, `reconstruct_spatial_draws`, `score_draw_block`, `verify_saved_probability`, `robust_trace_diagnostics` and `robust_flag_rows`; their MD5s are recorded.
- No package code, default or prior changes. Fit settings exactly as the spec: `n_factors = 0`, `n_lattrait = 0`, `n_supportpoints = 100`, `occCovariates = "environment"`, `spatCovariates = c("longitude", "latitude")`, `collCovariates = "collection"` for two-stage, defaults otherwise.
- Schedules: oracle 4 chains, 1,000 burn-in, 2,000 retained, doubled once on a flag; initial fits 2 chains, 3,000 burn-in, 5,000 retained; longer fits 4 chains, 6,000 burn-in, 12,000 retained. Selection by the spec's rule before any outcome comparison.
- Reading rules for the lesson vocabulary exactly as the spec: informative = binary-control centred field RMSE at least 20% below zero-field baseline and correlation at least 0.5 in all three communities; uninformative = reduction below 10% or correlation below 0.3 in all three; else intermediate; range recovered = at least half the posterior mass within one grid step of the generating value in all three.
- At most eight concurrent single-threaded fit workers; oracle with two workers.

---

## File structure

Created under `dev/simstudy/spatial-design-sweep/`:

- `generator.R`: landscape, arrangements, surveys, design statistics. No package internals.
- `oracle.R`: oracle covariance, sampler driver, oracle scoring, oracle lattice prediction.
- `score.R`: full-fit reconstruction (sites and lattice), scoring tables, verification differences.
- `run.R`: modes `prepare`, `oracle`, `freeze`, `pilot`, `initial`, `long`.
- `summarise.R`: modes `select` and `final`; aggregation, paired contrasts, reading labels.
- `verify.R`: independent audit of every selected fit and oracle result.
- `plot.R`: figures for the study report and the lesson.
- `export-teaching.R`: writes `vignettes/teaching-data/spatial-lesson.rds`.
- `verify-lesson.R`: checks the bundle against the committed compact results (archive mode always available).
- `test-generator.R`, `test-oracle.R`, `test-score.R`: research tests.
- `README.md`: reproduction commands and record.
- `results/`: committed compact CSVs and PNGs.

Modified: `vignettes/occJSDM-lesson-2.Rmd` (worked sections replace the outline), `vignettes/LESSON-PLAN.md`, `README.md` (Lesson 2 bullet), `TODO.md` (Lesson 2 item), `vignettes/occJSDM.Rmd` (Lesson 2 bullet), `dev/simstudy/spatial-design-sweep/PLAN.md` (amendments only).

---

### Task 1: Generator and design statistics

**Files:**
- Create: `dev/simstudy/spatial-design-sweep/generator.R`
- Create: `dev/simstudy/spatial-design-sweep/test-generator.R`
- Modify: `dev/simstudy/spatial-design-sweep/PLAN.md` (append amendment 1)

**Interfaces:**
- Produces: `sweep_seed(label, replicate)`; `sq_exp(a, b, range)`; `make_lattice(m = 40L)`; `make_arrangements(seed)` returning a named list `spread`, `pairs`, `clustered`, `grid` of 100 x 2 matrices with columns `x`, `y`, with attribute `cluster_centres` (10 x 2) and attribute `rejections` (integer); `make_landscape(replicate)`; `make_detection_parameters(replicate)`; `make_survey(land, arrangement, replicate, det)`; `arrangement_statistics(xy, range, environment_sites, environment_lattice)`; `make_sweep_input(replicate)` returning `list(replicate, landscape, detection, surveys, statistics, settings)` where `surveys[[arrangement]]` has `binary`, `two_stage` and `truth`.

- [ ] **Step 1: Append amendment 1 to the protocol**

In `dev/simstudy/spatial-design-sweep/PLAN.md`, replace the line `None yet. The production revision used for the full fits will be appended here before the first fit.` with:

```markdown
- **1 October 2026, amendment 1, before any simulation.** The fitter standardises each coordinate axis by its own standard deviation, so a field that is isotropic in raw coordinates is represented with a slightly anisotropic kernel whenever the two axes of an arrangement have unequal spreads. To keep this mismatch small and equal across designs, the arrangement generator rejects and redraws a spread design or a set of cluster centres whose 100 sites have axis standard deviations differing by more than 10 percent, counting the rejections. The grid is exact. The design table records the ratio of axis standard deviations for every arrangement, and the oracle uses the true isotropic kernel at the sites, so the oracle-versus-binary-control contrast includes this representational difference alongside parameter estimation.
- The production revision used for the full fits will be appended here by the freeze step before the first fit.
```

- [ ] **Step 2: Write the failing generator tests**

Create `dev/simstudy/spatial-design-sweep/test-generator.R`:

```r
# Run from the repository root: Rscript dev/simstudy/spatial-design-sweep/test-generator.R
library(testthat)
source("dev/simstudy/spatial-design-sweep/generator.R")

test_that("seeds are deterministic integers distinct across labels", {
  expect_identical(sweep_seed("spatial-design-data", 1L), sweep_seed("spatial-design-data", 1L))
  expect_true(sweep_seed("spatial-design-data", 1L) != sweep_seed("spatial-design-sites", 1L))
  expect_true(sweep_seed("spatial-design-data", 2L) == sweep_seed("spatial-design-data", 1L) + 1L)
})

test_that("arrangements have 100 sites each inside the unit square", {
  a <- make_arrangements(sweep_seed("spatial-design-sites", 1L))
  expect_identical(names(a), c("spread", "pairs", "clustered", "grid"))
  for (nm in names(a)) {
    expect_identical(dim(a[[nm]]), c(100L, 2L))
    expect_identical(colnames(a[[nm]]), c("x", "y"))
    expect_true(all(a[[nm]] >= 0 & a[[nm]] <= 1))
    expect_false(anyDuplicated(a[[nm]]) > 0)
    ratio <- sd(a[[nm]][, 1]) / sd(a[[nm]][, 2])
    expect_true(abs(ratio - 1) <= 0.1)
  }
})

test_that("pairs extend the first 80 spread sites with partners at distance 0.01", {
  a <- make_arrangements(sweep_seed("spatial-design-sites", 1L))
  expect_identical(a$pairs[1:80, ], a$spread[1:80, ])
  d <- as.matrix(dist(a$pairs))[81:100, 1:80]
  expect_equal(apply(d, 1, min), rep(0.01, 20), tolerance = 1e-12)
})

test_that("cluster centres are separated and points lie within 0.02 of a centre", {
  a <- make_arrangements(sweep_seed("spatial-design-sites", 1L))
  centres <- attr(a, "cluster_centres")
  expect_identical(dim(centres), c(10L, 2L))
  dc <- as.matrix(dist(centres)); diag(dc) <- Inf
  expect_true(min(dc) >= 0.2)
  for (k in 1:10) {
    pts <- a$clustered[(k - 1) * 10 + 1:10, ]
    d <- sqrt(colSums((t(pts) - centres[k, ])^2))
    expect_true(all(d <= 0.02 + 1e-12))
  }
})

test_that("the grid is the exact 10 by 10 lattice", {
  a <- make_arrangements(1L)
  g <- (1:10 - 0.5) / 10
  expect_equal(unname(a$grid), unname(as.matrix(expand.grid(x = g, y = g))))
})

test_that("the landscape has the right pieces and target prevalence on the lattice", {
  land <- make_landscape(1L)
  expect_identical(dim(land$points), c(1920L, 2L))
  expect_identical(land$index$lattice, 321:1920)
  expect_identical(land$index$pairs, c(1:80, 101:120))
  expect_identical(dim(land$field), c(1920L, 8L))
  expect_equal(colMeans(land$psi[land$index$lattice, ]),
               c(.05, .05, .25, .25, .25, .75, .75, .75), tolerance = 1e-6)
  expect_equal(land$range, 0.03)
  expect_true(all(land$z %in% 0:1))
})

test_that("the generator is reproducible", {
  expect_identical(make_sweep_input(2L), make_sweep_input(2L))
})

test_that("surveys map rows to sites, samples and primers consistently", {
  input <- make_sweep_input(1L)
  for (arr in c("spread", "pairs", "clustered", "grid")) {
    s <- input$surveys[[arr]]
    expect_identical(dim(s$two_stage$OTU), c(2400L, 8L))
    expect_identical(nrow(s$two_stage$info), 2400L)
    expect_identical(unname(s$two_stage$info$Site), rep(1:100, each = 24))
    expect_identical(unname(s$two_stage$info$Sample), rep(1:200, each = 12))
    expect_identical(unname(s$two_stage$info$Primer), rep(rep(1:2, each = 6), times = 200))
    expect_identical(unname(s$binary$OTU), unname(input$landscape$z[input$landscape$index[[arr]], ]))
    expect_identical(colnames(s$binary$OTU), input$landscape$species)
    expect_equal(s$binary$info$longitude, input$landscape$points[input$landscape$index[[arr]], 1])
  }
})

test_that("design statistics are computed for every arrangement", {
  input <- make_sweep_input(1L)
  st <- input$statistics
  expect_identical(st$arrangement, c("spread", "pairs", "clustered", "grid"))
  expect_true(all(c("mean_nearest_neighbour", "fraction_with_half_neighbour",
                    "effective_rank", "standardised_range_x", "standardised_range_y",
                    "axis_sd_ratio", "environment_span", "in_grid") %in% names(st)))
  expect_true(st$fraction_with_half_neighbour[st$arrangement == "clustered"] >
              st$fraction_with_half_neighbour[st$arrangement == "spread"])
  expect_equal(st$mean_nearest_neighbour[st$arrangement == "grid"], 0.1)
})
cat("Generator tests passed.\n")
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep && Rscript dev/simstudy/spatial-design-sweep/test-generator.R`
Expected: FAIL with `cannot open file 'dev/simstudy/spatial-design-sweep/generator.R'`.

- [ ] **Step 4: Write the generator**

Create `dev/simstudy/spatial-design-sweep/generator.R`:

```r
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
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep && Rscript dev/simstudy/spatial-design-sweep/test-generator.R`
Expected: `Generator tests passed.` with no failures. If `make_landscape` is slow, it is the two 1920 by 1920 Cholesky factorisations; a few seconds is normal.

- [ ] **Step 6: Print the design table for the three communities and check it against the spec**

Run:

```sh
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep && Rscript -e 'source("dev/simstudy/spatial-design-sweep/generator.R"); st <- do.call(rbind, lapply(1:3, function(r) make_sweep_input(r)$statistics)); print(st, digits = 3)'
```

Expected: `in_grid` TRUE for all twelve rows; `fraction_with_half_neighbour` near 0.3 for spread, near 0.5 for pairs, 1.0 for clustered, 0 for grid; `axis_sd_ratio` within 0.9 to 1.1. If any `in_grid` is FALSE, stop and report; the spec flags such a design rather than adjusting it.

- [ ] **Step 7: Commit**

```bash
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
git add dev/simstudy/spatial-design-sweep/generator.R dev/simstudy/spatial-design-sweep/test-generator.R dev/simstudy/spatial-design-sweep/PLAN.md
git commit -m "Add the sweep generator with design statistics and the anisotropy amendment

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: Oracle module, prepare and oracle run modes

**Files:**
- Create: `dev/simstudy/spatial-design-sweep/oracle.R`
- Create: `dev/simstudy/spatial-design-sweep/test-oracle.R`
- Create: `dev/simstudy/spatial-design-sweep/run.R` (modes `prepare` and `oracle`; Task 3 adds the rest)

**Interfaces:**
- Consumes: `make_sweep_input()` from Task 1; `ellipse_draws_cpp(lower, y, trials, offset, intercept, initial, nburn, niter)` from the diagnosis; `robust_trace_diagnostics(x)` and `robust_flag_rows(d)` from `spatial-amplitude-prior/robust.R`.
- Produces: `oracle_source_hashes(repo)`; `compile_oracle(repo, cache = NULL)`; `oracle_covariance(xy, range, jitter = 1e-8)`; `run_oracle(input, arrangement, species, seed, nburn, niter, nchain = 4L)` returning `list(draws, mean_proposals, offset, lower)`; `score_oracle(draws, truth_field, psi, offset)` returning `list(metrics, diagnostics, field_median, probability_mean)`; `oracle_lattice(draws, xy_sites, xy_lattice, range, offset_lattice, thin = 4L, seed)` returning `list(probability_mean, field_mean, draws_used)`; `run.R --mode=prepare` writing `STUDY/inputs/repNN.rds`, `STUDY/design-statistics.csv`, `STUDY/input-manifest.csv`; `run.R --mode=oracle` writing `STUDY/oracle/repNN-ARR-speciesSS-PHASE.rds` and `STUDY/oracle/oracle-selected.csv`.

- [ ] **Step 1: Write the failing oracle tests**

The MD5 of the unmodified oracle source on 1 October 2026 is `8154a611d1b7ab926eb32bf441e9af3e`; the test below pins it.

Create `dev/simstudy/spatial-design-sweep/test-oracle.R`:

```r
# Run from the repository root: Rscript dev/simstudy/spatial-design-sweep/test-oracle.R
library(testthat)
source("dev/simstudy/spatial-design-sweep/generator.R")
source("dev/simstudy/spatial-design-sweep/oracle.R")
ELLIPSE_MD5 <- "8154a611d1b7ab926eb32bf441e9af3e"

test_that("the oracle sampler source is the unmodified diagnosis file", {
  h <- oracle_source_hashes(".")
  expect_identical(unname(h[basename(names(h)) == "ellipse.cpp"]), ELLIPSE_MD5)
})

test_that("the oracle covariance is the true isotropic kernel plus jitter", {
  xy <- make_lattice(5L)
  Q <- oracle_covariance(xy, 0.3)
  expect_equal(Q, sq_exp(xy, xy, 0.3) + diag(1e-8, 25), tolerance = 1e-12)
  expect_true(min(eigen(Q, symmetric = TRUE, only.values = TRUE)$values) > 0)
})

test_that("the sampler reproduces a one-site logistic-normal posterior mean", {
  compile_oracle(".")
  y <- 1; offset <- -0.5; sd <- 1.2
  lower <- matrix(sd, 1, 1)
  set.seed(1)
  initial <- cbind(0, lower %*% matrix(rnorm(3), 1, 3))
  fit <- ellipse_draws_cpp(lower, y, 1L, offset, FALSE, initial, 2000L, 20000L)
  draws <- as.vector(fit$draws)
  numerator <- integrate(function(f) f * plogis(offset + f) * dnorm(f, 0, sd), -Inf, Inf)$value
  denominator <- integrate(function(f) plogis(offset + f) * dnorm(f, 0, sd), -Inf, Inf)$value
  expect_equal(mean(draws), numerator / denominator, tolerance = 0.02)
})

test_that("run_oracle and score_oracle return the agreed shapes", {
  compile_oracle(".")
  input <- make_sweep_input(1L)
  res <- run_oracle(input, "grid", 6L, seed = 5L, nburn = 50L, niter = 100L)
  expect_identical(dim(res$draws), c(100L, 100L, 4L))
  tr <- input$surveys$grid$truth
  sc <- score_oracle(res$draws, tr$field[, 6], tr$psi[, 6], res$offset)
  expect_true(all(c("centred_rmse", "centred_correlation", "centred_slope", "zero_field_rmse",
                    "occupancy_bias", "occupancy_mae", "flag_count") %in% names(sc$metrics)))
  expect_length(sc$field_median, 100L)
  lat <- oracle_lattice(res$draws, tr$xy, input$landscape$points[input$landscape$index$lattice, ],
                        input$landscape$range,
                        input$landscape$B0[6] + input$landscape$B[6] * input$landscape$environment[input$landscape$index$lattice],
                        thin = 10L, seed = 7L)
  expect_length(lat$probability_mean, 1600L)
  expect_true(all(lat$probability_mean > 0 & lat$probability_mean < 1))
})
cat("Oracle tests passed.\n")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep && Rscript dev/simstudy/spatial-design-sweep/test-oracle.R`
Expected: FAIL with `cannot open file 'dev/simstudy/spatial-design-sweep/oracle.R'`.

- [ ] **Step 3: Write the oracle module**

Create `dev/simstudy/spatial-design-sweep/oracle.R`:

```r
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
```

- [ ] **Step 4: Write run.R with the prepare and oracle modes**

Create `dev/simstudy/spatial-design-sweep/run.R`:

```r
#!/usr/bin/env Rscript
# Lesson 2 sweep driver. Modes: prepare, oracle, freeze, pilot, initial, long.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", name)
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }
  substring(hit, nchar(name) + 4L)
}
known <- c("repo", "study", "mode", "workers", "keys")
stopifnot(all(sub("^--([^=]+)=.*$", "\\1", args) %in% known))
repo <- normalizePath(option("repo", ".")); study <- option("study")
mode <- option("mode"); workers <- as.integer(option("workers", "2"))
stopifnot(mode %in% c("prepare", "oracle", "freeze", "pilot", "initial", "long"), workers %in% 1:8)
dir.create(study, showWarnings = FALSE, recursive = TRUE); study <- normalizePath(study)
scripts <- file.path(repo, "dev/simstudy/spatial-design-sweep")
source(file.path(scripts, "generator.R"))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
atomic <- function(object, path) {
  tmp <- paste0(path, ".tmp"); saveRDS(object, tmp); stopifnot(file.rename(tmp, path))
}
log <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n"); flush.console() }
arrangements <- c("spread", "pairs", "clustered", "grid")
generator_hash <- tools::md5sum(file.path(scripts, "generator.R"))

# ---- prepare --------------------------------------------------------------
dir.create(file.path(study, "inputs"), showWarnings = FALSE)
inputs <- character()
for (r in 1:3) {
  community <- sprintf("rep%02d", r)
  path <- file.path(study, "inputs", paste0(community, ".rds"))
  input <- make_sweep_input(r)
  if (file.exists(path)) stopifnot(identical(readRDS(path), input)) else atomic(input, path)
  inputs[community] <- path
}
statistics <- do.call(rbind, lapply(inputs, function(p) readRDS(p)$statistics))
rownames(statistics) <- NULL
stat_path <- file.path(study, "design-statistics.csv")
if (!file.exists(stat_path)) write.csv(statistics, stat_path, row.names = FALSE)
if (any(!statistics$in_grid)) log("WARNING: a design's standardised range lies outside the fitter's grid; see design-statistics.csv")
manifest <- expand.grid(replicate = 1:3, arrangement = arrangements, arm = c("binary", "two_stage"),
                        stringsAsFactors = FALSE)
manifest$community <- sprintf("rep%02d", manifest$replicate)
manifest$key <- sprintf("%s-%s-%s", manifest$community, manifest$arrangement, manifest$arm)
manifest$input_md5 <- unname(tools::md5sum(inputs[manifest$community]))
manifest$fit_seed <- mapply(function(a, r, arr) sweep_seed(paste0("spatial-design-fit-", a, "-", arr), r),
                            manifest$arm, manifest$replicate, manifest$arrangement)
manifest <- manifest[order(manifest$replicate, match(manifest$arrangement, arrangements), manifest$arm), ]
rownames(manifest) <- NULL
manifest_path <- file.path(study, "input-manifest.rds")
if (file.exists(manifest_path)) stopifnot(identical(readRDS(manifest_path), manifest)) else {
  atomic(manifest, manifest_path); write.csv(manifest, file.path(study, "input-manifest.csv"), row.names = FALSE)
}
if (mode == "prepare") { log("Prepared 3 communities, 4 arrangements, 24 fit jobs; generator md5", generator_hash); quit(status = 0) }

# ---- oracle ---------------------------------------------------------------
if (mode == "oracle") {
  source(file.path(scripts, "oracle.R"))
  out <- file.path(study, "oracle"); dir.create(out, showWarnings = FALSE)
  hashes <- compile_oracle(repo, cache = file.path(out, "cpp-cache"))
  jobs <- expand.grid(community = names(inputs), arrangement = arrangements, species = 1:8,
                      stringsAsFactors = FALSE)
  work <- function(i) {
    job <- jobs[i, ]; input <- readRDS(inputs[job$community]); tr <- input$surveys[[job$arrangement]]$truth
    key <- sprintf("%s-%s-species%02d", job$community, job$arrangement, job$species)
    seed <- sweep_seed(paste0("spatial-design-oracle-", job$arrangement, "-", job$species), input$replicate)
    perform <- function(phase) {
      dest <- file.path(out, paste0(key, "-", phase, ".rds"))
      nburn <- if (phase == "initial") 1000L else 2000L
      niter <- if (phase == "initial") 2000L else 4000L
      if (file.exists(dest)) { res <- readRDS(dest); stopifnot(identical(res$hashes, hashes), res$seed == seed); return(res) }
      started <- Sys.time()
      fit <- run_oracle(input, job$arrangement, job$species, seed, nburn, niter)
      scored <- score_oracle(fit$draws, tr$field[, job$species], tr$psi[, job$species], fit$offset)
      lattice_index <- input$landscape$index$lattice
      offset_lattice <- input$landscape$B0[job$species] + input$landscape$B[job$species] * input$landscape$environment[lattice_index]
      lat <- oracle_lattice(fit$draws, tr$xy, input$landscape$points[lattice_index, ], input$landscape$range,
                            offset_lattice, thin = 4L, seed = seed + 1L)
      res <- c(list(job = job, key = key, phase = phase, seed = seed, nburn = nburn, niter = niter, nchain = 4L,
                    hashes = hashes, generator_hash = generator_hash, input_md5 = unname(tools::md5sum(inputs[job$community])),
                    elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")),
                    draws = fit$draws, mean_proposals = fit$mean_proposals, lattice = lat), scored)
      atomic(res, dest); log(key, phase, "seconds", round(res$elapsed), "flags", res$metrics$flag_count); res
    }
    a <- perform("initial"); initial_flags <- a$metrics$flag_count
    if (initial_flags > 0L) a <- perform("long")
    data.frame(job, key = key, phase = a$phase, initial_flag_count = initial_flags, a$metrics,
               file = file.path(out, paste0(key, "-", a$phase, ".rds")))
  }
  results <- if (workers == 1L) lapply(seq_len(nrow(jobs)), work) else
    parallel::mclapply(seq_len(nrow(jobs)), work, mc.cores = min(workers, 2L), mc.preschedule = FALSE, mc.set.seed = FALSE)
  stopifnot(!any(vapply(results, inherits, logical(1), "try-error")))
  table <- do.call(rbind, results); table$result_md5 <- unname(tools::md5sum(table$file))
  write.csv(table, file.path(out, "oracle-selected.csv"), row.names = FALSE)
  stopifnot(nrow(table) == 96L)
  log("oracle finished; 96 posteriors;", sum(table$flag_count > 0), "still flagged after the doubled rerun")
  quit(status = 0)
}
stop("Mode ", mode, " is added in Task 3")
```

- [ ] **Step 5: Run the oracle tests to verify they pass**

Run: `cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep && Rscript dev/simstudy/spatial-design-sweep/test-oracle.R`
Expected: `Oracle tests passed.` The first run compiles the sampler; allow a minute.

- [ ] **Step 6: Prepare the study directory and run the oracle on all cells**

Run:

```sh
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
STUDY=$PWD/dev/simstudy/results/spatial-design-20261001
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=prepare
nohup Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=oracle --workers=2 > $STUDY/oracle.log 2>&1 &
```

Then wait with `until grep -q "oracle finished" $STUDY/oracle.log; do sleep 60; done` in a background Bash call. Expected: `oracle finished; 96 posteriors; N still flagged`. Inspect `$STUDY/oracle/oracle-selected.csv`: `centred_rmse / zero_field_rmse` should be near 1 for spread and grid and clearly below 1 for clustered and pairs, lower for the 25% and 75% species than the 5% ones. Report the per-arrangement means to Doug before continuing; this is the first checkpoint the spec promised.

- [ ] **Step 7: Commit**

```bash
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
git add dev/simstudy/spatial-design-sweep/oracle.R dev/simstudy/spatial-design-sweep/test-oracle.R dev/simstudy/spatial-design-sweep/run.R
git commit -m "Add the oracle module and the prepare and oracle run modes

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: Scoring module, freeze and the full fits

**Files:**
- Create: `dev/simstudy/spatial-design-sweep/score.R`
- Create: `dev/simstudy/spatial-design-sweep/test-score.R`
- Modify: `dev/simstudy/spatial-design-sweep/run.R` (replace the final `stop()` with the `freeze`, `pilot`, `initial`, `long` modes)

**Interfaces:**
- Consumes: from `dev/simstudy/spatial-targeted-recheck/score.R`: `trace_diagnostics(x)`, `independent_bases(fit)`, `reconstruct_spatial_draws(fit)`, `score_draw_block(draws, truth, metric, group)`, `verify_saved_probability(estimate, stored)`. From the frozen package: `occJSDM:::transform_new_covariates(df, list_matrix, remove_intercept = TRUE)` (defined in `R/jsdmfun.R`), `occJSDM:::createSpatialPredMatrix(Xs, l_s_grid, X_tilde, list_Xs_mat)`, `occJSDM:::precomputeSORmatrices`, `occJSDM:::KsBproduct`, `occJSDM::predictNewSites`.
- Produces: `score_source_hashes(repo)`; `reconstruct_field_draws(fit, thin = 4L)` returning array `n x S x draws`; `reconstruct_lattice_draws(fit, lattice_environment, lattice_xy, thin = 4L)` returning `list(with, without, draws_used, Ks, X_new)`; `sweep_groups(psi, prevalence)`; `distance_bins(lattice_xy, site_xy)`; `score_sweep_fit(fit, input, arrangement, arm, lattice_subset = 50L)` returning `list(groups, elements, species, field, range, range_summary, amplitude, lattice, probability, field_median, lattice_mean, lattice_mean_nospatial, lattice_distance, verification)`; `run.R` modes `freeze`, `pilot`, `initial`, `long` writing `STUDY/<mode>/KEY-fit.rds` and `KEY-result.rds`.

- [ ] **Step 1: Write the failing score tests**

Create `dev/simstudy/spatial-design-sweep/test-score.R`:

```r
# Run from the repository root with the installed occJSDM available.
library(testthat)
suppressPackageStartupMessages(library(occJSDM))
source("dev/simstudy/spatial-design-sweep/generator.R")
source("dev/simstudy/spatial-design-sweep/score.R")

test_that("sweep groups partition sites by band and species by prevalence", {
  psi <- matrix(c(.1, .5, .9, .3), 2, 2); prevalence <- c(.05, .75)
  g <- sweep_groups(psi, prevalence)
  expect_true(all(g$low + g$medium + g$high == 1))
  expect_identical(unname(g$prevalence_5pct[, 1]), c(TRUE, TRUE))
  expect_identical(unname(g$prevalence_5pct[, 2]), c(FALSE, FALSE))
})

test_that("a tiny binary fit scores and its lattice prediction matches the native median", {
  input <- make_sweep_input(1L)
  data <- input$surveys$grid$binary
  set.seed(3)
  fit <- suppressMessages(suppressWarnings(runOccJSDM(data,
    listParams = list(n_factors = 0L, n_lattrait = 0L, n_supportpoints = 100L),
    occCovariates = "environment", spatCovariates = c("longitude", "latitude"),
    MCMCparams = list(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L))))
  expect_identical(fit$infos$ps, 100L)
  sc <- score_sweep_fit(fit, input, "grid", "binary")
  expect_true(all(c("occupancy", "intercept", "environment_slope", "range", "spatial_sd") %in% sc$groups$metric))
  expect_identical(dim(sc$field_median), c(100L, 8L))
  expect_identical(dim(sc$lattice_mean), c(1600L, 8L))
  expect_true(sc$verification$lattice < 1e-8)
  expect_true(sc$verification$basis < 1e-10)
  expect_true(all(c("bin", "mae", "bias") %in% names(sc$lattice)))
  expect_identical(sort(unique(sc$lattice$bin)), c("0.02 to 0.05", "0.05 to 0.1", "above 0.1", "all", "up to 0.02"))
})
cat("Score tests passed.\n")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep && Rscript dev/simstudy/spatial-design-sweep/test-score.R`
Expected: FAIL with `cannot open file 'dev/simstudy/spatial-design-sweep/score.R'`.

- [ ] **Step 3: Write the scoring module**

Create `dev/simstudy/spatial-design-sweep/score.R`:

```r
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
```

- [ ] **Step 4: Run the score tests to verify they pass**

Run: `cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep && Rscript dev/simstudy/spatial-design-sweep/test-score.R`
Expected: `Score tests passed.` The tiny fit takes under a minute. `transform_new_covariates()` is defined in `R/jsdmfun.R` on the frozen main.

- [ ] **Step 5: Add the freeze, pilot, initial and long modes to run.R**

Replace the final line `stop("Mode ", mode, " is added in Task 3")` in `run.R` with:

```r
# ---- freeze ---------------------------------------------------------------
if (mode == "freeze") {
  stopifnot(!dir.exists(file.path(study, "library")))
  dirty <- system2("git", c("-C", shQuote(repo), "status", "--porcelain", "--", "R", "src", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)
  stopifnot(length(dirty) == 0L)
  rev <- system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)
  src <- file.path(study, "source-main"); dir.create(src)
  for (d in c("R", "src", "man", "data", "inst")) if (dir.exists(file.path(repo, d))) file.copy(file.path(repo, d), src, recursive = TRUE)
  file.copy(file.path(repo, c("DESCRIPTION", "NAMESPACE")), src)
  unlink(list.files(file.path(src, "src"), pattern = "\\.(o|so|dll)$", full.names = TRUE))
  dir.create(file.path(study, "library"))
  status <- system2("R", c("CMD", "INSTALL", "--preclean", paste0("--library=", shQuote(file.path(study, "library"))), shQuote(src)),
                    stdout = file.path(study, "install.log"), stderr = file.path(study, "install.log"))
  stopifnot(status == 0L)
  writeLines(rev, file.path(study, "source-revision.txt"))
  cat(sprintf("- **%s, amendment 2.** Production code frozen at main revision `%s` and installed into the study library before the first full fit.\n",
              format(Sys.Date(), "%d %B %Y"), rev), file = file.path(scripts, "PLAN.md"), append = TRUE)
  log("Frozen revision", rev, "installed into", file.path(study, "library"))
  quit(status = 0)
}

# ---- fits -----------------------------------------------------------------
.libPaths(c(file.path(study, "library"), .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == normalizePath(file.path(study, "library/occJSDM")))
RcppParallel::setThreadOptions(numThreads = 1)
source(file.path(scripts, "score.R"))
production <- c(file.path(study, "source-main", c("DESCRIPTION", "NAMESPACE")),
                list.files(file.path(study, "source-main/R"), pattern = "\\.R$", full.names = TRUE),
                list.files(file.path(study, "source-main/src"), pattern = "\\.(cpp|h)$|^Makevars", full.names = TRUE),
                file.path(find.package("occJSDM"), "libs/occJSDM.so"))
fit_hashes <- tools::md5sum(c(production, file.path(scripts, c("generator.R", "run.R"))))
score_hash <- score_source_hashes(repo)
jobs <- split(manifest, manifest$key)[manifest$key]
keys <- option("keys", "")
if (nzchar(keys)) { keys <- strsplit(keys, ",", fixed = TRUE)[[1]]; stopifnot(all(keys %in% names(jobs))); jobs <- jobs[keys] } else
if (mode == "pilot") jobs <- jobs["rep01-clustered-two_stage"] else
if (mode == "long") stop("Long runs require explicit --keys from summarise.R --mode=select")
mcmc <- switch(mode, pilot = list(nchain = 2L, nburn = 40L, niter = 60L, nthin = 1L),
               initial = list(nchain = 2L, nburn = 3000L, niter = 5000L, nthin = 1L),
               long = list(nchain = 4L, nburn = 6000L, niter = 12000L, nthin = 1L))
out <- file.path(study, mode); dir.create(out, showWarnings = FALSE)
settings <- list(revision = readLines(file.path(study, "source-revision.txt")), fit_hashes = fit_hashes,
                 score_hash = score_hash, mcmc = mcmc, mode = mode, workers = workers, threads_per_fit = 1L,
                 priors = "Unchanged defaults", session = sessionInfo())
settings_file <- file.path(out, "settings.rds")
if (file.exists(settings_file)) { old <- readRDS(settings_file)
  stopifnot(identical(old$fit_hashes, fit_hashes), identical(old$mcmc, mcmc), identical(old$score_hash, score_hash)) } else atomic(settings, settings_file)
run_job <- function(job) {
  input_file <- inputs[job$community]
  stopifnot(identical(unname(tools::md5sum(input_file)), job$input_md5))
  dest <- file.path(out, paste0(job$key, "-result.rds")); fitfile <- file.path(out, paste0(job$key, "-fit.rds"))
  if (file.exists(dest)) { old <- readRDS(dest); stopifnot(identical(old$job, job), identical(old$mcmc, mcmc)); return(job$key) }
  input <- readRDS(input_file); data <- input$surveys[[job$arrangement]][[job$arm]]
  log(job$key, "started")
  if (file.exists(fitfile)) { saved <- readRDS(fitfile); stopifnot(identical(saved$job, job), identical(saved$mcmc, mcmc))
    fit <- saved$fit; warnings <- saved$warnings; started <- saved$started; finished <- saved$finished } else {
    set.seed(job$fit_seed); warnings <- character(); started <- Sys.time()
    fit <- withCallingHandlers(suppressMessages(occJSDM::runOccJSDM(data,
      listParams = list(n_factors = 0L, n_lattrait = 0L, n_supportpoints = 100L), listPriors = list(), threshold = 1,
      occCovariates = "environment", collCovariates = if (job$arm == "binary") NULL else "collection",
      spatCovariates = c("longitude", "latitude"), MCMCparams = mcmc)),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    finished <- Sys.time()
    atomic(list(fit = fit, job = job, mcmc = mcmc, fit_hashes = fit_hashes, warnings = warnings, started = started, finished = finished), fitfile)
  }
  scores <- score_sweep_fit(fit, input, job$arrangement, job$arm)
  atomic(c(list(job = job, mcmc = mcmc, fit_hashes = fit_hashes, score_hash = score_hash, warnings = warnings,
                started = started, finished = finished), scores), dest)
  log(job$key, "complete; fit seconds", round(as.numeric(difftime(finished, started, units = "secs"))),
      "; warnings", length(warnings), "; max group Rhat", round(max(scores$groups$rhat, na.rm = TRUE), 3))
  job$key
}
work <- function(job) tryCatch(run_job(job), error = function(e) { log(job$key, "ERROR:", conditionMessage(e)); list(key = job$key, error = conditionMessage(e)) })
status <- if (workers == 1L) lapply(jobs, work) else parallel::mclapply(jobs, work, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
atomic(status, file.path(out, paste0("status-", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")))
stopifnot(all(vapply(status, is.character, logical(1))))
log("Completed", length(status), mode, "fits.")
```

- [ ] **Step 6: Freeze the production code and run the pilot**

Run:

```sh
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
STUDY=$PWD/dev/simstudy/results/spatial-design-20261001
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=freeze
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=pilot --workers=1
```

Expected: `Frozen revision <hash> installed`, then `rep01-clustered-two_stage complete` within a few minutes and `Completed 1 pilot fits.` Confirm `git diff dev/simstudy/spatial-design-sweep/PLAN.md` shows only the appended amendment 2 line. The pilot's scores are not evidence; they check dimensions, truth alignment and runtime.

- [ ] **Step 7: Commit the scoring module, the run modes and the amendment**

```bash
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
git add dev/simstudy/spatial-design-sweep/score.R dev/simstudy/spatial-design-sweep/test-score.R dev/simstudy/spatial-design-sweep/run.R dev/simstudy/spatial-design-sweep/PLAN.md
git commit -m "Add fit scoring with lattice prediction and the freeze, pilot, initial and long modes

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

- [ ] **Step 8: Launch the 24 initial fits detached**

Run:

```sh
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
STUDY=$PWD/dev/simstudy/results/spatial-design-20261001
nohup Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=initial --workers=8 > $STUDY/initial.log 2>&1 &
```

Wait with a Monitor on `$STUDY/initial.log` for lines matching `complete;|ERROR|Completed 24 initial fits` with a 30-minute timeout, re-armed as needed. Expected: 24 `complete` lines and `Completed 24 initial fits.` Report the per-arrangement binary-control field correlations from the result files to Doug; this is the second checkpoint. If a fit errors, read `$STUDY/initial.log`, fix the cause in the script (not the data), delete only that job's partial files and rerun with `--keys=KEY`.

---

### Task 4: Selection, longer runs, audit, figures and study README

**Files:**
- Create: `dev/simstudy/spatial-design-sweep/summarise.R`
- Create: `dev/simstudy/spatial-design-sweep/verify.R`
- Create: `dev/simstudy/spatial-design-sweep/plot.R`
- Create: `dev/simstudy/spatial-design-sweep/README.md`
- Create: `dev/simstudy/spatial-design-sweep/results/` (compact CSVs and PNGs)

**Interfaces:**
- Consumes: `STUDY/initial/*-result.rds`, `STUDY/long/*-result.rds`, `STUDY/oracle/*.rds`, `STUDY/input-manifest.rds`, `STUDY/design-statistics.csv`.
- Produces: `summarise.R --mode=select` writing `STUDY/long-keys.txt` and `STUDY/long-selection.csv`; `summarise.R --mode=final` writing `STUDY/summary-final/{groups,species,field,range,amplitude,lattice,oracle,aggregate,paired,reading,range-reading,selected-fits,long-run-sensitivity}.csv` and `selection.rds`; `verify.R` writing `STUDY/audit/audit.csv` with every difference and a `passed` column; `plot.R` writing PNGs and copying CSVs into `dev/simstudy/spatial-design-sweep/results/`; `diagnostic_reasons(result)`; `read_sweep_selection(study, final)`; `reading_labels(field_rows)`.

- [ ] **Step 1: Write summarise.R**

Create `dev/simstudy/spatial-design-sweep/summarise.R`:

```r
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", name)
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }
  substring(hit, nchar(name) + 4L)
}
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
mode <- option("mode", "select"); stopifnot(mode %in% c("select", "final"))

diagnostic_reasons <- function(a) {
  reasons <- character()
  if (any(a$groups$rhat > 1.05, na.rm = TRUE)) reasons <- c(reasons, "group Rhat > 1.05")
  if (any(a$species$rhat > 1.05, na.rm = TRUE)) reasons <- c(reasons, "species Rhat > 1.05")
  if (any(a$elements$rhat > 1.05, na.rm = TRUE)) reasons <- c(reasons, "element Rhat > 1.05")
  primary <- a$groups$metric == "occupancy"
  if (any(a$groups$ess_mean[primary] < 100, na.rm = TRUE)) reasons <- c(reasons, "occupancy group ESS < 100")
  if (any(a$species$ess_mean < 100, na.rm = TRUE)) reasons <- c(reasons, "occupancy species ESS < 100")
  if (any(!is.finite(a$elements$rhat[a$elements$metric != "range"]))) reasons <- c(reasons, "non-range element Rhat unavailable")
  if (any(grepl("Convergence|Low ESS", a$warnings))) reasons <- c(reasons, "native convergence warning")
  reasons
}

read_sweep_selection <- function(study, final = FALSE) {
  manifest <- readRDS(file.path(study, "input-manifest.rds")); stopifnot(nrow(manifest) == 24L)
  selected <- list(); initial <- list(); rows <- list()
  for (key in manifest$key) {
    f <- file.path(study, "initial", paste0(key, "-result.rds"))
    if (!file.exists(f)) stop("Initial batch incomplete: ", key)
    a <- readRDS(f); stopifnot(identical(a$mcmc, list(nchain = 2L, nburn = 3000L, niter = 5000L, nthin = 1L)))
    initial[[key]] <- a; reasons <- diagnostic_reasons(a); needs_long <- length(reasons) > 0L
    b <- a; phase <- "initial"; selected_file <- f
    if (final && needs_long) {
      lf <- file.path(study, "long", paste0(key, "-result.rds"))
      if (!file.exists(lf)) stop("Prespecified longer check missing: ", key)
      b <- readRDS(lf); phase <- "long"; selected_file <- lf
      stopifnot(identical(b$mcmc, list(nchain = 4L, nburn = 6000L, niter = 12000L, nthin = 1L)), identical(b$job, a$job))
    }
    selected[[key]] <- b
    rows[[key]] <- data.frame(key = key, community = a$job$community, arrangement = a$job$arrangement, arm = a$job$arm,
      replicate = a$job$replicate, needs_long = needs_long, initial_reasons = paste(reasons, collapse = "; "),
      phase = phase, selected_file = selected_file, selected_md5 = unname(tools::md5sum(selected_file)),
      selected_reasons = paste(diagnostic_reasons(b), collapse = "; "),
      max_group_rhat = max(b$groups$rhat, na.rm = TRUE), warnings = paste(b$warnings, collapse = " | "))
  }
  list(manifest = do.call(rbind, rows), selected = selected, initial = initial)
}

collect <- function(results, table) do.call(rbind, lapply(names(results), function(key) {
  r <- results[[key]]; d <- r[[table]]; if (is.null(d)) return(NULL)
  data.frame(key = key, community = r$job$community, arrangement = r$job$arrangement, arm = r$job$arm,
             replicate = r$job$replicate, phase = if (identical(r$mcmc$nchain, 4L)) "long" else "initial", d) }))

reading_labels <- function(field_rows) {
  field_rows$reduction <- 1 - field_rows$centred_rmse / field_rows$zero_field_rmse
  cells <- split(field_rows, interaction(field_rows$arrangement, field_rows$group, drop = TRUE))
  do.call(rbind, lapply(cells, function(d) {
    stopifnot(nrow(d) == 3L)
    label <- if (all(d$reduction >= .2 & d$centred_correlation >= .5)) "informative" else
             if (all(d$reduction < .1 | d$centred_correlation < .3)) "uninformative" else "intermediate"
    data.frame(arrangement = d$arrangement[1], group = d$group[1], label = label,
               min_reduction = min(d$reduction), min_correlation = min(d$centred_correlation)) }))
}

selection <- read_sweep_selection(study, final = mode == "final")
if (mode == "select") {
  write.csv(selection$manifest, file.path(study, "long-selection.csv"), row.names = FALSE)
  keys <- selection$manifest$key[selection$manifest$needs_long]
  writeLines(keys, file.path(study, "long-keys.txt"))
  cat(length(keys), "of 24 fits require the prespecified longer checks.\n"); quit(status = 0)
}
out <- file.path(study, "summary-final")
if (dir.exists(out) && length(list.files(out))) stop("Use an empty output directory: ", out)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tables <- list(groups = collect(selection$selected, "groups"), species = collect(selection$selected, "species"),
               field = collect(selection$selected, "field"), range = collect(selection$selected, "range_summary"),
               amplitude = collect(selection$selected, "amplitude"), lattice = collect(selection$selected, "lattice"))
for (t in names(tables)) write.csv(tables[[t]], file.path(out, paste0(t, ".csv")), row.names = FALSE)
oracle <- read.csv(file.path(study, "oracle/oracle-selected.csv"))
oracle$group <- paste0("prevalence_", c(5, 5, 25, 25, 25, 75, 75, 75)[oracle$species], "pct")
write.csv(oracle, file.path(out, "oracle.csv"), row.names = FALSE)
field <- tables$field; field$group <- paste0("prevalence_", field$target * 100, "pct")
field_cells <- aggregate(cbind(centred_rmse, centred_correlation, zero_field_rmse) ~ community + arrangement + arm + group, field, mean)
oracle_cells <- aggregate(cbind(centred_rmse, centred_correlation, zero_field_rmse) ~ community + arrangement + group, oracle, mean)
oracle_cells$arm <- "oracle"
cells <- rbind(field_cells, oracle_cells[names(field_cells)])
aggregate_cells <- do.call(rbind, lapply(split(cells, interaction(cells$arrangement, cells$arm, cells$group, drop = TRUE)), function(d)
  data.frame(arrangement = d$arrangement[1], arm = d$arm[1], group = d$group[1], n = nrow(d),
             centred_rmse = mean(d$centred_rmse), centred_rmse_min = min(d$centred_rmse), centred_rmse_max = max(d$centred_rmse),
             correlation = mean(d$centred_correlation), correlation_min = min(d$centred_correlation), correlation_max = max(d$centred_correlation),
             zero_field_rmse = mean(d$zero_field_rmse))))
write.csv(aggregate_cells, file.path(out, "aggregate.csv"), row.names = FALSE)
reading <- reading_labels(field_cells[field_cells$arm == "binary", ])
range_reading <- aggregate(mass_within_one_step ~ arrangement + arm, tables$range, function(v) all(v >= .5))
names(range_reading)[3] <- "range_recovered"
write.csv(reading, file.path(out, "reading.csv"), row.names = FALSE)
write.csv(range_reading, file.path(out, "range-reading.csv"), row.names = FALSE)
pairs <- merge(cells[cells$arm == "oracle", ], cells[cells$arm == "binary", ], by = c("community", "arrangement", "group"), suffixes = c("_oracle", "_binary"))
pairs <- merge(pairs, cells[cells$arm == "two_stage", ], by = c("community", "arrangement", "group"))
paired <- data.frame(pairs[c("community", "arrangement", "group")],
                     estimation_cost = pairs$centred_rmse_binary - pairs$centred_rmse_oracle,
                     detection_cost = pairs$centred_rmse - pairs$centred_rmse_binary)
write.csv(paired, file.path(out, "paired.csv"), row.names = FALSE)
write.csv(selection$manifest, file.path(out, "selected-fits.csv"), row.names = FALSE)
initial_groups <- collect(selection$initial, "groups")
sens <- merge(initial_groups, tables$groups, by = c("key", "metric", "group"), suffixes = c("_initial", "_selected"))
sens$bias_change <- sens$bias_selected - sens$bias_initial; sens$mae_change <- sens$mae_selected - sens$mae_initial
write.csv(sens[c("key", "metric", "group", "bias_change", "mae_change")], file.path(out, "long-run-sensitivity.csv"), row.names = FALSE)
saveRDS(selection$manifest, file.path(out, "selection.rds"))
saveRDS(list(scripts = tools::md5sum(file.path(repo, "dev/simstudy/spatial-design-sweep", c("summarise.R", "score.R", "generator.R", "oracle.R"))),
             session = sessionInfo(), time = Sys.time()), file.path(out, "provenance.rds"))
cat("Wrote final summary to", out, "\n")
```

- [ ] **Step 2: Run selection, launch the longer runs, then the final summary**

Run:

```sh
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
STUDY=$PWD/dev/simstudy/results/spatial-design-20261001
Rscript dev/simstudy/spatial-design-sweep/summarise.R --repo=. --study=$STUDY --mode=select
KEYS=$(paste -sd, $STUDY/long-keys.txt)
[ -n "$KEYS" ] && nohup Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=long --workers=8 --keys=$KEYS > $STUDY/long.log 2>&1 &
```

Wait for `Completed N long fits.` in `$STUDY/long.log` (Monitor, re-armed), then:

```sh
Rscript dev/simstudy/spatial-design-sweep/summarise.R --repo=. --study=$STUDY --mode=final
```

Expected: `Wrote final summary to .../summary-final`. Read `reading.csv`, `range-reading.csv`, `aggregate.csv` and `paired.csv` and report them to Doug before writing any lesson text.

- [ ] **Step 3: Write verify.R and run the audit**

Create `dev/simstudy/spatial-design-sweep/verify.R`:

```r
#!/usr/bin/env Rscript
# Independent audit: rebuild probabilities, field medians, range table and lattice
# predictions from saved draws and compare with the scored values.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) { hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }; substring(hit, nchar(name) + 4L) }
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
.libPaths(c(file.path(study, "library"), .libPaths())); suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == normalizePath(file.path(study, "library/occJSDM")))
source(file.path(repo, "dev/simstudy/spatial-design-sweep/generator.R"))
source(file.path(repo, "dev/simstudy/spatial-design-sweep/score.R"))
sel <- read.csv(file.path(study, "summary-final/selected-fits.csv"))
rows <- list()
for (i in seq_len(nrow(sel))) {
  r <- readRDS(sel$selected_file[i]); stopifnot(unname(tools::md5sum(sel$selected_file[i])) == sel$selected_md5[i])
  fitfile <- sub("-result.rds$", "-fit.rds", sel$selected_file[i]); fit <- readRDS(fitfile)$fit
  input <- readRDS(file.path(study, "inputs", paste0(r$job$community, ".rds")))
  stopifnot(unname(tools::md5sum(file.path(study, "inputs", paste0(r$job$community, ".rds")))) == r$job$input_md5)
  js <- fit$results_output$jsdm_output; H <- independent_bases(fit); n <- nrow(fit$Xs); S <- 8L
  ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]; acc <- matrix(0, n, S)
  for (ch in seq_len(nc)) for (it in seq_len(ni)) acc <- acc + plogis(sweep(fit$X_psi %*% matrix(js$B_output[, , it, ch], 1, S) +
    H[[js$idx_ls_output[it, ch]]] %*% matrix(js$Bs_output[, , it, ch], 100L, S), 2, js$B0_output[, it, ch], "+"))
  prob_diff <- max(abs(acc / (ni * nc) - r$probability))
  fm <- apply(reconstruct_field_draws(fit, thin = 4L), c(1, 2), median); field_diff <- max(abs(fm - r$field_median))
  grid <- fit$infos$l_s_grid
  range_diff <- max(abs(vapply(seq_len(nc), function(ch) tabulate(js$idx_ls_output[, ch], length(grid)) / ni, numeric(length(grid))) -
                      matrix(r$range$frequency, length(grid), nc)))
  land <- input$landscape; li <- land$index$lattice
  lat <- reconstruct_lattice_draws(fit, land$environment[li], land$points[li, ], thin = 4L)
  lattice_diff <- max(abs(lat$with - r$lattice_mean), abs(lat$without - r$lattice_mean_nospatial))
  rows[[i]] <- data.frame(key = sel$key[i], phase = sel$phase[i], probability = prob_diff, field_median = field_diff,
    range = range_diff, lattice = lattice_diff, native_lattice = r$verification$lattice, basis = r$verification$basis,
    passed = prob_diff < 1e-10 & field_diff < 1e-10 & range_diff < 1e-12 & lattice_diff < 1e-10 & r$verification$lattice < 1e-8)
  cat(sel$key[i], if (rows[[i]]$passed) "ok" else "MISMATCH", "\n")
}
oracle <- read.csv(file.path(study, "oracle/oracle-selected.csv"))
oracle_ok <- all(unname(tools::md5sum(oracle$file)) == oracle$result_md5)
audit <- do.call(rbind, rows); dir.create(file.path(study, "audit"), showWarnings = FALSE)
write.csv(audit, file.path(study, "audit/audit.csv"), row.names = FALSE)
writeLines(c(paste("oracle result hashes match:", oracle_ok), paste("fits passed:", sum(audit$passed), "of", nrow(audit))),
           file.path(study, "audit/summary.txt"))
stopifnot(oracle_ok, all(audit$passed)); cat("Audit passed for", nrow(audit), "fits and 96 oracle results.\n")
```

Run: `Rscript dev/simstudy/spatial-design-sweep/verify.R --repo=. --study=$STUDY`. Expected: `Audit passed for 24 fits and 96 oracle results.`

- [ ] **Step 4: Write plot.R and generate the figures**

Create `dev/simstudy/spatial-design-sweep/plot.R`:

```r
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) { hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }; substring(hit, nchar(name) + 4L) }
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
out <- file.path(repo, "dev/simstudy/spatial-design-sweep/results"); dir.create(out, showWarnings = FALSE)
summary <- file.path(study, "summary-final")
for (f in c("aggregate.csv", "reading.csv", "range-reading.csv", "paired.csv", "selected-fits.csv", "groups.csv",
            "species.csv", "field.csv", "range.csv", "amplitude.csv", "lattice.csv", "oracle.csv", "long-run-sensitivity.csv"))
  stopifnot(file.copy(file.path(summary, f), out, overwrite = TRUE))
stopifnot(file.copy(file.path(study, "design-statistics.csv"), out, overwrite = TRUE))
stopifnot(file.copy(file.path(study, "audit/audit.csv"), file.path(out, "audit.csv"), overwrite = TRUE))
order_arr <- c("spread", "pairs", "clustered", "grid")
agg <- read.csv(file.path(out, "aggregate.csv")) |> mutate(arrangement = factor(arrangement, order_arr),
  arm = factor(arm, c("oracle", "binary", "two_stage"), c("Oracle: true states, known parameters", "occJSDM: true states", "occJSDM: eDNA survey")),
  reduction = 1 - centred_rmse / zero_field_rmse)
p <- ggplot(agg, aes(arrangement, reduction, colour = arm)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(position = position_dodge(.5), size = 2.5) +
  facet_wrap(~ group) + theme_bw(base_size = 12) +
  labs(x = "Site arrangement, 100 sites each", y = "Field RMSE reduction relative to a flat field",
       colour = NULL, title = "How much of the spatial field each design recovers",
       caption = "Mean of three communities per point. Zero means no better than assuming no field.") +
  theme(legend.position = "bottom")
ggsave(file.path(out, "field-recovery.png"), p, width = 10, height = 6, dpi = 150)
lat <- read.csv(file.path(out, "lattice.csv")) |> filter(bin != "all") |>
  mutate(arrangement = factor(arrangement, order_arr), bin = factor(bin, c("up to 0.02", "0.02 to 0.05", "0.05 to 0.1", "above 0.1"))) |>
  group_by(arrangement, arm, spatial_term, bin) |> summarise(mae = mean(mae), .groups = "drop")
p <- ggplot(lat, aes(bin, 100 * mae, colour = spatial_term, group = spatial_term)) + geom_line() + geom_point() +
  facet_grid(arm ~ arrangement) + theme_bw(base_size = 11) +
  labs(x = "Distance from the nearest surveyed site", y = "Mean absolute prediction error (points)", colour = "Spatial term",
       title = "Prediction at unsurveyed locations by distance from the survey") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(out, "lattice-prediction.png"), p, width = 11, height = 6, dpi = 150)
rng <- read.csv(file.path(out, "range.csv")) |> mutate(arrangement = factor(arrangement, order_arr))
p <- ggplot(rng, aes(arrangement, mass_within_one_step, colour = arm)) + geom_hline(yintercept = .5, linetype = "dashed") +
  geom_point(position = position_dodge(.4), size = 2.5) + theme_bw(base_size = 12) + ylim(0, 1) +
  labs(x = NULL, y = "Posterior mass within one grid step of the true range", colour = "Arm", title = "Range recovery by arrangement")
ggsave(file.path(out, "range-recovery.png"), p, width = 8, height = 5, dpi = 150)
cat("Figures and compact results written to", out, "\n")
```

Run: `Rscript dev/simstudy/spatial-design-sweep/plot.R --repo=. --study=$STUDY`. Open the three PNGs with the Read tool and check axes, labels and that the arrangements appear in the expected order.

- [ ] **Step 5: Write the study README**

Create `dev/simstudy/spatial-design-sweep/README.md` with these sections, one line per paragraph, no pipe tables: a two-sentence statement of the question and the date completed; a results paragraph quoting the reading labels per arrangement and species group, the range-recovery outcome, the oracle-versus-binary and binary-versus-two-stage contrasts, and the lattice result by distance bin, all from `results/`; a paragraph on convergence qualifications naming any fit with remaining flags; a Reproduction section listing the commands from `PLAN.md` in order with `STUDY` and `nohup` as used; a Files section listing each script and each compact result; and a final paragraph recording that the raw archive lives at `dev/simstudy/results/spatial-design-20261001/`, is gitignored, and is needed only for re-auditing draws.

- [ ] **Step 6: Commit**

```bash
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
git add dev/simstudy/spatial-design-sweep/summarise.R dev/simstudy/spatial-design-sweep/verify.R dev/simstudy/spatial-design-sweep/plot.R dev/simstudy/spatial-design-sweep/README.md dev/simstudy/spatial-design-sweep/results
git commit -m "Select, audit, summarise and plot the site-arrangement sweep

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: Teaching bundle, Lesson 2 text, lesson verifier, renders and documentation

**Files:**
- Create: `dev/simstudy/spatial-design-sweep/export-teaching.R`
- Create: `dev/simstudy/spatial-design-sweep/verify-lesson.R`
- Create: `vignettes/teaching-data/spatial-lesson.rds`
- Modify: `vignettes/occJSDM-lesson-2.Rmd` (replace everything from `## 2A.` to the end of the file's old sections), `vignettes/occJSDM-lesson-2.md` (render), `vignettes/teaching-data/lesson-2-*.png` (render), `vignettes/LESSON-PLAN.md`, `README.md`, `TODO.md`, `vignettes/occJSDM.Rmd`, `vignettes/occJSDM.md`.

**Interfaces:**
- Consumes: `STUDY/summary-final/*.csv`, `STUDY/oracle/*.rds`, `STUDY/initial|long/*-result.rds`, `STUDY/inputs/rep01.rds`.
- Produces: bundle `spatial-lesson.rds` with elements `landscape`, `arrangements`, `statistics`, `oracle`, `fits` (groups, species, field, range, amplitude, lattice), `aggregate`, `reading`, `range_reading`, `paired`, `field_maps`, `lattice_maps`, `selected_fits`, `audit`, `provenance`.

- [ ] **Step 1: Write export-teaching.R and export the bundle**

Create `dev/simstudy/spatial-design-sweep/export-teaching.R`:

```r
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) { hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }; substring(hit, nchar(name) + 4L) }
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
summary <- file.path(study, "summary-final"); results <- file.path(repo, "dev/simstudy/spatial-design-sweep/results")
read <- function(f) read.csv(file.path(summary, f))
input1 <- readRDS(file.path(study, "inputs/rep01.rds")); land <- input1$landscape
arrangements <- do.call(rbind, lapply(1:3, function(r) { inp <- readRDS(file.path(study, "inputs", sprintf("rep%02d.rds", r)))
  do.call(rbind, lapply(names(inp$surveys), function(a) data.frame(community = sprintf("rep%02d", r), arrangement = a,
    site = seq_len(100), x = inp$surveys[[a]]$truth$xy[, 1], y = inp$surveys[[a]]$truth$xy[, 2]))) }))
sel <- read("selected-fits.csv"); map_species <- c("species06", "species02")
field_maps <- list(); lattice_maps <- list()
for (i in seq_len(nrow(sel))) {
  r <- readRDS(sel$selected_file[i]); inp <- readRDS(file.path(study, "inputs", paste0(r$job$community, ".rds")))
  for (sp in map_species) { s <- match(sp, inp$landscape$species)
    field_maps[[length(field_maps) + 1L]] <- data.frame(community = r$job$community, arrangement = r$job$arrangement,
      source = r$job$arm, site = seq_len(100), species = sp, value = r$field_median[, s]) }
  if (r$job$community == "rep01") { li <- inp$landscape$index$lattice; s <- match("species06", inp$landscape$species)
    lattice_maps[[length(lattice_maps) + 1L]] <- data.frame(community = "rep01", arrangement = r$job$arrangement,
      source = r$job$arm, cell = seq_along(li), x = inp$landscape$points[li, 1], y = inp$landscape$points[li, 2],
      truth = inp$landscape$psi[li, s], with = r$lattice_mean[, s], without = r$lattice_mean_nospatial[, s],
      distance = r$lattice_distance) }
}
oracle <- read.csv(file.path(study, "oracle/oracle-selected.csv"))
for (i in seq_len(nrow(oracle))) { sp <- sprintf("species%02d", oracle$species[i]); if (!sp %in% map_species) next
  o <- readRDS(oracle$file[i])
  field_maps[[length(field_maps) + 1L]] <- data.frame(community = oracle$community[i], arrangement = oracle$arrangement[i],
    source = "oracle", site = seq_len(100), species = sp, value = o$field_median) }
for (r in 1:3) { inp <- readRDS(file.path(study, "inputs", sprintf("rep%02d.rds", r)))
  for (a in names(inp$surveys)) for (sp in map_species) { s <- match(sp, inp$landscape$species)
    field_maps[[length(field_maps) + 1L]] <- data.frame(community = sprintf("rep%02d", r), arrangement = a, source = "truth",
      site = seq_len(100), species = sp, value = inp$surveys[[a]]$truth$field[, s]) } }
bundle <- list(
  landscape = list(points = land$points, index = land$index, environment = land$environment, field = land$field, psi = land$psi,
                   range = land$range, env_range = land$env_range, species = land$species, prevalence = land$prevalence,
                   B0 = land$B0, B = land$B, cluster_centres = land$cluster_centres),
  arrangements = arrangements, statistics = read.csv(file.path(study, "design-statistics.csv")),
  oracle = read("oracle.csv"), fits = list(groups = read("groups.csv"), species = read("species.csv"), field = read("field.csv"),
    range = read("range.csv"), amplitude = read("amplitude.csv"), lattice = read("lattice.csv")),
  aggregate = read("aggregate.csv"), reading = read("reading.csv"), range_reading = read("range-reading.csv"),
  paired = read("paired.csv"), field_maps = do.call(rbind, field_maps), lattice_maps = do.call(rbind, lattice_maps),
  selected_fits = sel, audit = read.csv(file.path(study, "audit/audit.csv")),
  provenance = list(revision = readLines(file.path(study, "source-revision.txt")),
    compact_hashes = tools::md5sum(list.files(results, full.names = TRUE)),
    script_hashes = tools::md5sum(list.files(file.path(repo, "dev/simstudy/spatial-design-sweep"), pattern = "\\.R$", full.names = TRUE)),
    exported = Sys.time()))
out <- file.path(repo, "vignettes/teaching-data/spatial-lesson.rds")
saveRDS(bundle, out, compress = "xz"); cat("Saved", out, round(file.size(out) / 1e6, 2), "MB\n")
```

Run: `Rscript dev/simstudy/spatial-design-sweep/export-teaching.R --repo=. --study=$STUDY`. Expected: a bundle under 5 MB.

- [ ] **Step 2: Write verify-lesson.R (needs no raw archive) and run it**

Create `dev/simstudy/spatial-design-sweep/verify-lesson.R`:

```r
#!/usr/bin/env Rscript
# Checks the Lesson 2 bundle against the committed compact results. Needs no raw archive.
repo <- normalizePath(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ".")
results <- file.path(repo, "dev/simstudy/spatial-design-sweep/results")
b <- readRDS(file.path(repo, "vignettes/teaching-data/spatial-lesson.rds"))
stopifnot(identical(unname(tools::md5sum(list.files(results, full.names = TRUE))), unname(b$provenance$compact_hashes)))
for (t in c("groups", "species", "field", "range", "amplitude", "lattice"))
  stopifnot(identical(b$fits[[t]], read.csv(file.path(results, paste0(t, ".csv")))))
stopifnot(identical(b$aggregate, read.csv(file.path(results, "aggregate.csv"))),
          identical(b$reading, read.csv(file.path(results, "reading.csv"))))
f <- b$fits$field; f$group <- paste0("prevalence_", f$target * 100, "pct"); f <- f[f$arm == "binary", ]
cells <- aggregate(cbind(centred_rmse, centred_correlation, zero_field_rmse) ~ community + arrangement + group, f, mean)
cells$reduction <- 1 - cells$centred_rmse / cells$zero_field_rmse
for (i in seq_len(nrow(b$reading))) { d <- cells[cells$arrangement == b$reading$arrangement[i] & cells$group == b$reading$group[i], ]
  label <- if (all(d$reduction >= .2 & d$centred_correlation >= .5)) "informative" else
           if (all(d$reduction < .1 | d$centred_correlation < .3)) "uninformative" else "intermediate"
  stopifnot(identical(label, b$reading$label[i])) }
stopifnot(all(abs(colMeans(b$landscape$psi[b$landscape$index$lattice, ]) - b$landscape$prevalence) < 1e-6),
          all(b$audit$passed), nrow(b$selected_fits) == 24L)
cat("Lesson 2 bundle verified against the committed results.\n")
```

Run: `Rscript dev/simstudy/spatial-design-sweep/verify-lesson.R .` Expected: `Lesson 2 bundle verified against the committed results.`

- [ ] **Step 3: Replace the outline sublessons in Lesson 2 with the worked sections**

In `vignettes/occJSDM-lesson-2.Rmd`: change the title to `"Lesson 2: Spatial landscapes and survey design"` in both `title:` and `%\VignetteIndexEntry{}`; add `fig_width: 8` and `fig_height: 5` under `rmarkdown::html_vignette:`; extend the setup chunk to `knitr::opts_chunk$set(echo = TRUE, collapse = FALSE, comment = "#>", fig.width = 8, fig.height = 5, dpi = 150, dev = "png", fig.path = "teaching-data/lesson-2-")` followed by `source("lesson-links.R")`; rewrite the first paragraph of `## What this lesson will add` to `**The concepts come first; the worked sections below test them on a controlled simulation.** The sweep behind them fixed its design before any fit; the protocol and audited results are in the repository's development folder. Dispersal and same-scale environmental confounding are not simulated here; the closing section says what remains.`; keep the second paragraph of that section and the whole `## Spatial effects, inference and sampling design` section unchanged; then replace everything from `## 2A. Smooth environmental gradients shape species distributions` to the end of the file with the following.

````markdown
```{r lesson-2-load, message=FALSE}
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

sweep <- readRDS("teaching-data/spatial-lesson.rds")

arrangement_order <- c("spread", "pairs", "clustered", "grid")
arrangement_labels <- c(spread = "Spread at random", pairs = "Spread plus close pairs",
                        clustered = "Ten clusters of ten", grid = "Regular grid (control)")
arm_labels <- c(oracle = "Oracle: true states, known parameters",
                binary = "occJSDM: true states", two_stage = "occJSDM: eDNA survey")
group_labels <- c(prevalence_5pct = "5% species", prevalence_25pct = "25% species",
                  prevalence_75pct = "75% species")

theme_set(theme_bw(base_size = 12))
```

## 2A. One landscape, four surveys

Everything in the worked sections comes from a simulation designed to isolate one question: with the budget fixed at 100 sites in a fixed study area, how does the arrangement of those sites change what the spatial submodel can learn? The landscape has a single broad environmental gradient, which every arrangement estimates equally well, and one independent spatial field per species with a short range, 3% of the side of the area, on which the arrangements differ. All species share that range and a field standard deviation of 1 on the log-odds scale, which is what occJSDM assumes, so the sweep tests information, not model mismatch. Three independent communities were generated; the maps below show the first.

```{r landscape-maps, fig.height=4.2}
lat <- sweep$landscape$index$lattice
landscape_cells <- tibble(
  x = sweep$landscape$points[lat, 1], y = sweep$landscape$points[lat, 2],
  environment = sweep$landscape$environment[lat],
  field = sweep$landscape$field[lat, "species06"],
  probability = sweep$landscape$psi[lat, "species06"]
) |>
  pivot_longer(-c(x, y), names_to = "layer", values_to = "value") |>
  mutate(layer = factor(layer, c("environment", "field", "probability"),
                        c("Environmental gradient", "Spatial field (species 6)", "True occupancy probability (species 6)")))

ggplot(landscape_cells, aes(x, y, fill = value)) +
  geom_raster() +
  facet_wrap(~ layer, scales = "free") +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = NULL, title = "One simulated landscape, community 1")
```

The gradient alone would give species 6 a smooth trend across the area. The field adds patches a few percent of the side wide, and the probability map is their sum on the log-odds scale. A geographically structured distribution is therefore expected even where no spatial process acts, and the spatial submodel's job is only the patches.

```{r arrangement-maps, fig.height=3.2}
sites <- sweep$arrangements |>
  filter(community == "rep01") |>
  mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels))

ggplot(sites, aes(x, y)) +
  geom_point(size = 1.1, colour = "#0072B2") +
  facet_wrap(~ arrangement, nrow = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "Four ways to place 100 sites")
```

The design table records what each arrangement gives the model. A site can inform the field only if another site lies within about one range of it, where their field values are correlated above 0.5.

```{r design-table}
sweep$statistics |>
  group_by(arrangement) |>
  summarise(nearest_neighbour = mean(mean_nearest_neighbour),
            with_correlated_neighbour = mean(fraction_with_half_neighbour),
            standardised_range = mean((standardised_range_x + standardised_range_y) / 2),
            axis_ratio = mean(axis_sd_ratio), .groups = "drop") |>
  mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels)) |>
  arrange(arrangement) |>
  knitr::kable(digits = c(0, 3, 2, 3, 2),
               col.names = c("Arrangement", "Mean nearest-neighbour distance", "Fraction of sites with a correlated neighbour",
                             "Range on the fitter's standardised scale", "Ratio of axis spreads"),
               caption = "Means over three communities; the study area has side 1 and the field range is 0.03")
```

One trap is worth naming. occJSDM standardises each coordinate axis before fitting, so shrinking the whole study area changes nothing: the standardised range column is about 0.10 for every arrangement, inside the fitter's grid of 0.01 to 0.30. Closer spacing means more sites within one range of each other, which at a fixed budget means clustering some of them.

## 2B. What the survey data contain

Before asking what occJSDM recovers, ask what the data allow. An oracle sampler was handed the true occupied states, intercept, slope, range and amplitude and asked only for the field. Nothing can do better from the same states. Its error is compared with the error of assuming a flat field.

```{r oracle-recovery, fig.height=4.8}
recovery <- sweep$aggregate |>
  mutate(reduction = 1 - centred_rmse / zero_field_rmse,
         arrangement = factor(arrangement, arrangement_order, arrangement_labels),
         arm = factor(arm, names(arm_labels), arm_labels),
         group = factor(group, names(group_labels), group_labels))

ggplot(filter(recovery, arm == arm_labels["oracle"]), aes(arrangement, reduction)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = 1 - centred_rmse_max / zero_field_rmse, ymax = 1 - centred_rmse_min / zero_field_rmse),
                  colour = "#0072B2") +
  facet_wrap(~ group) +
  labs(x = NULL, y = "Field error reduction relative to a flat field",
       title = "The ceiling: what the occupancy states contain about the field",
       caption = "Point: mean of three communities. Bar: their range. Zero means the states say nothing about the field.") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
```

```{r reading-labels}
sweep$reading |>
  mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels),
         group = factor(group, names(group_labels), group_labels)) |>
  arrange(group, arrangement) |>
  select(group, arrangement, label, min_reduction, min_correlation) |>
  knitr::kable(digits = 2, col.names = c("Species group", "Arrangement", "Reading", "Smallest error reduction", "Smallest correlation"),
               caption = "Reading rules fixed before fitting, applied to occJSDM's true-state fits across all three communities")
```

The reading labels were defined before any result existed: informative means at least a 20% error reduction and a correlation of at least 0.5 in every community, uninformative means under 10% or under 0.3 in every community, and anything else is intermediate. Read the table rather than the prose for the result; the prose below describes the pattern the table shows.

`r if (all(sweep$reading$label[sweep$reading$arrangement == "spread"] != "informative")) "With sites spread at random, the field is not recoverable for any species group: most sites have no correlated neighbour, so one occupancy state per site tells the model nothing about the patches." else "Spread sites are informative for at least one species group in this landscape, which was not expected from the September studies; the table and the maps below show how far."` `r if (any(sweep$reading$label[sweep$reading$arrangement %in% c("clustered", "pairs")] == "informative")) "Placing sites within one range of each other changes that for the common species. Rarity remains a separate limit: a species at 5% occupancy has about five occupied sites among 100, too few to reveal where its patches are, however the sites are arranged." else "Even the clustered designs do not reach the informative threshold for any group in all three communities; the intermediate cells show how much is recovered."`

## 2C. What occJSDM delivers

The full model must also estimate the intercepts, slopes, range and amplitude, and in the survey arm it must see the field through two field samples, two primers and six PCRs per sample. The oracle-to-true-state gap is the cost of estimation; the true-state-to-survey gap is the cost of detection.

```{r fit-recovery, fig.height=5}
ggplot(recovery, aes(arrangement, reduction, colour = arm)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(position = position_dodge(.5), size = 2.4) +
  facet_wrap(~ group) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#7B3294")) +
  labs(x = NULL, y = "Field error reduction relative to a flat field", colour = NULL,
       title = "Oracle ceiling, true-state fit and eDNA-survey fit",
       caption = "Means of three communities. Each point's community range is in the saved tables.") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
```

```{r field-maps, fig.height=6.5}
maps <- sweep$field_maps |>
  filter(community == "rep01", species == "species06", arrangement %in% c("spread", "clustered")) |>
  left_join(filter(sweep$arrangements, community == "rep01"), by = c("community", "arrangement", "site")) |>
  mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels),
         source = factor(source, c("truth", "oracle", "binary", "two_stage"),
                         c("Truth", arm_labels)))

ggplot(maps, aes(x, y, colour = value)) +
  geom_point(size = 1.8) +
  facet_grid(arrangement ~ source) +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, colour = "Field (log-odds)",
       title = "Species 6, community 1: the field at the surveyed sites",
       caption = "Fitted values are posterior medians. Compare the fitted maps with the truth column in each row.")
```

```{r range-amplitude, fig.height=4.2}
sweep$fits$range |>
  mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels),
         arm = factor(arm, names(arm_labels)[-1], arm_labels[-1])) |>
  ggplot(aes(arrangement, mass_within_one_step, colour = arm)) +
  geom_hline(yintercept = .5, linetype = "dashed") +
  geom_point(position = position_dodge(.4), size = 2.4) +
  scale_colour_manual(values = c("#D55E00", "#7B3294")) +
  ylim(0, 1) +
  labs(x = NULL, y = "Posterior mass within one grid step of the true range", colour = NULL,
       title = "Does the fit find the range?", caption = "One point per community. Dashed line: the reading rule.") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
```

```{r occupancy-bands}
sweep$fits$groups |>
  filter(metric == "occupancy", group %in% c("low", "medium", "high")) |>
  group_by(arrangement, arm, group) |>
  summarise(signed = 100 * mean(bias), absolute = 100 * mean(mae), .groups = "drop") |>
  mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels),
         arm = factor(arm, names(arm_labels)[-1], arm_labels[-1]),
         group = factor(group, c("low", "medium", "high"), c("Below 20%", "20% to 80%", "Above 80%"))) |>
  arrange(arm, arrangement, group) |>
  knitr::kable(digits = 1, col.names = c("Arrangement", "Arm", "True probability", "Signed error (points)", "Absolute error (points)"),
               caption = "Occupancy error at the surveyed sites, mean of three communities")
```

Convergence qualifications are part of the result. `r sum(sweep$selected_fits$needs_long)` of the 24 initial fits met the prespecified rule for a longer run; `r sum(nzchar(sweep$selected_fits$selected_reasons))` selected fits retain a flag after it. Every flagged fit is kept and named in the saved tables.

## 2D. Predicting unsurveyed locations

Clustering buys neighbours and spends coverage. The lattice of 1,600 unsurveyed locations shows the trade. Prediction error is plotted against the distance from each location to its nearest surveyed site, with and without the spatial term.

```{r lattice-prediction, fig.height=6}
sweep$fits$lattice |>
  filter(bin != "all") |>
  group_by(arrangement, arm, spatial_term, bin) |>
  summarise(mae = 100 * mean(mae), .groups = "drop") |>
  mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels),
         arm = factor(arm, names(arm_labels)[-1], arm_labels[-1]),
         bin = factor(bin, c("up to 0.02", "0.02 to 0.05", "0.05 to 0.1", "above 0.1")),
         spatial_term = factor(spatial_term, c("with", "without"), c("Environment and field", "Environment only"))) |>
  ggplot(aes(bin, mae, colour = spatial_term, group = spatial_term)) +
  geom_line() + geom_point() +
  facet_grid(arm ~ arrangement) +
  scale_colour_manual(values = c("#0072B2", "grey45")) +
  labs(x = "Distance from the nearest surveyed site", y = "Mean absolute error (points)", colour = NULL,
       title = "Prediction error at unsurveyed locations",
       caption = "Where the field was learned, it helps only near surveyed sites. Far from any site both lines meet: the field is unknown there.") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
```

```{r lattice-maps, fig.height=4.5}
sweep$lattice_maps |>
  filter(source == "binary", arrangement %in% c("spread", "clustered")) |>
  mutate(error = 100 * (with - truth),
         arrangement = factor(arrangement, arrangement_order, arrangement_labels)) |>
  ggplot(aes(x, y, fill = error)) +
  geom_raster() +
  geom_point(data = filter(sweep$arrangements, community == "rep01", arrangement %in% c("spread", "clustered")) |>
               mutate(arrangement = factor(arrangement, arrangement_order, arrangement_labels)),
             aes(x, y), inherit.aes = FALSE, size = .6, colour = "black") +
  facet_wrap(~ arrangement) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Error (points)",
       title = "Species 6, community 1: prediction error over the area from true-state fits",
       caption = "Black dots are the surveyed sites.")
```

This is the design decision the conceptual section anticipated. Spread sites give an even but featureless map; clusters give a sharp local map and no information between clusters; close pairs added to a spread design are the hedge, and the saved tables show where they land between the two. The right choice depends on whether the map must be local or regional, which is a question about the study's purpose, not about the model.

## What this establishes, and what it does not

The sweep is a controlled, model-matched simulation: one broad gradient, one field range shared by all species, no dispersal, no species-specific ranges, and no environmental covariate at the field's scale. Within that, three communities support the reading labels above, not confidence intervals. Rare species remain unrecoverable at every arrangement. The two-stage survey arm shows what detection error costs on top of estimation. Nothing here validates spatial prediction on real data, establishes interval coverage, or shows what happens when species disperse at different scales; those are the separate contrasts still to come, with same-scale environmental confounding the first candidate.

Try these with the saved tables, without refitting:

1. From `sweep$statistics`, compute the fraction of sites with a correlated neighbour for a design of your own, by replacing the coordinates with yours and the range with a plausible one for your system.
2. `sweep$lattice_maps` has every lattice cell for community 1 with its distance to the nearest site. Recompute the distance bins at 0.01, 0.03 and 0.06 and redraw the prediction-error figure for the two arrangements it contains.
3. Pick the species group and arrangement whose reading is intermediate and explain, from the three communities' values, what made the rule withhold the informative label.

## Reproduction record

The protocol, scripts, compact results, audit and figures are under `dev/simstudy/spatial-design-sweep/` in the source repository, with the frozen production revision recorded in the protocol's amendments. The compact bundle `teaching-data/spatial-lesson.rds` carries everything this lesson renders; `verify-lesson.R` checks it against the committed results without the raw archive. No fit is rerun while knitting.

Return to the [Quickstart and lesson guide](occJSDM.md), [Lesson 1](occJSDM-lesson-1.md) or [Lesson 3](occJSDM-lesson-3.md).
````

The old `## What will stay consistent across the lessons` section is removed; its content is now covered by "What this establishes". After editing, run `grep -n "—" vignettes/occJSDM-lesson-2.Rmd` and expect no output.

- [ ] **Step 4: Render Lesson 2, run the link check, inspect the figures**

Run:

```sh
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-2.Rmd", quiet=TRUE); rmarkdown::render("vignettes/occJSDM-lesson-2.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE), quiet=TRUE)'
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
```

Expected: both renders succeed; `Shared lesson-link regression checks passed.` Open each `vignettes/teaching-data/lesson-2-*.png` with the Read tool and check the figures: raster maps fill the square, the four arrangements appear in order, the recovery figure has three facets, the lattice figure has eight panels. Fix any figure problem in the Rmd, not in the bundle.

- [ ] **Step 5: Update the lesson plan, README, TODO and Quickstart**

- In `vignettes/LESSON-PLAN.md`: change the Lesson 2 status cell to `Built: concept section plus the worked site-arrangement sweep; dispersal and confounding deferred`; in the `### Lesson 2: space, habitat and dispersal` section, replace the four numbered items with one paragraph saying the sweep replaced 2A to 2D on 1 October 2026, naming the arrangements, the oracle, the two arms and the reading rules, and listing dispersal, species-specific ranges and same-scale confounding as the deferred contrasts; append a decisions-log entry dated with the completion date naming the protocol path and the headline reading labels.
- In `README.md`: replace the Lesson 2 bullet's second sentence with `Its worked sections test four arrangements of 100 sites, spread, spread with close pairs, clustered, and a grid control, against an oracle ceiling and occJSDM fits on true states and on eDNA surveys, and show the coverage cost of clustering when predicting unsurveyed locations.`
- In `TODO.md`: strike the Lesson 2 study-design item with `~~` and append ` **DONE <date>:** see [the sweep](dev/simstudy/spatial-design-sweep/README.md) and Lesson 2.` without changing other items.
- In `vignettes/occJSDM.Rmd`: change the Lesson 2 bullet to `[Lesson 2: Spatial landscapes and survey design](occJSDM-lesson-2.md). What the spatial field learns, and a worked sweep of four ways to place 100 sites, with an oracle ceiling, true-state and eDNA-survey fits, and prediction at unsurveyed locations.` Then re-render the Quickstart Markdown: `Rscript -e 'rmarkdown::render("vignettes/occJSDM.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE), quiet=TRUE)'`.

- [ ] **Step 6: Run the full checks and commit**

Run:

```sh
cd ~/src/occJSDM-worktrees/lesson-2-spatial-sweep
Rscript dev/simstudy/spatial-design-sweep/verify-lesson.R .
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
git grep -n "—" -- vignettes dev/simstudy/spatial-design-sweep ':!*.png' ':!*.rds' | wc -l
git status --short
```

Expected: verifier passes, link check passes, zero em-dashes, and the status lists only the intended files. Then:

```bash
git add vignettes/occJSDM-lesson-2.Rmd vignettes/occJSDM-lesson-2.md vignettes/teaching-data/spatial-lesson.rds vignettes/teaching-data/lesson-2-*.png vignettes/LESSON-PLAN.md README.md TODO.md vignettes/occJSDM.Rmd vignettes/occJSDM.md dev/simstudy/spatial-design-sweep/export-teaching.R dev/simstudy/spatial-design-sweep/verify-lesson.R
git commit -m "Build Lesson 2 from the site-arrangement sweep

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

- [ ] **Step 7: Open the pull request**

Push the branch and open a PR to main with `gh pr create`, titled `Lesson 2: site-arrangement sweep and worked spatial sections`, whose body summarises the protocol, the reading labels per arrangement and group, the range-recovery outcome, the lattice trade-off, the convergence qualifications, and the deferred contrasts, and ends with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`. Report the PR link to Doug.

---

## Self-review notes

- Spec coverage: simulation (Task 1), oracle (Task 2), freeze and full fits with the exact settings and schedules (Task 3), five outcomes and reading rules (Tasks 3 and 4), validation tests and audit (Tasks 1 to 4), compute budget and detached launch (Tasks 2 to 4), lesson mapping and bundle with the archive-free verifier (Task 5), reproduction commands (Task 4 README and Task 5), decision log and amendments (Task 1 and Task 3 freeze).
- The oracle's lattice prediction (outcome 4 for the oracle) is computed in Task 2 and stored in each oracle result; `summarise.R` does not tabulate it, and the lesson's lattice figure uses only the two fitted arms. If the oracle line is wanted there, add an `oracle_lattice.csv` in `summarise.R` from each result's `lattice$probability_mean` against `psi` at the lattice, binned with `distance_bins()`.
- Type consistency checked: `score_sweep_fit` returns `range_summary`, which `collect(..., "range_summary")` reads into the `range` table; `field` rows carry `target`, which `summarise.R` and `verify-lesson.R` turn into `group`; `selected_fits` columns `needs_long` and `selected_reasons` are what the lesson's inline R reads; the score test expects the lattice table's `bin` values including `all`.
