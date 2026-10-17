# Refit replicates of named arms, optionally changed by a variant, and keep per-species summaries of the spread of
# each posterior draw's occupancy linear predictor around B0, for README.md.
# Usage: Rscript refit_spread.R <package root> <output directory> <cores> <jobs> <replicates> <nburn> <niter>
#   jobs: comma-separated labels, each an arm or arm+variant (variants: sb0, d4, gt4); replicates: first-last.
# Each job saves <output directory>/<job>-<replicate>.rds as it finishes, so a killed run keeps its finished jobs;
# a rerun skips those. Workers are set up as dev/simstudy/run_study.R sets up its own: clusterSetRNGStream(cl,
# 20260727L), one RcppParallel thread, then devtools::load_all(). Under that L'Ecuyer generator the set.seed() in
# draw_truth() gives the runner's datasets and chain seeds, so at the study's chain length a refit without a variant
# reproduces its arm's stored rows exactly. A variant keeps its arm's seed, so it fits the same dataset.
args <- commandArgs(trailingOnly = TRUE)
pkg_root <- normalizePath(args[1], mustWork = TRUE)
out_dir <- args[2]
cores <- as.integer(args[3])
labels <- strsplit(args[4], ",", fixed = TRUE)[[1]]
reps <- as.integer(strsplit(args[5], "-", fixed = TRUE)[[1]])
mcmc <- list(nchain = 2, nburn = as.integer(args[6]), niter = as.integer(args[7]), nthin = 1)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
library(parallel)

jobs <- expand.grid(replicate = seq(reps[1], reps[2]), job = labels, stringsAsFactors = FALSE)
jobs$file <- file.path(out_dir, sprintf("%s-%02d.rds", jobs$job, jobs$replicate))
jobs <- jobs[!file.exists(jobs$file), ]
message(sprintf("%s: %d jobs to run at %d + %d", format(Sys.time()), nrow(jobs), mcmc$nburn, mcmc$niter))

cl <- makeCluster(cores, outfile = "")
clusterSetRNGStream(cl, 20260727L)
clusterExport(cl, c("pkg_root"))
invisible(clusterEvalQ(cl, {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
  RcppParallel::setThreadOptions(numThreads = 1)
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))
  NULL
}))

# What each variant changes. fit_d goes on the scenario, so simstudy_fit_args() and the rows see it; listPriors is
# merged into the scenario's; listParams is merged into what simstudy_fit_args() hands runOccJSDM(), since the
# scenario's gt also sets the simulated latent traits.
variants <- list(
  sb0 = list(listPriors = list(sigma_b0 = 2)),
  d4 = list(fit_d = 4L),
  gt4 = list(listParams = list(n_lattrait = 4L)))

# The intercept a, one per column of the centred spreads rc (n x S), at which the site mean of plogis(a + rc) equals
# target. Newton steps on the logit of that mean, which is close to linear in a; uniroot for any column that fails.
solve_intercept <- function(rc, target) {
  lt <- stats::qlogis(target)
  a <- lt
  g <- rep(Inf, length(a))
  for (k in 1:60) {
    p <- stats::plogis(sweep(rc, 2, a, "+"))
    f <- colMeans(p)
    g <- stats::qlogis(f) - lt
    if (all(is.finite(g)) && max(abs(g)) < 1e-10) break
    dg <- colMeans(p * (1 - p)) / (f * (1 - f))
    step <- g / dg
    step[!is.finite(step)] <- 0
    a <- a - pmax(pmin(step, 5), -5)
  }
  bad <- which(!is.finite(g) | abs(g) >= 1e-8)
  for (s in bad) {
    a[s] <- stats::uniroot(function(x) mean(stats::plogis(x + rc[, s])) - target[s], c(-60, 60),
                           tol = 1e-12)$root
  }
  a
}

col_sd <- function(x) sqrt(pmax(colSums(sweep(x, 2, colMeans(x))^2) / (nrow(x) - 1), 0))
col_cov <- function(x, y) colSums(sweep(x, 2, colMeans(x)) * sweep(y, 2, colMeans(y))) / (nrow(x) - 1)

# Per-draw statistics of the fitted linear predictor, summarised per species, against the truth.
spread_summaries <- function(fit, sim, a) {
  jo <- fit$results_output$jsdm_output
  inf <- fit$infos
  n <- inf$n; S <- inf$S; ps <- inf$ps
  X <- as.matrix(process_covariates(inf$data_info, a$occCovariates, "Site", n, remove_intercept = TRUE,
                                    spline_vars = FALSE)$df)
  Ks_all <- NULL
  if (ps > 0) invisible(utils::capture.output(
    Ks_all <- suppressMessages(precomputeSORmatrices(inf$l_s_grid, inf$list_Xs))$Ks_all))

  jp <- sim$true_params$jsdmParams_true
  B0_t <- jp$B0
  UL_t <- jp$U %*% jp$L
  F_t <- jp$spatField
  r_t <- jp$eta - matrix(B0_t, n, S, byrow = TRUE)
  XB_t <- r_t - UL_t - F_t
  site_t <- UL_t + F_t
  m_t <- colMeans(r_t)
  sd_t <- col_sd(r_t)
  rc_t <- sweep(r_t, 2, m_t)
  psibar_t <- colMeans(stats::plogis(jp$eta))
  truth <- data.frame(
    species = seq_len(S), B0_true = B0_t, m_true = m_t, a_true = B0_t + m_t, sd_r_true = sd_t,
    sd_XB_true = col_sd(XB_t), sd_UL_true = col_sd(UL_t), sd_SP_true = col_sd(F_t), sd_site_true = col_sd(site_t),
    var_XB_true = col_sd(XB_t)^2, var_site_true = col_sd(site_t)^2, cov_XB_site_true = col_cov(XB_t, site_t),
    m_UL_true = colMeans(UL_t), m_SP_true = colMeans(F_t), psibar_true = psibar_t,
    zbar_true = colMeans(sim$true_params$z_true))

  niter <- dim(jo$B0_output)[2]; nchain <- dim(jo$B0_output)[3]
  stat_names <- c("B0", "m_r", "a_fit", "sd_r", "sd_XB", "sd_UL", "sd_SP", "sd_site", "var_XB", "var_site",
                  "cov_XB_site", "m_UL", "m_SP", "psibar", "aJ", "aS", "cor_r")
  draws <- array(NA_real_, c(niter, nchain, S, length(stat_names)), dimnames = list(NULL, NULL, NULL, stat_names))
  sum_r <- sum_site <- sum_XB <- sum_psi <- matrix(0, n, S)
  for (ch in seq_len(nchain)) for (t in seq_len(niter)) {
    B0 <- jo$B0_output[, t, ch]
    XB <- X %*% matrix(jo$B_output[, , t, ch], ncol = S)
    UL <- matrix(jo$U_output[, , t, ch], n) %*% matrix(jo$L_output[, , t, ch], ncol = S)
    SP <- if (ps > 0) {
      Ks <- matrix(Ks_all[, , jo$idx_ls_output[t, ch]], nrow = n)
      KsBproduct(Ks, matrix(jo$Bs_output[, , t, ch], ncol = S), inf$list_Xs$Xs_centers)
    } else matrix(0, n, S)
    site <- UL + SP
    r <- XB + site
    psi <- stats::plogis(sweep(r, 2, B0, "+"))
    m_r <- colMeans(r)
    sd_r <- col_sd(r)
    sum_r <- sum_r + r; sum_site <- sum_site + site; sum_XB <- sum_XB + XB; sum_psi <- sum_psi + psi
    draws[t, ch, , ] <- cbind(
      B0, m_r, B0 + m_r, sd_r, col_sd(XB), col_sd(UL), col_sd(SP), col_sd(site), col_sd(XB)^2, col_sd(site)^2,
      col_cov(XB, site), colMeans(UL), colMeans(SP), colMeans(psi),
      solve_intercept(sweep(r, 2, m_r), psibar_t),
      solve_intercept(sweep(rc_t, 2, sd_r / sd_t, "*"), psibar_t),
      col_cov(r, r_t) / (sd_r * sd_t))
  }
  nd <- niter * nchain
  rebuild_error <- max(abs(sum_psi / nd - fit$results_output$psi_output))

  summ <- truth
  for (k in stat_names) {
    x <- matrix(draws[, , , k], nd, S)
    summ[[paste0(k, "_mean")]] <- colMeans(x)
    summ[[paste0(k, "_sd")]] <- col_sd(x)
    summ[[paste0(k, "_lo")]] <- apply(x, 2, stats::quantile, 0.025, names = FALSE)
    summ[[paste0(k, "_hi")]] <- apply(x, 2, stats::quantile, 0.975, names = FALSE)
  }
  for (k in c("B0", "sd_r", "aJ", "m_r")) {
    summ[[paste0(k, "_rhat")]] <- vapply(seq_len(S), function(s) posterior::rhat(draws[, , s, k]), 0)
    summ[[paste0(k, "_ess")]] <- vapply(seq_len(S), function(s) posterior::ess_bulk(draws[, , s, k]), 0)
  }

  hyper_one <- function(x) if (is.null(x)) NULL else
    c(mean = mean(x), lo = stats::quantile(x, 0.025, names = FALSE), hi = stats::quantile(x, 0.975, names = FALSE),
      rhat = posterior::rhat(x), ess = posterior::ess_bulk(x))
  hyper <- list(
    sigma_b = hyper_one(jo$sigmab_output), sigma_bs = hyper_one(jo$sigmabs_output),
    sigma_h = hyper_one(jo$sigmah_output),
    l_s = hyper_one(matrix(inf$l_s_grid[jo$idx_ls_output], niter, nchain)),
    mu_b0 = hyper_one(jo$mu_b0_output))
  hyper <- hyper[!vapply(hyper, is.null, TRUE)]

  list(species = summ, hyper = hyper, rebuild_error = rebuild_error,
       site_means = list(r = sum_r / nd, site = sum_site / nd, XB = sum_XB / nd),
       site_truth = list(r = r_t, site = site_t, XB = XB_t))
}

one <- function(k) {
  job <- jobs[k, ]
  parts <- strsplit(job$job, "+", fixed = TRUE)[[1]]
  arm <- parts[1]
  variant <- if (length(parts) > 1) parts[2] else ""
  sc <- simstudy_scenarios()
  s <- sc[[which(vapply(sc, `[[`, "", "label") == arm)]]
  v <- if (nzchar(variant)) variants[[variant]] else list()
  if (nzchar(variant) && is.null(v)) stop("unknown variant ", variant)
  started <- proc.time()[["elapsed"]]
  truth <- draw_truth(s, simstudy_seed_for(s, job$replicate))
  sim <- simstudy_simulate(truth)

  s_fit <- s
  if (!is.null(v$fit_d)) s_fit$fit_d <- v$fit_d
  if (!is.null(v$listPriors)) s_fit$listPriors <- utils::modifyList(s$listPriors %||% list(), v$listPriors)
  a <- simstudy_fit_args(sim, s_fit)
  if (!is.null(v$listParams)) a$listParams <- utils::modifyList(a$listParams, v$listParams)
  # As simstudy_fit(), with the variant's arguments.
  fit <- suppressMessages(suppressWarnings(
    runOccJSDM(a$data, listParams = a$listParams, listPriors = a$listPriors, occCovariates = a$occCovariates,
               collCovariates = a$collCovariates, spatCovariates = a$spatCovariates, MCMCparams = mcmc)))

  s_rows <- s_fit
  s_rows$label <- job$job
  sp <- spread_summaries(fit, sim, a)
  res <- list(
    job = job$job, arm = arm, variant = variant, replicate = job$replicate, mcmc = mcmc,
    listParams = a$listParams, listPriors = a$listPriors,
    intercept_prior = fit$infos$intercept_prior,
    rows = simstudy_rows(fit, sim, truth, s_rows, job$replicate),
    occupancy = simstudy_occupancy_rows(fit, sim, s_rows, job$replicate),
    species = sp$species, hyper = sp$hyper, rebuild_error = sp$rebuild_error,
    site_means = sp$site_means, site_truth = sp$site_truth,
    elapsed = proc.time()[["elapsed"]] - started)
  saveRDS(res, job$file)
  msg <- sprintf("%s %s replicate %d in %.0f s (rebuild error %.1e)", format(Sys.time(), "%H:%M:%S"),
                 job$job, job$replicate, res$elapsed, res$rebuild_error)
  message(msg)
  msg
}
clusterExport(cl, c("jobs", "one", "mcmc", "variants", "solve_intercept", "col_sd", "col_cov",
                    "spread_summaries"))
done <- parLapplyLB(cl, seq_len(nrow(jobs)), function(k) {
  tryCatch(one(k), error = function(e) paste("FAILED", jobs$job[k], jobs$replicate[k], conditionMessage(e)))
})
message(paste(unlist(done), collapse = "\n"))
message(sprintf("%s: finished, %d files in %s", format(Sys.time()),
                length(list.files(out_dir, "[.]rds$")), out_dir))
stopCluster(cl)
