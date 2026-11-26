# Tests Alex Diana's suggestion on PR 46 (25 November 2026): that sigma_b falls
# near zero because 4 latent traits over-fit the 5 measured covariates (here 6
# design columns, since one covariate is categorical with 3 levels), and
# that fewer latent traits would stop it. The five simulated data sets and
# seeds of repro.R are fitted with a library that has no sigma_b move (built
# from a2a64a1), once with 2 latent traits (label "gt2") and once with 4
# (label "gt4", the control, which should repeat the "base" files of repro.R).
# Every data set comes from simulateOccJSDMData(); no dataset content is read.
# The truth is simulated with 4 latent traits and 3 embedding blends in both
# arms, so the 2-trait fit has fewer latent traits than the truth.
#
# Usage:
#   Rscript --vanilla lattrait.R --lib PATH --lattrait 2|4 --datasets 1:5 --workers 5
#     [--commit SHA] [--nburn 5000] [--niter 20000] [--outdir DIR]
# --nburn, --niter and --outdir exist for the short check of the script only.
#
# Writes dev/simstudy/results/sigma-b-move-20261122/gt<N>-d<k>.rds per data set.
#
# Reading rules, fixed on 26 November 2026 before any run:
# 1. Low stays. A kept sweep is low when sigma_b is under 0.02, and a stay is
#    long at 500 sweeps or more, as in summarise.R. With 2 latent traits the
#    reading is "long stays gone" if no chain of the five fits has a long stay,
#    and "long stays remain" otherwise. It bears on the suggestion only if the
#    4-trait arm has at least one long stay on the same seeds.
# 2. Measured coefficients. A fit fails if any measured coefficient has an
#    R-hat over 1.1 at the first 5,000 kept sweeps. Counts are given per arm.
# 3. Recovery. Per data set, the root mean squared error of the posterior mean
#    with 2 latent traits over that with 4, for (a) the embedding's
#    contribution to the linear predictor and (b) the measured coefficients.
#    Each reads "no harm" if the ratio is at most 1.10 in at least 4 of the 5
#    data sets, and "harm" otherwise. The coverage of the measured
#    coefficients' 95% intervals is given beside it, with no rule.

Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing --", name)
    return(default)
  }
  args[i + 1L]
}
lib <- normalizePath(arg("lib"), mustWork = TRUE)
gt_fit <- as.integer(arg("lattrait"))
if (!gt_fit %in% c(2L, 4L)) stop("--lattrait must be 2 or 4")
label <- paste0("gt", gt_fit)
# "1:5" or "1,3,5"; parsed without evaluating code.
parse_datasets <- function(x) {
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    ends <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq(ends[1], ends[2]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}
datasets <- parse_datasets(arg("datasets"))
workers <- as.integer(arg("workers"))
nburn <- as.integer(arg("nburn", "5000"))
niter <- as.integer(arg("niter", "20000"))

library(occJSDM, lib.loc = lib)
pkg_path <- normalizePath(dirname(system.file(package = "occJSDM")))
if (pkg_path != lib) stop("occJSDM was loaded from ", pkg_path, ", not from ", lib)
RcppParallel::setThreadOptions(numThreads = 1)

# The study directory, found from this script's own path.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
study_dir <- dirname(normalizePath(script))
out_dir <- arg("outdir", file.path(study_dir, "..", "results", "sigma-b-move-20261122"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)

# The commit the library was built from, as given on the command line.
commit <- arg("commit", NA_character_)

now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
power_state <- function() {
  paste(c(system("pmset -g batt", intern = TRUE)[1],
          grep("lowpowermode", system("pmset -g", intern = TRUE), value = TRUE)),
        collapse = "; ")
}

S <- 48L; n <- 256L
datasettings <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L,
                     K = rep(1L, n), ncov_psi = 5L, ncov_theta = 1L)
jsdmParams <- list(gt = 4L, d = 2L, ds = 0, sigma_b = 0.1, sigma_bs = 0.5, sigma_ts = 0.5,
                   sigma_h = 1, sigma_s = 0.5, l_s = 0.3, tau = rep(1, S),
                   useSpatField = FALSE, embedding = list(n_emb = 64L))
MCMCparams <- list(nchain = 4L, nburn = nburn, niter = niter, nthin = 1L)

# The largest unsplit coda R-hat among the measured rows of B_output, over all
# chains, and the number of measured coefficients over 1.1. As in repro.R.
measured_rhat <- function(B_output, is_emb, keep) {
  meas <- which(!is_emb)
  nchain <- dim(B_output)[4]
  chains <- lapply(seq_len(nchain), function(ch) {
    draws <- array(B_output[meas, , keep, ch, drop = FALSE],
                   dim = c(length(meas), dim(B_output)[2], length(keep)))
    m <- matrix(aperm(draws, c(3, 1, 2)), nrow = length(keep))   # iterations by coefficients
    coda::mcmc(m)
  })
  r <- coda::gelman.diag(coda::mcmc.list(chains), autoburnin = FALSE,
                         multivariate = FALSE)$psrf[, 1]
  c(max = max(r, na.rm = TRUE), over_1.1 = sum(r > 1.1, na.rm = TRUE))
}

fit_one <- function(k) {
  set.seed(1000L + k)
  sim <- simulateOccJSDMData(datasettings, list(), jsdmParams, model = "binary")
  dat <- sim$data_list
  emb <- dat$embeddingVars
  meas <- setdiff(grep("^X_psi", names(dat$info), value = TRUE), emb)
  occ <- c(meas, emb)
  tp <- sim$true_params$jsdmParams_true
  is_emb_true <- tp$covariate_group == "embedding"
  if (sum(is_emb_true) != length(emb) || nrow(tp$B) != length(tp$covariate_group) ||
      ncol(tp$B) != S)
    stop("the true B does not have one embedding row per embedding column")
  B_emb_true <- tp$B[is_emb_true, , drop = FALSE]
  B_meas_true <- tp$B[!is_emb_true, , drop = FALSE]

  power_start <- power_state()
  started <- now()
  cat(sprintf("%s-d%d started %s; power %s\n", label, k, started, power_start))
  set.seed(2000L + k)
  occJSDM:::setOccJSDMSeed(2000L + k)
  warns <- character()
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    suppressMessages(runOccJSDM(dat,
      listParams = list(n_factors = 2L, n_lattrait = gt_fit, embeddingVars = emb),
      occCovariates = occ, MCMCparams = MCMCparams)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  finished <- now()
  power_end <- power_state()

  jo <- fit$results_output$jsdm_output
  cp <- fit$infos$coef_prior
  if (!identical(cp$type, "grouped")) stop("the fit's prior is not grouped")
  if (!is.null(fit$infos$axis_move$sigma_b_acceptance))
    stop("this library has the sigma_b move; the test needs one without it")
  is_emb <- cp$is_emb
  if (!identical(colnames(fit$X_psi)[is_emb], emb))
    stop("the embedding columns of X_psi are not the simulated ones in order")
  # The measured block is 6 design columns: 4 numeric covariates and the 2
  # indicator columns of the categorical one. The true B has one row per design
  # column. With the fit's own design matrix, eta - B0 - X B must be the two
  # residual factors alone, so of rank 2; this checks rows, scaling and coding.
  if (!identical(is_emb, is_emb_true)) stop("the fit's columns are not the truth's rows in order")
  X_all <- as.matrix(fit$X_psi)
  left <- tp$eta - matrix(as.numeric(tp$B0), nrow(X_all), S, byrow = TRUE) - X_all %*% tp$B
  sv <- svd(left, nu = 0, nv = 0)$d
  if (sv[3] > 1e-8 * sv[1]) stop("the fit's design matrix does not match the true B")
  if (!identical(as.character(fit$infos$speciesNames), colnames(dat$OTU)))
    stop("the fit's species are not the simulated ones in order")
  if (!identical(dim(jo$B_output)[1:2], c(ncol(fit$X_psi), S)) ||
      dim(jo$B_output)[3] != niter || dim(jo$B_output)[4] != MCMCparams$nchain)
    stop("B_output is not design columns by species by sweeps by chains")
  X_emb <- as.matrix(fit$X_psi[, is_emb, drop = FALSE])

  # B_output is design columns by species by kept sweeps by chains.
  B_emb_mean <- apply(jo$B_output[is_emb, , , , drop = FALSE], c(1, 2), mean)
  rmse_contribution <- sqrt(mean((X_emb %*% (B_emb_mean - B_emb_true))^2))
  B_meas_draws <- jo$B_output[!is_emb, , , , drop = FALSE]
  B_meas_mean <- apply(B_meas_draws, c(1, 2), mean)
  B_meas_lo <- apply(B_meas_draws, c(1, 2), stats::quantile, probs = 0.025)
  B_meas_hi <- apply(B_meas_draws, c(1, 2), stats::quantile, probs = 0.975)
  B_meas_chain_means <- apply(B_meas_draws, c(1, 2, 4), mean)
  rmse_measured <- sqrt(mean((B_meas_mean - B_meas_true)^2))
  coverage_measured <- mean(B_meas_lo <= B_meas_true & B_meas_true <= B_meas_hi)
  width_measured <- mean(B_meas_hi - B_meas_lo)

  first <- seq_len(min(5000L, niter))
  rhat <- rbind(first_5000 = measured_rhat(jo$B_output, is_emb, first),
                all = measured_rhat(jo$B_output, is_emb, seq_len(niter)))
  conv <- tryCatch(returnGroupedConvergence(fit), error = function(e) conditionMessage(e))

  res <- list(
    label = label, dataset = k, n_lattrait = gt_fit,
    sigma_b = matrix(jo$sigmab_output, niter, MCMCparams$nchain),
    sigma_emb = matrix(jo$sigmaemb_output, niter, MCMCparams$nchain),
    rhat = rhat, convergence = conv,
    rmse_contribution = rmse_contribution, rmse_measured = rmse_measured,
    coverage_measured = coverage_measured, width_measured = width_measured,
    B_emb_true = B_emb_true, B_emb_mean = B_emb_mean,
    B_meas_true = B_meas_true, B_meas_mean = B_meas_mean,
    B_meas_lo = B_meas_lo, B_meas_hi = B_meas_hi, B_meas_chain_means = B_meas_chain_means,
    axis_move = fit$infos$axis_move, warnings = warns,
    minutes = minutes, started = started, finished = finished,
    power_start = power_start, power_end = power_end,
    lib = lib, commit = commit, R = as.character(getRversion()),
    n_measured = sum(!is_emb), n_emb = sum(is_emb),
    settings = list(datasettings = datasettings, jsdmParams = jsdmParams,
                    MCMCparams = MCMCparams, data_seed = 1000L + k, fit_seed = 2000L + k,
                    n_factors = 2L, n_lattrait = gt_fit))
  rm(fit, jo, B_meas_draws); gc()
  out_file <- file.path(out_dir, sprintf("%s-d%d.rds", label, k))
  tmp <- paste0(out_file, ".part")
  saveRDS(res, tmp)
  file.rename(tmp, out_file)
  sprintf("%s-d%d finished %s: %.1f minutes, RMSE embedding %.4f, measured %.4f, coverage %.3f, %d warnings; power %s",
          label, k, finished, minutes, rmse_contribution, rmse_measured, coverage_measured,
          length(warns), power_end)
}

cat("occJSDM", as.character(utils::packageVersion("occJSDM", lib.loc = lib)), "from", lib, "\n")
cat("arm", label, "; data sets", paste(datasets, collapse = ","), "; workers", workers,
    "; sweeps", nburn, "+", niter, "; commit", commit, "; output", out_dir, "\n")
cat("started", now(), "\n")
cat("power at start:", power_state(), "\n")
out <- parallel::mclapply(datasets, function(k) {
  tryCatch(fit_one(k), error = function(e) paste0(label, "-d", k, " FAILED: ", conditionMessage(e)))
}, mc.cores = workers, mc.preschedule = FALSE)
cat(unlist(out), sep = "\n")
cat("power at end:", power_state(), "\n")
cat("ended", now(), "\n")
