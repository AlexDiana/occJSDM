# The same simulated data with and without the cap (evidence step 3 of
# dev/superpowers/specs/2026-11-25-axis-slab-design.md).
#
# The five data sets of the sigma_b move's simulation
# (dev/simstudy/sigma-b-move-20261122/repro.R): simulateOccJSDMData() with 256
# sites, 48 species, 5 measured covariates, 64 embedding columns, 2 traits, 4
# latent traits, model "binary", data set k simulated with seed 1000 + k and
# fitted with seed 2000 + k (set.seed and setOccJSDMSeed). Each fit has 4
# chains of 5,000 burn-in and 5,000 kept sweeps, run one after another in its
# worker (ncores = 1), with n_factors = 2, n_lattrait = 4 and the embedding
# columns under the grouped prior. --slab 2 fits with listPriors =
# list(slab_scale = 2), --slab Inf with listPriors = list(slab_scale = Inf);
# any other number is passed the same way. The fit's prior record is checked
# to hold the value asked for. No dataset content is read.
#
# Usage, from anywhere:
#   Rscript --vanilla compare.R --lib PATH --slab 2|Inf --datasets 1:5
#     --workers 3 [--commit SHA] [--outdir DIR] [--nburn 5000] [--niter 5000]
# --workers is the number of data sets fitted at once (forked workers, one
# fit each). --nburn, --niter and --outdir exist for the short check of the
# script; --outdir defaults to dev/simstudy/results/axis-slab-20261125/.
#
# Writes slab<VALUE>-d<k>.rds per data set, where VALUE is the --slab value
# as given: a list with
#   - B_emb_true: the true embedding rows of B, embedding columns by species
#     (simulateOccJSDMData()'s jsdmParams_true$B, rows labelled "embedding");
#   - B_emb_mean: their posterior means over the kept sweeps of all chains,
#     in the same orientation and with the same dimnames (checked), and
#     B_emb_chain_means, the same per chain;
#   - rmse_contribution: the root mean squared error, over sites and species,
#     of the posterior mean of the embedding's contribution X_emb %*% B_emb
#     against the truth's X_emb %*% B_emb_true. The posterior mean is computed
#     as X_emb %*% B_emb_mean, which equals the mean over draws of
#     X_emb %*% B_emb by linearity. X_emb is the fit's own standardised
#     embedding block, checked to equal the simulated columns;
#   - convergence: returnGroupedConvergence() in full, with checks and
#     verdict also stored on their own;
#   - axis_scale_median: the pooled median, over kept sweeps and chains, of
#     the stored axis scale of each embedding column (lambda_output, which
#     under a cap is the effective axis scale);
#   - axis_move: the fit's infos$axis_move;
#   - minutes, the power state and the full `pmset -g` output at the start and
#     end of the fit, the commit given, the library path, and the settings.

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
slab_arg <- arg("slab")
slab <- suppressWarnings(as.numeric(slab_arg))
if (length(slab) != 1L || is.na(slab) || slab <= 0)
  stop("--slab must be a positive number or Inf")
label <- paste0("slab", slab_arg)
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
if (is.na(workers) || workers < 1L || workers > 3L)
  stop("--workers must be 1, 2 or 3 for these runs")
nburn <- as.integer(arg("nburn", "5000"))
niter <- as.integer(arg("niter", "5000"))
commit <- arg("commit", NA_character_)

library(occJSDM, lib.loc = lib)
pkg_path <- normalizePath(dirname(system.file(package = "occJSDM")))
if (pkg_path != lib) stop("occJSDM was loaded from ", pkg_path, ", not from ", lib)
RcppParallel::setThreadOptions(numThreads = 1)

# The study directory, found from this script's own path.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
study_dir <- dirname(normalizePath(script))
out_dir <- arg("outdir", file.path(study_dir, "..", "results", "axis-slab-20261125"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
script_head <- tryCatch(
  system2("git", c("-C", shQuote(study_dir), "rev-parse", "--short", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_)

power_state <- function() {
  paste(c(system("pmset -g batt", intern = TRUE)[1],
          grep("lowpowermode", system("pmset -g", intern = TRUE), value = TRUE)),
        collapse = "; ")
}
pmset_full <- function() system("pmset -g", intern = TRUE)
now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "GMT", usetz = TRUE)

S <- 48L; n <- 256L
datasettings <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L,
                     K = rep(1L, n), ncov_psi = 5L, ncov_theta = 1L)
jsdmParams <- list(gt = 4L, d = 2L, ds = 0, sigma_b = 0.1, sigma_bs = 0.5, sigma_ts = 0.5,
                   sigma_h = 1, sigma_s = 0.5, l_s = 0.3, tau = rep(1, S),
                   useSpatField = FALSE, embedding = list(n_emb = 64L))
MCMCparams <- list(nchain = 4L, nburn = nburn, niter = niter, nthin = 1L, ncores = 1L)
listPriors <- list(slab_scale = slab)

fit_one <- function(k) {
  set.seed(1000L + k)
  sim <- simulateOccJSDMData(datasettings, list(), jsdmParams, model = "binary")
  dat <- sim$data_list
  emb <- dat$embeddingVars
  occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), emb), emb)
  tp <- sim$true_params$jsdmParams_true
  is_emb_true <- tp$covariate_group == "embedding"
  if (sum(is_emb_true) != length(emb) || nrow(tp$B) != length(tp$covariate_group) ||
      ncol(tp$B) != S)
    stop("the true B does not have one embedding row per embedding column")
  B_emb_true <- tp$B[is_emb_true, , drop = FALSE]
  dimnames(B_emb_true) <- list(emb, colnames(dat$OTU))

  power_start <- power_state()
  pmset_start <- pmset_full()
  started <- now()
  cat(sprintf("%s-d%d started %s; power %s\n", label, k, started, power_start))
  set.seed(2000L + k)
  occJSDM:::setOccJSDMSeed(2000L + k)
  warns <- character()
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    suppressMessages(runOccJSDM(dat,
      listParams = list(n_factors = 2L, n_lattrait = 4L, embeddingVars = emb),
      listPriors = listPriors,
      occCovariates = occ, MCMCparams = MCMCparams)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  finished <- now()
  power_end <- power_state()
  pmset_end <- pmset_full()

  jo <- fit$results_output$jsdm_output
  cp <- fit$infos$coef_prior
  if (!identical(cp$type, "grouped")) stop("the fit's prior is not grouped")
  if (!identical(unname(cp$slab_scale), slab))
    stop("the fit's slab_scale is ", format(cp$slab_scale), ", not ", format(slab))
  if (!identical(fit$infos$model, "binary")) stop("the fit's model type is not binary")
  is_emb <- cp$is_emb
  if (!identical(colnames(fit$X_psi)[is_emb], emb))
    stop("the embedding columns of X_psi are not the simulated ones in order")
  if (!identical(as.character(fit$infos$speciesNames), colnames(dat$OTU)))
    stop("the fit's species are not the simulated ones in order")
  X_emb <- as.matrix(fit$X_psi[, is_emb, drop = FALSE])
  X_sim <- as.matrix(dat$info[, emb, drop = FALSE])
  dimnames(X_emb) <- dimnames(X_sim) <- NULL
  if (!isTRUE(all.equal(X_emb, X_sim)))
    stop("the fit's standardised X_emb differs from the simulated columns")
  if (!identical(dim(jo$B_output)[1:2], c(ncol(fit$X_psi), S)) ||
      dim(jo$B_output)[3] != niter || dim(jo$B_output)[4] != MCMCparams$nchain)
    stop("B_output is not design columns by species by sweeps by chains")

  # B_output is design columns by species by kept sweeps by chains.
  B_emb_draws <- jo$B_output[is_emb, , , , drop = FALSE]
  B_emb_chain_means <- apply(B_emb_draws, c(1, 2, 4), mean)
  B_emb_mean <- apply(B_emb_draws, c(1, 2), mean)
  dimnames(B_emb_mean) <- dimnames(B_emb_true)
  dimnames(B_emb_chain_means) <- c(dimnames(B_emb_true), list(NULL))
  if (!identical(dim(B_emb_mean), dim(B_emb_true))) stop("B_emb_mean and B_emb_true differ in shape")
  contrib_err <- X_emb %*% (B_emb_mean - B_emb_true)
  rmse_contribution <- sqrt(mean(contrib_err^2))

  conv <- returnGroupedConvergence(fit)
  # lambda_output is embedding columns by kept sweeps by chains.
  axis_scale_median <- apply(jo$lambda_output, 1, stats::median)
  names(axis_scale_median) <- emb

  res <- list(
    label = label, dataset = k, slab_scale = cp$slab_scale, listPriors = listPriors,
    B_emb_true = B_emb_true, B_emb_mean = B_emb_mean,
    B_emb_chain_means = B_emb_chain_means,
    X_emb = X_emb, rmse_contribution = rmse_contribution,
    contribution_mean = "X_emb %*% B_emb_mean (equal by linearity to the mean over draws of X_emb %*% B_emb)",
    convergence = conv, checks = conv$checks, verdict = conv$verdict,
    axis_scale_median = axis_scale_median,
    sigma_emb = jo$sigmaemb_output, sigma_b = jo$sigmab_output,
    axis_move = fit$infos$axis_move,
    warnings = warns,
    minutes = minutes, started = started, finished = finished,
    power_start = power_start, power_end = power_end,
    pmset_start = pmset_start, pmset_end = pmset_end,
    lib = lib, commit = commit, script_head = script_head,
    R = as.character(getRversion()),
    settings = list(datasettings = datasettings, jsdmParams = jsdmParams,
                    MCMCparams = MCMCparams, data_seed = 1000L + k,
                    fit_seed = 2000L + k, n_factors = 2L, n_lattrait = 4L))
  rm(fit, jo, B_emb_draws); gc()
  out_file <- file.path(out_dir, sprintf("%s-d%d.rds", label, k))
  tmp <- paste0(out_file, ".part")
  saveRDS(res, tmp)
  file.rename(tmp, out_file)
  sprintf("%s-d%d finished %s: %.1f minutes, RMSE %.4f, verdict %s, %d warnings; power %s",
          label, k, finished, minutes, rmse_contribution, conv$verdict, length(warns), power_end)
}

cat("occJSDM", as.character(utils::packageVersion("occJSDM", lib.loc = lib)), "from", lib, "\n")
cat("arm", label, ": listPriors =", deparse(listPriors), "; data sets",
    paste(datasets, collapse = ","), "; workers", workers, "; sweeps", nburn, "+", niter, "\n")
cat("commit", commit, "; script at HEAD", script_head, "; output", out_dir, "\n")
cat("started", now(), "\n")
cat("power at start:", power_state(), "\n")
out <- parallel::mclapply(datasets, function(k) {
  tryCatch(fit_one(k), error = function(e) paste0(label, "-d", k, " FAILED: ", conditionMessage(e)))
}, mc.cores = workers, mc.preschedule = FALSE)
cat(unlist(out), sep = "\n")
cat("power at end:", power_state(), "\n")
cat("ended", now(), "\n")
