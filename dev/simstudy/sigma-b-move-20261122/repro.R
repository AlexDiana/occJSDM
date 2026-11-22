# Looks for chains whose sigma_b (the scale of the measured covariates'
# species leftovers under the grouped prior) falls near zero and stays there.
# Evidence steps 1 and 2 of dev/superpowers/specs/2026-11-22-sigma-b-move-design.md:
# the same simulated data sets and seeds are fitted once with a library built
# from a2a64a1 (label "base") and once with a library holding the sigma_b move
# (label "move"). Every data set comes from simulateOccJSDMData(); no dataset
# content is read.
#
# Usage:
#   Rscript --vanilla repro.R --lib PATH --label LABEL --datasets 1:5 --workers 3
#     [--commit SHA] [--nburn 5000] [--niter 20000]
# --nburn and --niter exist for the short check of the script only.
#
# Writes dev/simstudy/results/sigma-b-move-20261122/LABEL-d<k>.rds per data
# set: a list with sigma_b (kept sweeps by chains), rhat, axis_move (the fit's
# infos$axis_move), minutes, power, lib, commit and settings.

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
label <- arg("label")
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
out_dir <- file.path(study_dir, "..", "results", "sigma-b-move-20261122")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)

# The commit the library was built from, as given on the command line.
commit <- arg("commit", NA_character_)

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
# chains, and the number of measured coefficients over 1.1.
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
  occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars),
           dat$embeddingVars)
  power <- power_state()
  set.seed(2000L + k)
  occJSDM:::setOccJSDMSeed(2000L + k)
  t0 <- Sys.time()
  fit <- suppressMessages(runOccJSDM(dat,
    listParams = list(n_factors = 2L, n_lattrait = 4L, embeddingVars = dat$embeddingVars),
    occCovariates = occ, MCMCparams = MCMCparams))
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  jo <- fit$results_output$jsdm_output
  is_emb <- fit$infos$coef_prior$is_emb
  first <- seq_len(min(5000L, niter))
  rhat <- rbind(first_5000 = measured_rhat(jo$B_output, is_emb, first),
                all = measured_rhat(jo$B_output, is_emb, seq_len(niter)))
  res <- list(sigma_b = matrix(jo$sigmab_output, niter, MCMCparams$nchain),
              rhat = rhat, axis_move = fit$infos$axis_move, minutes = minutes,
              power = power, lib = lib, commit = commit, dataset = k,
              n_measured = sum(!is_emb), n_emb = sum(is_emb),
              settings = list(datasettings = datasettings, jsdmParams = jsdmParams,
                              MCMCparams = MCMCparams))
  rm(fit, jo); gc()
  saveRDS(res, file.path(out_dir, sprintf("%s-d%d.rds", label, k)))
  sprintf("%s-d%d: %.1f minutes", label, k, minutes)
}

cat("occJSDM", as.character(utils::packageVersion("occJSDM", lib.loc = lib)), "from", lib, "\n")
cat("power at start:", power_state(), "\n")
out <- parallel::mclapply(datasets, function(k) {
  tryCatch(fit_one(k), error = function(e) paste0(label, "-d", k, " FAILED: ", conditionMessage(e)))
}, mc.cores = workers, mc.preschedule = FALSE)
cat(unlist(out), sep = "\n")
