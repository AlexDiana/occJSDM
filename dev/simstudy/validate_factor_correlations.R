# Paired post-processing diagnostic; no nominal-coverage claim.
# From the package root:
# Rscript dev/simstudy/validate_factor_correlations.R [replicates=3] [outdir] [before_ref=90026d6]
# Fits once per dataset, capturing the actual raw factors at the final
# reparameterisation boundary, then applies old and new transforms to the
# same draws. The sampler is identical by construction.
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args)) as.integer(args[1]) else 3L
out_dir <- if (length(args) >= 2) args[2] else tempfile("factor-correlations-")
before_ref <- if (length(args) >= 3) args[3] else "90026d6"
stopifnot(length(n_reps) == 1, !is.na(n_reps), n_reps > 0)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
devtools::load_all(quiet = TRUE)
RcppParallel::setThreadOptions(numThreads = 1)

before_source <- system2("git", c("show", shQuote(paste0(before_ref, ":R/jsdmfun.R"))),
                         stdout = TRUE)
stopifnot(is.null(attr(before_source, "status")))
before_env <- new.env(parent = asNamespace("occJSDM"))
eval(parse(text = before_source), envir = before_env)
after_source <- readLines("R/jsdmfun.R", warn = FALSE)
fit_env <- new.env(parent = asNamespace("occJSDM"))
sys.source("R/runOccJSDM.R", envir = fit_env)
raw_factors <- list()
fit_env$reparamFactorModel <- function(U_output, L_output) {
  raw_factors[[length(raw_factors) + 1L]] <<- list(U_output = U_output, L_output = L_output)
  occJSDM:::reparamFactorModel(U_output, L_output)
}

# Deliberately dense, non-boundary correlations, unlike the sparse simulator
# default. All draws otherwise follow its U ~ N(0, sigma_h^2 I), continuous
# Normal(eta, tau^2) and binary Bernoulli(plogis(eta)) generative equations.
L_true <- rbind(c(2, 1, 0.5, -1.5, -0.5, 1.5),
                c(0, 1, -1.5, 0.5, -1, 0.75))
S <- ncol(L_true)
sigma_h <- 1
B0 <- c(-0.4, 0.3, -0.2, 0.5, 0.2, -0.3)
B <- c(0.2, -0.4, 0.3, -0.2, 0.4, -0.3)
truth <- cov2cor(crossprod(L_true))
pairs <- which(lower.tri(truth), arr.ind = TRUE)
stopifnot(all(is.finite(truth)), max(abs(truth[lower.tri(truth)])) < 0.99)
mcmc <- list(nchain = 2, nburn = 1000, niter = 1000, nthin = 1)
rows <- checks <- list()
started <- Sys.time()

for (model in c("continuous", "binary")) {
  for (replicate in seq_len(n_reps)) {
    seed <- 61900L + 100L * match(model, c("continuous", "binary")) + replicate
    set.seed(seed)
    n <- if (model == "continuous") 500L else 1000L
    X <- as.numeric(scale(rnorm(n)))
    U <- matrix(rnorm(n * 2, sd = sigma_h), n, 2)
    eta <- matrix(B0, n, S, byrow = TRUE) + X %o% B + U %*% L_true
    y <- if (model == "continuous") eta + matrix(rnorm(n * S, sd = 0.5), n, S) else
      matrix(rbinom(n * S, 1, plogis(eta)), n, S)
    colnames(y) <- paste0("species", seq_len(S))
    dat <- list(info = data.frame(Site = seq_len(n), X = X), OTU = y, traits = NULL)
    raw_factors <- list()
    fit <- suppressMessages(fit_env$runOccJSDM(
      dat, listParams = list(n_factors = 2), occCovariates = "X", MCMCparams = mcmc))
    stopifnot(length(raw_factors) >= 1L)
    raw <- raw_factors[[1]]
    old <- do.call(before_env$reparamFactorModel, raw)
    new <- fit$results_output$jsdm_output[c("U_output", "L_output")]
    versions <- list(raw = raw, before = old, after = new)
    for (version in names(versions)) {
      current <- versions[[version]]
      current_fit <- fit
      current_fit$results_output$jsdm_output$L_output <- current$L_output
      ci <- returnResidualCorrelationMatrix(current_fit)
      key <- paste(model, replicate, version, sep = ":")
      rows[[key]] <- data.frame(
        model, replicate, version, species1 = pairs[, 1], species2 = pairs[, 2],
        truth = truth[pairs], estimate = ci[2, , ][pairs],
        lower = ci[1, , ][pairs], upper = ci[3, , ][pairs])
      predictor_error <- covariance_error <- correlation_error <- 0
      for (chain in seq_len(mcmc$nchain)) for (iter in seq_len(mcmc$niter)) {
        before_L <- raw$L_output[, , iter, chain]
        after_L <- current$L_output[, , iter, chain]
        predictor_error <- max(predictor_error, abs(
          current$U_output[, , iter, chain] %*% after_L -
            raw$U_output[, , iter, chain] %*% before_L))
        covariance_error <- max(covariance_error, abs(crossprod(after_L) - crossprod(before_L)))
        correlation_error <- max(correlation_error, abs(
          cov2cor(crossprod(after_L)) - cov2cor(crossprod(before_L))))
      }
      checks[[key]] <- data.frame(model, replicate, version,
                                  predictor_error, covariance_error, correlation_error)
    }
    saveRDS(list(data = dat, truth = list(L = L_true, sigma_h = sigma_h, eta = eta,
                                         B0 = B0, B = B), fit = fit, raw = raw, seed = seed),
             file.path(out_dir, sprintf("%s-%d.rds", model, replicate)))
    cat(sprintf("%s replicate %d/%d complete\n", model, replicate, n_reps))
    saveRDS(list(rows = do.call(rbind, rows), checks = do.call(rbind, checks)),
             file.path(out_dir, "paired-results.rds"))
  }
}

result <- do.call(rbind, rows)
result$bias <- result$estimate - result$truth
result$abs_error <- abs(result$bias)
result$squared_error <- result$bias^2
result$width <- result$upper - result$lower
result$sign_correct <- sign(result$estimate) == sign(result$truth)
summary <- aggregate(cbind(bias, abs_error, squared_error, width, sign_correct) ~ model + version,
                       result, mean)
summary$rmse <- sqrt(summary$squared_error)
summary$squared_error <- NULL
check_rows <- do.call(rbind, checks)
stopifnot(max(subset(check_rows, version == "after")$predictor_error) < 1e-10,
          max(subset(check_rows, version == "after")$covariance_error) < 1e-10,
          max(subset(check_rows, version == "after")$correlation_error) < 1e-10)
write.csv(result, file.path(out_dir, "correlation-pairs.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(check_rows, file.path(out_dir, "invariants.csv"), row.names = FALSE)
saveRDS(list(results = result, summary = summary, checks = check_rows,
              mcmc = mcmc, n_reps = n_reps, before_ref = before_ref,
              before_source = before_source, after_source = after_source,
              fitting_source = readLines("R/runOccJSDM.R", warn = FALSE),
              started = started, finished = Sys.time(), sessionInfo = sessionInfo()),
         file.path(out_dir, "validation.rds"))
print(summary)
