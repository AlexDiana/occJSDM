# Does storing every nthin-th draw, with z and psi kept, change the chain?
# Fits the same dataset twice in one process: once as the main run fitted it
# (nthin = 1, summarisedLatentPresences = TRUE) and once as refit.R fits it
# (nthin = 10, a tenth of the kept draws, summarisedLatentPresences = FALSE),
# both for the same number of iterations, and compares the draws at every
# 10th iteration. Also checks that the thinned fit's stored psi equals
# logistic(eta) rebuilt from its stored terms (terms.R). Two settings: the
# 10-species rare arm, replicate 1, at 50 + 100 iterations, and the
# 100-species rare arm, replicate 3, at 20 + 40. Under L'Ecuyer-CMRG with one
# thread, as refit.R's workers run.
# Usage: Rscript check-storage.R <package root> <output csv>
args <- commandArgs(trailingOnly = TRUE)
pkg_root <- normalizePath(args[1], mustWork = TRUE)
out_csv <- args[2]
run_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
RNGkind("L'Ecuyer-CMRG")
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
RcppParallel::setThreadOptions(numThreads = 1)
suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))
source(file.path(run_dir, "terms.R"))

fit_once <- function(arm, replicate, nburn, niter, nthin, summarised) {
  sc <- simstudy_scenarios()
  s <- sc[[which(vapply(sc, `[[`, "", "label") == arm)]]
  truth <- draw_truth(s, simstudy_seed_for(s, replicate))
  sim <- simstudy_simulate(truth)
  a <- simstudy_fit_args(sim, s)
  fit <- suppressMessages(suppressWarnings(
    runOccJSDM(a$data, listParams = a$listParams, listPriors = a$listPriors,
               occCovariates = a$occCovariates, collCovariates = a$collCovariates,
               spatCovariates = a$spatCovariates,
               MCMCparams = list(nchain = 2, nburn = nburn, niter = niter, nthin = nthin),
               summarisedLatentPresences = summarised)))
  list(fit = fit, sim = sim)
}
# Every keep-th index along the draw dimension, the second to last.
thin_draws <- function(a, keep) {
  d <- length(dim(a))
  idx <- seq(keep, dim(a)[d - 1], by = keep)
  if (d == 2) a[idx, , drop = FALSE] else if (d == 3) a[, idx, , drop = FALSE] else
    a[, , idx, , drop = FALSE]
}
compare <- function(label, full, thinned, keep) {
  jf <- full$results_output$jsdm_output; jt <- thinned$results_output$jsdm_output
  q <- list(mu_b0 = list(jf$mu_b0_output, jt$mu_b0_output),
            B0 = list(jf$B0_output, jt$B0_output),
            B = list(jf$B_output, jt$B_output),
            Bs = list(jf$Bs_output, jt$Bs_output),
            U = list(jf$U_output, jt$U_output),
            L = list(jf$L_output, jt$L_output),
            idx_ls = list(jf$idx_ls_output, jt$idx_ls_output),
            sigmah = list(jf$sigmah_output, jt$sigmah_output),
            beta_theta = list(full$results_output$beta_theta_output,
                              thinned$results_output$beta_theta_output),
            theta0 = list(full$results_output$theta0_output, thinned$results_output$theta0_output),
            p = list(full$results_output$p_output, thinned$results_output$p_output),
            chain_seeds = list(full$infos$chain_seeds, thinned$infos$chain_seeds))
  do.call(rbind, lapply(names(q), function(nm) {
    a <- q[[nm]][[1]]; b <- q[[nm]][[2]]
    if (nm != "chain_seeds") a <- thin_draws(a, keep)
    data.frame(setting = label, quantity = nm, n = length(b),
               max_abs_difference = max(abs(a - b)),
               identical = identical(dim(a), dim(b)) && all(a == b))
  }))
}
# The stored psi of the thinned fit against logistic(eta) rebuilt from its
# stored terms, over every stored draw of both chains.
psi_check <- function(label, kept) {
  d <- kept$draws
  diffs <- vapply(seq_len(dim(d$B0)[3]), function(ch) {
    max(vapply(seq_len(dim(d$B0)[2]), function(t) {
      max(abs(plogis(eta_terms(d, kept$design, t, ch)$eta) - d$psi[, , t, ch]))
    }, 0))
  }, 0)
  data.frame(setting = label, quantity = "psi_vs_rebuilt_eta", n = length(d$psi),
             max_abs_difference = max(diffs), identical = max(diffs) < 1e-10)
}
settings <- list(
  list(label = "rare_m3_rep1_50+100", arm = "base_rare_cm_m3", replicate = 1L, nburn = 50L, niter = 100L),
  list(label = "rare_s100_rep3_20+40", arm = "base_rare_cm_s100", replicate = 3L, nburn = 20L, niter = 40L))
res <- do.call(rbind, lapply(settings, function(st) {
  message(sprintf("%s: %s", format(Sys.time(), "%H:%M:%S"), st$label))
  full <- fit_once(st$arm, st$replicate, st$nburn, st$niter, 1L, TRUE)
  thinned <- fit_once(st$arm, st$replicate, st$nburn, st$niter %/% 10L, 10L, FALSE)
  kept <- keep_fit(thinned$fit, thinned$sim)
  rbind(compare(st$label, full$fit, thinned$fit, 10L), psi_check(st$label, kept))
}))
write.csv(res, out_csv, row.names = FALSE)
print(res, row.names = FALSE)
message(sprintf("%s: all identical: %s", format(Sys.time(), "%H:%M:%S"), all(res$identical)))
