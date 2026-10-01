#!/usr/bin/env Rscript
# Convergence of the site field values in every selected fit, from the saved draws.
# The protocol's convergence rule covers occupancy probabilities and model parameters but not the
# site field, so this summary was computed after the fits (PLAN.md amendment 3, clause f). For each
# selected fit it rebuilds the 100 sites x 8 species field draws with score.R's reconstruct_field_draws()
# at thin = 1, reshapes each element to iterations x chains, and takes posterior::rhat and
# posterior::ess_bulk (and ess_mean) per element. Changes no fit and no score. Run after verify.R.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) { hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }; substring(hit, nchar(name) + 4L) }
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
.libPaths(c(file.path(study, "library"), .libPaths())); suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == normalizePath(file.path(study, "library/occJSDM")))
source(file.path(repo, "dev/simstudy/spatial-design-sweep/score.R"))

text <- c("initial_reasons", "selected_reasons", "warnings")
sel <- read.csv(file.path(study, "summary-final/selected-fits.csv"),
                colClasses = setNames(rep("character", length(text)), text))
rows <- list()
for (i in seq_len(nrow(sel))) {
  stopifnot(unname(tools::md5sum(sel$selected_file[i])) == sel$selected_md5[i])
  fit <- readRDS(sub("-result.rds$", "-fit.rds", sel$selected_file[i]))$fit
  js <- fit$results_output$jsdm_output; ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  d <- reconstruct_field_draws(fit, thin = 1L)  # sites x species x (iterations * chains), chain-major
  n <- dim(d)[1]; S <- dim(d)[2]; stopifnot(dim(d)[3] == ni * nc)
  rh <- bulk <- em <- matrix(NA_real_, n, S)
  for (s in seq_len(S)) for (j in seq_len(n)) {
    m <- matrix(d[j, s, ], ni, nc)
    rh[j, s] <- posterior::rhat(m); bulk[j, s] <- posterior::ess_bulk(m); em[j, s] <- posterior::ess_mean(m)
  }
  rows[[i]] <- data.frame(key = sel$key[i], community = sel$community[i], arrangement = sel$arrangement[i],
    arm = sel$arm[i], phase = sel$phase[i], n_elements = length(rh), max_rhat = max(rh),
    min_ess_bulk = min(bulk), min_ess_mean = min(em), n_rhat_over_1.05 = sum(rh > 1.05))
  print(rows[[i]]); rm(fit, d); invisible(gc())
}
out <- do.call(rbind, rows); dir.create(file.path(study, "audit"), showWarnings = FALSE)
write.csv(out, file.path(study, "audit/field-convergence.csv"), row.names = FALSE)
cat("max field-element Rhat", max(out$max_rhat), "min bulk ESS", min(out$min_ess_bulk), "\n")
