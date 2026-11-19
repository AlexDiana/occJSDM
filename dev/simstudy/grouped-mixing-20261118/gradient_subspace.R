#!/usr/bin/env Rscript
# Are the habitat gradients themselves recovered, and do the chains agree on them? Scores the fits
# of rotation_test.R and orientation_test.R on the gradients as a set, not one by one: the
# three-dimensional subspace spanned by the leading directions of a chain's posterior mean of the
# embedding's contribution (sites by species), in site space and in species space. Three, because
# the simulated truth has three blends. Reads the per-run _summary.rds files (kept locally, not
# tracked), so it runs only where the fits were made.
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/gradient_subspace.R [--write]
d <- "dev/simstudy/grouped-mixing-20261118/results"
files <- c(list.files(file.path(d, "rotation"), pattern = "_summary.rds$", full.names = TRUE),
           list.files(file.path(d, "orientation"), pattern = "_summary.rds$", full.names = TRUE))
r <- 3L
top <- function(M) svd(scale(M, scale = FALSE), nu = r, nv = r)
agree <- function(U, V) sum((t(U) %*% V)^2) / r        # 1 = the same r-dimensional subspace
out <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f); truth <- x$true_contribution; st <- top(truth); nch <- dim(x$contribution_mean)[3]
  est <- lapply(seq_len(nch), function(k) top(x$contribution_mean[, , k]))
  pairs <- utils::combn(nch, 2)
  site_chain <- apply(pairs, 2, function(p) agree(est[[p[1]]]$u, est[[p[2]]]$u))
  species_chain <- apply(pairs, 2, function(p) agree(est[[p[1]]]$v, est[[p[2]]]$v))
  tc <- scale(truth, scale = FALSE)
  # the share of the true contribution that lies in a chain's three estimated site gradients
  captured <- sapply(est, function(e) sum((t(e$u) %*% tc)^2) / sum(tc^2))
  data.frame(fit = gsub("_b5000_summary.rds", "", basename(f), fixed = TRUE),
             truth_in_its_top3 = sum(st$d[1:r]^2) / sum(st$d^2),
             truth_direction_shares = paste(round(st$d[1:r]^2 / sum(st$d^2), 2), collapse = " "),
             captured_min = min(captured), captured_max = max(captured),
             site_between_chains_min = min(site_chain), site_between_chains_max = max(site_chain),
             species_between_chains_min = min(species_chain), species_between_chains_max = max(species_chain))
}))
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], round, 3)
print(out, row.names = FALSE)
if ("--write" %in% commandArgs(TRUE)) utils::write.csv(out, file.path(d, "gradient_subspace.csv"), row.names = FALSE)
