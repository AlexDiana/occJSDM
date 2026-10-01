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

# The stored lattice verification shares the package's basis function, so this check builds
# the whitened lattice basis from the kernel alone: standardise each axis with the site mean
# and sd of the fit's raw coordinates, take the squared-exponential kernel to the knots, and
# multiply by the inverse transpose of chol(knot kernel + 1e-5 I), as independent_bases does for sites.
independent_lattice_basis_difference <- function(fit, site_xy, lattice_xy, n_cells = 50L) {
  set.seed(1L); cells <- sort(sample(nrow(lattice_xy), n_cells))
  raw <- lattice_xy[cells, , drop = FALSE]
  scaled <- sweep(sweep(raw, 2, colMeans(site_xy), "-"), 2, apply(site_xy, 2, sd), "/")
  knots <- fit$infos$list_Xs$X_tilde
  kernel <- function(a, b, l) exp(-(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2) / (2 * l^2))
  mine <- lapply(fit$infos$l_s_grid, function(l) {
    lower <- t(chol(kernel(knots, knots, l) + diag(1e-5, nrow(knots))))
    t(solve(lower, t(kernel(scaled, knots, l))))
  })
  native <- occJSDM:::createSpatialPredMatrix(data.frame(longitude = raw[, 1], latitude = raw[, 2]),
    fit$infos$l_s_grid, knots, fit$infos$list_Xs_mat)
  max(vapply(seq_along(mine), function(g) max(abs(mine[[g]] - native[, , g])), numeric(1)))
}

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
  basis_diff <- independent_lattice_basis_difference(fit, input$surveys[[r$job$arrangement]]$truth$xy, land$points[li, ])
  rows[[i]] <- data.frame(key = sel$key[i], phase = sel$phase[i], probability = prob_diff, field_median = field_diff,
    range = range_diff, lattice = lattice_diff, independent_basis = basis_diff,
    native_lattice = r$verification$lattice, basis = r$verification$basis,
    passed = prob_diff < 1e-10 & field_diff < 1e-10 & range_diff < 1e-12 & lattice_diff < 1e-10 &
      basis_diff < 1e-8 & r$verification$lattice < 1e-8)
  cat(sel$key[i], if (rows[[i]]$passed) "ok" else "MISMATCH", "\n")
}
oracle <- read.csv(file.path(study, "oracle/oracle-selected.csv"))
oracle_ok <- all(unname(tools::md5sum(oracle$file)) == oracle$result_md5)
audit <- do.call(rbind, rows); dir.create(file.path(study, "audit"), showWarnings = FALSE)
write.csv(audit, file.path(study, "audit/audit.csv"), row.names = FALSE)
writeLines(c(paste("oracle result hashes match:", oracle_ok), paste("fits passed:", sum(audit$passed), "of", nrow(audit))),
           file.path(study, "audit/summary.txt"))
stopifnot(oracle_ok, all(audit$passed)); cat("Audit passed for", nrow(audit), "fits and 96 oracle results.\n")
