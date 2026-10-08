#!/usr/bin/env Rscript
# Planning feasibility audit only. Does not save inputs, fit models or alter the old study.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[1] else ".")
bundle_path <- file.path(repo, "vignettes/teaching-data/spatial-lesson.rds")
generator_path <- file.path(repo, "dev/simstudy/spatial-design-sweep/generator.R")
b <- readRDS(bundle_path)
stopifnot(unname(tools::md5sum(generator_path)) == unname(b$provenance$script_hashes["generator.R"]))
source(generator_path)
tolerance <- 1e-12
checks <- list()
replay <- list()
record <- function(community, check, n, difference = 0, exact = FALSE) {
  stopifnot(is.finite(difference), difference <= if (exact) 0 else tolerance)
  checks[[length(checks) + 1L]] <<- data.frame(community, check, n,
    max_absolute_difference = difference, criterion = if (exact) "exact" else "absolute <= 1e-12")
}
for (r in 1:3) {
  inp <- make_sweep_input(r)
  id <- sprintf("rep%02d", r)
  if (r == 1) for (nm in names(b$landscape)) {
    old <- b$landscape[[nm]]; new <- inp$landscape[[nm]]
    if (is.numeric(old)) {
      stopifnot(identical(dim(old), dim(new)))
      record(id, paste0("landscape_", nm), length(old), max(abs(old - new)))
    } else {
      stopifnot(identical(old, new))
      record(id, paste0("landscape_", nm), length(old), exact = TRUE)
    }
  }
  ar <- b$arrangements[b$arrangements$community == id, ]
  stopifnot(nrow(ar) == 400L)
  xy <- t(mapply(function(a, i) inp$surveys[[a]]$truth$xy[i, ], ar$arrangement, ar$site))
  stopifnot(identical(unname(as.matrix(ar[c("x", "y")])), unname(xy)))
  record(id, "survey_coordinates", length(xy), exact = TRUE)
  fm <- b$field_maps[b$field_maps$community == id & b$field_maps$source == "truth", ]
  stopifnot(nrow(fm) == 800L)
  values <- mapply(function(a, i, s) inp$surveys[[a]]$truth$field[i, s], fm$arrangement, fm$site, fm$species)
  record(id, "saved_two_species_fields", length(values), max(abs(values - fm$value)))
  oracle <- b$oracle[b$oracle$community == id, ]
  stopifnot(nrow(oracle) == 32L)
  null_rmse <- mapply(function(a, s) {
    f <- inp$surveys[[a]]$truth$field[, s]; sqrt(mean((f - mean(f))^2))
  }, oracle$arrangement, oracle$species)
  record(id, "all_species_zero_field_rmse", length(null_rmse), max(abs(null_rmse - oracle$zero_field_rmse)))
  stats <- b$statistics[b$statistics$replicate == r, ]
  stats <- stats[match(inp$statistics$arrangement, stats$arrangement), ]
  stopifnot(nrow(stats) == 4L)
  for (nm in names(stats)) {
    if (is.numeric(stats[[nm]])) {
      record(id, paste0("design_", nm), nrow(stats), max(abs(inp$statistics[[nm]] - stats[[nm]])))
    } else stopifnot(identical(inp$statistics[[nm]], stats[[nm]]))
  }
  ss <- b$fits$species[b$fits$species$community == id, ]
  stopifnot(nrow(ss) == 64L)
  occupied <- mapply(function(a, s) sum(inp$surveys[[a]]$truth$z[, s]), ss$arrangement, ss$species)
  detections <- mapply(function(a, arm, s) sum(inp$surveys[[a]][[arm]]$OTU[, s]), ss$arrangement, ss$arm, ss$species)
  mean_psi <- mapply(function(a, s) mean(inp$surveys[[a]]$truth$psi[, s]), ss$arrangement, ss$species)
  record(id, "occupied_counts_both_arms", length(occupied), max(abs(occupied - ss$occupied)), TRUE)
  record(id, "detection_counts_both_arms", length(detections), max(abs(detections - ss$detections)), TRUE)
  record(id, "mean_true_probability", length(mean_psi), max(abs(mean_psi - ss$truth)))

  # Replay the generator's data stream through the two environmental normal
  # vectors and species-field normal matrix, then recover its occupancy uniforms.
  N <- nrow(inp$landscape$points); S <- length(inp$landscape$species)
  set.seed(sweep_seed("spatial-design-data", r))
  invisible(rnorm(N)); invisible(rnorm(N)); invisible(rnorm(N * S))
  u <- matrix(runif(N * S), N, S)
  z <- matrix(as.integer(u < inp$landscape$psi), N, S)
  stopifnot(identical(z, unname(inp$landscape$z)))
  record(id, "replayed_uniforms_match_regenerated_z", length(z), exact = TRUE)
  gap <- min(abs(u - inp$landscape$psi))
  stopifnot(is.finite(gap), gap > tolerance)
  replay[[r]] <- data.frame(community = id, n = length(u),
    minimum_uniform_threshold_distance = gap, comparison_tolerance = tolerance)
  if (r == 1L) {
    z_saved_psi <- matrix(as.integer(u < b$landscape$psi), N, S)
    stopifnot(identical(z_saved_psi, unname(inp$landscape$z)))
    record(id, "saved_original_psi_and_replayed_uniforms_match_z", length(z), exact = TRUE)
  }
}
cat("R: ", R.version.string, "\n", sep = "")
cat("RNG: ", paste(RNGkind(), collapse = ", "), "\n", sep = "")
cat("Generator MD5: ", unname(tools::md5sum(generator_path)), "\n", sep = "")
cat("Bundle MD5: ", unname(tools::md5sum(bundle_path)), "\n", sep = "")
cat("Original package revision: ", b$provenance$revision, "\n", sep = "")
write.csv(do.call(rbind, checks), stdout(), row.names = FALSE)
cat("Replay threshold margins (current stream; full historical probabilities survive only for rep01):\n")
write.csv(do.call(rbind, replay), stdout(), row.names = FALSE)
cat("PASS: reconstruction agrees with all checked surviving evidence.\n")
cat("LIMIT: original full inputs and posterior draws are absent; no original-file or full-response-matrix identity is certified.\n")
