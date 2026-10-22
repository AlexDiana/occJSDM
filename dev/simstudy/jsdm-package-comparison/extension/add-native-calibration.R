# Add native coefficient intervals to a verified saved-fit calibration bundle.
# Rscript add-native-calibration.R BASE_BUNDLE NATIVE_OUTPUT CODE OUTPUT_BUNDLE
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
base <- normalizePath(args[1]); root <- normalizePath(args[2]); code <- normalizePath(args[3])
source(file.path(code, "calibration-math.R"))
bundle <- readRDS(base)
stopifnot(bundle$complete, is.null(bundle$native_intervals))
jobs <- subset(bundle$manifest, scored & package %in% c("gllvm", "sjSDM") & scenario != "traits" &
                 !(scenario == "curved" & response == "linear"))
records <- vector("list", nrow(jobs)); paths <- character(nrow(jobs)); dependencies <- character(nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i, ]; paths[i] <- file.path(root, "jobs", paste0(j$job, ".rds"))
  r <- readRDS(paths[i])
  dependencies[i] <- file.path(dirname(dirname(names(r$source_hashes)[1])), "pilot-math.R")
  stopifnot(identical(r$job, j), !r$truth_used,
    identical(tools::md5sum(names(r$source_hashes)), r$source_hashes),
    identical(tools::md5sum(names(r$input_hashes)), r$input_hashes))
  idx <- which(bundle$rows$job == j$job & bundle$rows$target == "Environmental coefficient")
  stopifnot(length(idx) > 0, all(is.na(bundle$rows$lower[idx])), all(is.na(bundle$rows$upper[idx])))
  if (r$ok) {
    stopifnot(length(idx) == nrow(r$intervals),
      max(abs(bundle$rows$estimate[idx] - r$intervals$estimate)) < 1e-10)
    bundle$rows$lower[idx] <- r$intervals$lower
    bundle$rows$upper[idx] <- r$intervals$upper
  }
  records[[i]] <- data.frame(job=j$job, package=j$package, ok=r$ok,
    original_diagnostic_pass=j$diagnostic_pass, numerical_pass=r$numerical_pass,
    intervals=if (r$ok) sum(is.finite(r$intervals$lower) & is.finite(r$intervals$upper)) else 0L,
    sampling=if (!is.null(r$audit$sampling)) r$audit$sampling else NA_integer_,
    se_change=if (!is.null(r$audit$checks)) tail(r$audit$checks$max_relative_se_change, 1) else NA_real_,
    parameter_change=if (!is.null(r$audit$parameter_change)) r$audit$parameter_change else NA_real_,
    method=if (!is.null(r$method)) r$method else NA_character_,
    error=if (!is.null(r$error)) r$error else NA_character_,
    warnings=paste(r$warnings, collapse="; "))
}
bundle$summary <- calibration_group_summary(bundle$rows)
bundle$per_community <- do.call(rbind, lapply(split(bundle$rows,
  interaction(bundle$rows$job, bundle$rows$target, bundle$rows$term, drop=TRUE)), function(d) {
    cbind(d[1, c("package", "scenario", "n_sites", "n_species", "response", "use_traits", "target", "term", "replicate", "diagnostic_pass")], calibration_summary(d))
  }))
bundle$native_intervals <- do.call(rbind, records)
passed_native <- with(bundle$native_intervals, ok & numerical_pass & original_diagnostic_pass)
native_rows <- subset(bundle$rows, target == "Environmental coefficient" &
  job %in% bundle$native_intervals$job[passed_native])
bundle$native_passed_summary <- subset(calibration_group_summary(native_rows), population == "All scored fits")
bundle$native_passed_summary$population <- "Passed fit and native interval checks"
bundle$native_provenance <- list(base=tools::md5sum(base), jobs=tools::md5sum(paths),
  worker_dependencies=tools::md5sum(unique(dependencies)),
  sources=tools::md5sum(file.path(code, c("add-native-calibration.R", "calibration-math.R"))),
  exported_at=Sys.time(), session=sessionInfo())
saveRDS(bundle, args[4], compress="xz")
write.csv(bundle$native_intervals, file.path(dirname(args[4]), "native-interval-status.csv"), row.names=FALSE)
write.csv(bundle$summary, file.path(dirname(args[4]), "interval-summary.csv"), row.names=FALSE)
cat("Added native coefficient interval records for", nrow(jobs), "fits. Verify before installing the bundle.\n")
