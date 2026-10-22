#!/usr/bin/env Rscript
# Describe all selected fits, retaining unresolved flags and initial sensitivity.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  x <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(x)) sub(paste0("^--", name, "="), "", x[1]) else default
}
repo <- normalizePath(arg("repo", ".")); study <- normalizePath(arg("study"))
out <- file.path(repo, "dev/simstudy/current-main-recheck/results")
selection <- readRDS(file.path(study, "selection.rds"))
stopifnot(nrow(selection) == 110L, !anyDuplicated(selection$key))
rows <- list(); flags <- list(); warnings <- list()
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
for (i in seq_len(nrow(selection))) {
  s <- selection[i, ]; r <- readRDS(s$result_file)
  g <- r$groups
  primary <- if (s$family == "jsdm") g[g$scope == "original100", ] else g[g$metric == "occupancy_original_sites", ]
  stopifnot(nrow(primary) == 4L)
  rows[[s$key]] <- data.frame(s[c("key", "family", "arm", "scenario", "replicate", "schedule")],
    warnings = length(r$warnings), max_primary_rhat = safe_max(primary$rhat),
    max_group_rhat = r$diagnostics$max_group_rhat, max_element_rhat = r$diagnostics$max_element_rhat,
    unresolved_rhat = r$diagnostics$unresolved_rhat,
    min_primary_ess_mean = safe_min(primary$ess_mean), max_primary_mcse = safe_max(primary$mcse),
    max_primary_chain_gap = safe_max(primary$chain_gap),
    fit_seconds = as.numeric(difftime(r$finished, r$started, units = "secs")))
  for (block in c("groups", "elements", "additional_diagnostics")) {
    z <- r[[block]]
    keep <- !is.finite(z$rhat) | z$rhat <= 0 | z$rhat > 1.05
    z <- z[keep, , drop = FALSE]
    if (!nrow(z)) next
    if (!"metric" %in% names(z)) z$metric <- z$scope
    if (!"group" %in% names(z)) z$group <- NA_character_
    if (!"element" %in% names(z)) z$element <- NA_integer_
    flags[[paste(s$key, block)]] <- data.frame(key = s$key, block = block,
      z[c("metric", "group", "element", "rhat", "ess_mean", "mcse", "chain_gap")])
  }
  if (length(r$warnings)) warnings[[s$key]] <- data.frame(key = s$key, warning = r$warnings)
}
diagnostics <- do.call(rbind, rows)
flagged <- if (length(flags)) do.call(rbind, flags) else data.frame(key=character(), block=character(), metric=character(), group=character(), element=integer(), rhat=numeric(), ess_mean=numeric(), mcse=numeric(), chain_gap=numeric())
warning_rows <- if (length(warnings)) do.call(rbind, warnings) else data.frame(key=character(), warning=character())
write.csv(diagnostics, file.path(out, "diagnostic-summary.csv"), row.names = FALSE)
write.csv(flagged, file.path(out, "flagged-diagnostics.csv"), row.names = FALSE)
write.csv(warning_rows, file.path(out, "warnings.csv"), row.names = FALSE)

scores <- read.csv(file.path(out, "community-scores.csv"))
primary <- subset(scores, metric == "occupancy_original_sites")
keys <- c("key", "family", "arm", "scenario", "replicate", "group")
sensitivity <- merge(subset(primary, version == "current_initial"), subset(primary, version == "current"),
                     by = keys, suffixes = c("_initial", "_selected"))
stopifnot(nrow(sensitivity) == 440L)
for (m in c("bias", "mae")) sensitivity[[paste0("change_", m)]] <- sensitivity[[paste0(m, "_selected")]] - sensitivity[[paste0(m, "_initial")]]
write.csv(sensitivity, file.path(out, "long-run-community-sensitivity.csv"), row.names = FALSE)
changes <- aggregate(cbind(change_bias, change_mae) ~ family + arm + scenario + group, sensitivity, mean)
write.csv(changes, file.path(out, "long-run-summary-sensitivity.csv"), row.names = FALSE)
cat("Selected fits:", nrow(diagnostics), "; longer fits:", sum(diagnostics$schedule == "long"), "\n")
cat("Fitting warnings:", nrow(warning_rows), "; diagnostic rows above screen or undefined:", nrow(flagged), "\n")
cat("Maximum primary/group/element Rhat:", max(diagnostics$max_primary_rhat), max(diagnostics$max_group_rhat), max(diagnostics$max_element_rhat), "\n")
cat("Largest absolute ten-community sensitivity, bias/MAE:", max(abs(changes$change_bias)), max(abs(changes$change_mae)), "\n")
print(diagnostics[diagnostics$warnings > 0 | diagnostics$max_group_rhat > 1.05 | diagnostics$max_element_rhat > 1.05 | diagnostics$unresolved_rhat > 0, ], row.names = FALSE)
