# Records the two earlier readings of the map and how far the post hoc changes moved the point estimates.
# Reads the ignored study directory's earlier score sets and writes only new files into results/earlier-readings/.
#   first-run-*       scores-joint-sampler/ : the first truth 1 run, with the joint sampler
#   pre-resolution-*  scores-before-resolution/ : the truth 1 re-run with the blockwise sampler, before the resolution extension
#   agreement.csv     per cell, how much gate-species slope and contrast correlation changed between an earlier score set and the final one
# run-r -f earlier-readings.R --args --study=DIR
args <- commandArgs(trailingOnly = TRUE)
study <- sub("^--study=", "", grep("^--study=", args, value = TRUE))
res <- "dev/simstudy/spatial-recoverability-pilot/results"
out <- file.path(res, "earlier-readings"); dir.create(out, showWarnings = FALSE)
for (pair in list(c("scores-joint-sampler", "first-run"), c("scores-before-resolution", "pre-resolution")))
  for (f in c("gate.csv", "map.json")) file.copy(file.path(study, pair[1], f), file.path(out, paste0(pair[2], "-", f)), overwrite = TRUE)

final <- read.csv(file.path(res, "species.csv"), stringsAsFactors = FALSE)
key <- c("cell", "replicate", "arrangement", "species")
compare <- function(dir, label, truth_only = NULL, changed_only = FALSE) {
  old <- read.csv(file.path(study, dir, "species.csv"), stringsAsFactors = FALSE)
  m <- merge(old[c(key, "field_slope", "contrast_correlation")], final[c(key, "field_slope", "contrast_correlation", "truth", "occupancy")], by = key, suffixes = c("_old", "_new"))
  m$d_slope <- abs(m$field_slope_new - m$field_slope_old); m$d_corr <- abs(m$contrast_correlation_new - m$contrast_correlation_old)
  stopifnot(nrow(m) == nrow(final))
  if (!is.null(truth_only)) m <- m[m$truth == truth_only, ]
  if (changed_only) m <- m[m$d_slope > 0 | m$d_corr > 0, ]
  m <- m[m$occupancy %in% c(0.25, 0.5, 0.75), ]
  row <- function(d, cell) data.frame(comparison = label, cell = cell, gate_species_compared = nrow(d),
    median_abs_slope_change = median(d$d_slope), max_abs_slope_change = max(d$d_slope),
    median_abs_correlation_change = median(d$d_corr), max_abs_correlation_change = max(d$d_corr))
  rbind(do.call(rbind, lapply(split(m, m$cell), function(d) row(d, d$cell[1]))), row(m, "all"))
}
agreement <- rbind(
  compare("scores-joint-sampler", "first run (joint sampler) against final, truth 1", truth_only = "t1"),
  compare("scores-before-resolution", "pre-resolution against final, species whose fit was replaced", changed_only = TRUE))
rownames(agreement) <- NULL
write.csv(agreement, file.path(out, "agreement.csv"), row.names = FALSE)
print(agreement)
cat("Wrote", out, "\n")
