#!/usr/bin/env Rscript
# Checks the Lesson 2 bundle against the committed compact results. Needs no raw archive.
repo <- normalizePath(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ".")
results <- file.path(repo, "dev/simstudy/spatial-design-sweep/results")
b <- readRDS(file.path(repo, "vignettes/teaching-data/spatial-lesson.rds"))
# Same reader as export-teaching.R: the reason and warning columns stay character, with NA read as "".
read_selected <- function(path) {
  text <- c("initial_reasons", "selected_reasons", "warnings")
  d <- read.csv(path, colClasses = setNames(rep("character", length(text)), text))
  for (v in text) d[[v]][is.na(d[[v]])] <- ""
  d
}
files <- list.files(results, full.names = TRUE)
stopifnot(identical(setNames(unname(tools::md5sum(files)), basename(files)), b$provenance$compact_hashes))
for (t in c("groups", "species", "field", "range", "amplitude", "lattice"))
  stopifnot(identical(b$fits[[t]], read.csv(file.path(results, paste0(t, ".csv")))))
for (t in list(c("aggregate", "aggregate.csv"), c("reading", "reading.csv"), c("range_reading", "range-reading.csv"),
               c("paired", "paired.csv"), c("oracle", "oracle.csv"), c("statistics", "design-statistics.csv"),
               c("audit", "audit.csv")))
  stopifnot(identical(b[[t[1]]], read.csv(file.path(results, t[2]))))
stopifnot(identical(b$selected_fits, read_selected(file.path(results, "selected-fits.csv"))))
f <- b$fits$field; f$group <- paste0("prevalence_", f$target * 100, "pct"); f <- f[f$arm == "binary", ]
cells <- aggregate(cbind(centred_rmse, centred_correlation, zero_field_rmse) ~ community + arrangement + group, f, mean)
cells$reduction <- 1 - cells$centred_rmse / cells$zero_field_rmse
for (i in seq_len(nrow(b$reading))) { d <- cells[cells$arrangement == b$reading$arrangement[i] & cells$group == b$reading$group[i], ]
  stopifnot(nrow(d) == 3L)
  label <- if (all(d$reduction >= .2 & d$centred_correlation >= .5)) "informative" else
           if (all(d$reduction < .1 | d$centred_correlation < .3)) "uninformative" else "intermediate"
  stopifnot(identical(label, b$reading$label[i])) }
stopifnot(all(abs(colMeans(b$landscape$psi[b$landscape$index$lattice, ]) - b$landscape$prevalence) < 1e-6),
          all(b$audit$passed), nrow(b$selected_fits) == 24L, sum(b$selected_fits$needs_long) == 3L,
          !any(nzchar(b$selected_fits$selected_reasons)), !any(nzchar(b$selected_fits$warnings)))
# The maps' truth must be the landscape's own field at the surveyed sites, and the lattice map's truth its psi.
truth <- b$field_maps[b$field_maps$source == "truth" & b$field_maps$community == "rep01", ]
lat <- b$landscape$index$lattice
stopifnot(nrow(b$lattice_maps) == 8L * length(lat),
          all(abs(b$lattice_maps$truth - rep(b$landscape$psi[lat, "species06"], 8L)) < 1e-12),
          nrow(truth) == 800L, all(table(b$field_maps$source) == 2400L))
cat("Lesson 2 bundle verified against the committed results.\n")
