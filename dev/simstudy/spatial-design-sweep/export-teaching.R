#!/usr/bin/env Rscript
# Writes vignettes/teaching-data/spatial-lesson.rds, the compact bundle Lesson 2 renders from.
# Run after plot.R: the tables are read from the committed compact results, the maps from the raw study.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) { hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }; substring(hit, nchar(name) + 4L) }
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
summary <- file.path(study, "summary-final"); results <- file.path(repo, "dev/simstudy/spatial-design-sweep/results")
read <- function(f) read.csv(file.path(results, f))
# The reason and warning columns must stay character: an all-empty column otherwise reads as logical NA,
# and the committed copy writes those empty values as NA, so both are normalised to "".
read_selected <- function(path) {
  text <- c("initial_reasons", "selected_reasons", "warnings")
  d <- read.csv(path, colClasses = setNames(rep("character", length(text)), text))
  for (v in text) d[[v]][is.na(d[[v]])] <- ""
  d
}
# The tables carried by the bundle are the committed compact copies; check they match the study's own.
for (f in c("groups.csv", "species.csv", "field.csv", "range.csv", "amplitude.csv", "lattice.csv",
            "aggregate.csv", "reading.csv", "range-reading.csv", "paired.csv", "oracle-lattice.csv"))
  stopifnot(unname(tools::md5sum(file.path(results, f))) == unname(tools::md5sum(file.path(summary, f))))
stopifnot(unname(tools::md5sum(file.path(results, "design-statistics.csv"))) == unname(tools::md5sum(file.path(study, "design-statistics.csv"))),
          unname(tools::md5sum(file.path(results, "audit.csv"))) == unname(tools::md5sum(file.path(study, "audit/audit.csv"))))
input1 <- readRDS(file.path(study, "inputs/rep01.rds")); land <- input1$landscape
arrangements <- do.call(rbind, lapply(1:3, function(r) { inp <- readRDS(file.path(study, "inputs", sprintf("rep%02d.rds", r)))
  do.call(rbind, lapply(names(inp$surveys), function(a) data.frame(community = sprintf("rep%02d", r), arrangement = a,
    site = seq_len(100), x = inp$surveys[[a]]$truth$xy[, 1], y = inp$surveys[[a]]$truth$xy[, 2]))) }))
sel <- read_selected(file.path(summary, "selected-fits.csv")); map_species <- c("species06", "species02")
stopifnot(nrow(sel) == 24L, all(unname(tools::md5sum(sel$selected_file)) == sel$selected_md5))
field_maps <- list(); lattice_maps <- list()
for (i in seq_len(nrow(sel))) {
  r <- readRDS(sel$selected_file[i]); inp <- readRDS(file.path(study, "inputs", paste0(r$job$community, ".rds")))
  stopifnot(identical(r$job$key, sel$key[i]))
  for (sp in map_species) { s <- match(sp, inp$landscape$species)
    field_maps[[length(field_maps) + 1L]] <- data.frame(community = r$job$community, arrangement = r$job$arrangement,
      source = r$job$arm, site = seq_len(100), species = sp, value = r$field_median[, s]) }
  if (r$job$community == "rep01") { li <- inp$landscape$index$lattice; s <- match("species06", inp$landscape$species)
    lattice_maps[[length(lattice_maps) + 1L]] <- data.frame(community = "rep01", arrangement = r$job$arrangement,
      source = r$job$arm, cell = seq_along(li), x = inp$landscape$points[li, 1], y = inp$landscape$points[li, 2],
      truth = inp$landscape$psi[li, s], with = unname(r$lattice_mean[, s]), without = unname(r$lattice_mean_nospatial[, s]),
      distance = r$lattice_distance) }
}
oracle <- read.csv(file.path(study, "oracle/oracle-selected.csv"))
for (i in seq_len(nrow(oracle))) { sp <- sprintf("species%02d", oracle$species[i]); if (!sp %in% map_species) next
  o <- readRDS(oracle$file[i])
  field_maps[[length(field_maps) + 1L]] <- data.frame(community = oracle$community[i], arrangement = oracle$arrangement[i],
    source = "oracle", site = seq_len(100), species = sp, value = o$field_median) }
for (r in 1:3) { inp <- readRDS(file.path(study, "inputs", sprintf("rep%02d.rds", r)))
  for (a in names(inp$surveys)) for (sp in map_species) { s <- match(sp, inp$landscape$species)
    field_maps[[length(field_maps) + 1L]] <- data.frame(community = sprintf("rep%02d", r), arrangement = a, source = "truth",
      site = seq_len(100), species = sp, value = unname(inp$surveys[[a]]$truth$field[, s])) } }
named_md5 <- function(files) setNames(unname(tools::md5sum(files)), basename(files))
bundle <- list(
  landscape = list(points = land$points, index = land$index, environment = land$environment, field = land$field, psi = land$psi,
                   range = land$range, env_range = land$env_range, species = land$species, prevalence = land$prevalence,
                   B0 = land$B0, B = land$B, cluster_centres = land$cluster_centres),
  arrangements = arrangements, statistics = read("design-statistics.csv"),
  oracle = read("oracle.csv"), oracle_lattice = read("oracle-lattice.csv"), fits = list(groups = read("groups.csv"), species = read("species.csv"), field = read("field.csv"),
    range = read("range.csv"), amplitude = read("amplitude.csv"), lattice = read("lattice.csv")),
  aggregate = read("aggregate.csv"), reading = read("reading.csv"), range_reading = read("range-reading.csv"),
  paired = read("paired.csv"), field_maps = do.call(rbind, field_maps), lattice_maps = do.call(rbind, lattice_maps),
  selected_fits = read_selected(file.path(results, "selected-fits.csv")), audit = read("audit.csv"),
  provenance = list(revision = readLines(file.path(study, "source-revision.txt")),
    compact_hashes = named_md5(list.files(results, full.names = TRUE)),
    script_hashes = named_md5(list.files(file.path(repo, "dev/simstudy/spatial-design-sweep"), pattern = "\\.R$", full.names = TRUE)),
    exported = Sys.time()))
out <- file.path(repo, "vignettes/teaching-data/spatial-lesson.rds")
saveRDS(bundle, out, compress = "xz"); cat("Saved", out, round(file.size(out) / 1e6, 2), "MB\n")
