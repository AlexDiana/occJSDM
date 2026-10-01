#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", name)
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }
  substring(hit, nchar(name) + 4L)
}
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
mode <- option("mode", "select"); stopifnot(mode %in% c("select", "final"))

diagnostic_reasons <- function(a) {
  reasons <- character()
  if (any(a$groups$rhat > 1.05, na.rm = TRUE)) reasons <- c(reasons, "group Rhat > 1.05")
  if (any(a$species$rhat > 1.05, na.rm = TRUE)) reasons <- c(reasons, "species Rhat > 1.05")
  if (any(a$elements$rhat > 1.05, na.rm = TRUE)) reasons <- c(reasons, "element Rhat > 1.05")
  primary <- a$groups$metric == "occupancy"
  if (any(a$groups$ess_mean[primary] < 100, na.rm = TRUE)) reasons <- c(reasons, "occupancy group ESS < 100")
  if (any(a$species$ess_mean < 100, na.rm = TRUE)) reasons <- c(reasons, "occupancy species ESS < 100")
  if (any(!is.finite(a$elements$rhat[a$elements$metric != "range"]))) reasons <- c(reasons, "non-range element Rhat unavailable")
  if (any(grepl("Convergence|Low ESS", a$warnings))) reasons <- c(reasons, "native convergence warning")
  reasons
}

read_sweep_selection <- function(study, final = FALSE) {
  manifest <- readRDS(file.path(study, "input-manifest.rds")); stopifnot(nrow(manifest) == 24L)
  selected <- list(); initial <- list(); rows <- list()
  for (key in manifest$key) {
    f <- file.path(study, "initial", paste0(key, "-result.rds"))
    if (!file.exists(f)) stop("Initial batch incomplete: ", key)
    a <- readRDS(f); stopifnot(identical(a$mcmc, list(nchain = 2L, nburn = 3000L, niter = 5000L, nthin = 1L)))
    initial[[key]] <- a; reasons <- diagnostic_reasons(a); needs_long <- length(reasons) > 0L
    b <- a; phase <- "initial"; selected_file <- f
    if (final && needs_long) {
      lf <- file.path(study, "long", paste0(key, "-result.rds"))
      if (!file.exists(lf)) stop("Prespecified longer check missing: ", key)
      b <- readRDS(lf); phase <- "long"; selected_file <- lf
      stopifnot(identical(b$mcmc, list(nchain = 4L, nburn = 6000L, niter = 12000L, nthin = 1L)), identical(b$job, a$job))
    }
    selected[[key]] <- b
    rows[[key]] <- data.frame(key = key, community = a$job$community, arrangement = a$job$arrangement, arm = a$job$arm,
      replicate = a$job$replicate, needs_long = needs_long, initial_reasons = paste(reasons, collapse = "; "),
      phase = phase, selected_file = selected_file, selected_md5 = unname(tools::md5sum(selected_file)),
      selected_reasons = paste(diagnostic_reasons(b), collapse = "; "),
      max_group_rhat = max(b$groups$rhat, na.rm = TRUE), warnings = paste(b$warnings, collapse = " | "))
  }
  list(manifest = do.call(rbind, rows), selected = selected, initial = initial)
}

collect <- function(results, table) do.call(rbind, lapply(names(results), function(key) {
  r <- results[[key]]; d <- r[[table]]; if (is.null(d)) return(NULL)
  data.frame(key = key, community = r$job$community, arrangement = r$job$arrangement, arm = r$job$arm,
             replicate = r$job$replicate, phase = if (identical(r$mcmc$nchain, 4L)) "long" else "initial", d) }))

reading_labels <- function(field_rows) {
  field_rows$reduction <- 1 - field_rows$centred_rmse / field_rows$zero_field_rmse
  cells <- split(field_rows, interaction(field_rows$arrangement, field_rows$group, drop = TRUE))
  do.call(rbind, lapply(cells, function(d) {
    stopifnot(nrow(d) == 3L)
    label <- if (all(d$reduction >= .2 & d$centred_correlation >= .5)) "informative" else
             if (all(d$reduction < .1 | d$centred_correlation < .3)) "uninformative" else "intermediate"
    data.frame(arrangement = d$arrangement[1], group = d$group[1], label = label,
               min_reduction = min(d$reduction), min_correlation = min(d$centred_correlation)) }))
}

selection <- read_sweep_selection(study, final = mode == "final")
if (mode == "select") {
  write.csv(selection$manifest, file.path(study, "long-selection.csv"), row.names = FALSE)
  keys <- selection$manifest$key[selection$manifest$needs_long]
  writeLines(keys, file.path(study, "long-keys.txt"))
  cat(length(keys), "of 24 fits require the prespecified longer checks.\n"); quit(status = 0)
}
out <- file.path(study, "summary-final")
if (dir.exists(out) && length(list.files(out))) stop("Use an empty output directory: ", out)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tables <- list(groups = collect(selection$selected, "groups"), species = collect(selection$selected, "species"),
               field = collect(selection$selected, "field"), range = collect(selection$selected, "range_summary"),
               amplitude = collect(selection$selected, "amplitude"), lattice = collect(selection$selected, "lattice"))
for (t in names(tables)) write.csv(tables[[t]], file.path(out, paste0(t, ".csv")), row.names = FALSE)
oracle <- read.csv(file.path(study, "oracle/oracle-selected.csv"))
oracle$group <- paste0("prevalence_", c(5, 5, 25, 25, 25, 75, 75, 75)[oracle$species], "pct")
write.csv(oracle, file.path(out, "oracle.csv"), row.names = FALSE)
# Oracle lattice prediction (outcome 4 for the oracle), tabulated from the stored oracle results with the
# distance bins of score.R. "with" is the oracle's posterior mean probability; "without" is the
# environment-only probability from the true parameters, plogis(B0 + B * environment).
source(file.path(repo, "dev/simstudy/spatial-design-sweep/score.R"))
error_rows <- function(err, bins) rbind(
  data.frame(bin = "all", n = length(err), bias = mean(err), mae = mean(abs(err)), rmse = sqrt(mean(err^2))),
  do.call(rbind, lapply(levels(bins$bin), function(b) { e <- err[bins$bin == b]
    data.frame(bin = b, n = length(e), bias = mean(e), mae = mean(abs(e)), rmse = sqrt(mean(e^2))) })))
oracle_lattice <- do.call(rbind, lapply(seq_len(nrow(oracle)), function(i) {
  o <- oracle[i, ]; res <- readRDS(o$file)
  stopifnot(unname(tools::md5sum(o$file)) == o$result_md5, res$job$community == o$community,
            res$job$arrangement == o$arrangement, res$job$species == o$species)
  input_file <- file.path(study, "inputs", paste0(o$community, ".rds"))
  stopifnot(unname(tools::md5sum(input_file)) == res$input_md5)
  input <- readRDS(input_file); land <- input$landscape; li <- land$index$lattice; s <- o$species
  sites <- input$surveys[[o$arrangement]]$truth$xy
  bins <- distance_bins(land$points[li, ], sites); truth <- land$psi[li, s]
  with <- res$lattice$probability_mean; without <- plogis(land$B0[s] + land$B[s] * land$environment[li])
  stopifnot(length(with) == length(li), all(is.finite(with)))
  data.frame(community = o$community, arrangement = o$arrangement, species = s, group = o$group,
             rbind(data.frame(spatial_term = "with", error_rows(with - truth, bins)),
                   data.frame(spatial_term = "without", error_rows(without - truth, bins))))
}))
write.csv(oracle_lattice, file.path(out, "oracle-lattice.csv"), row.names = FALSE)
field <- tables$field; field$group <- paste0("prevalence_", field$target * 100, "pct")
field_cells <- aggregate(cbind(centred_rmse, centred_correlation, zero_field_rmse) ~ community + arrangement + arm + group, field, mean)
oracle_cells <- aggregate(cbind(centred_rmse, centred_correlation, zero_field_rmse) ~ community + arrangement + group, oracle, mean)
oracle_cells$arm <- "oracle"
cells <- rbind(field_cells, oracle_cells[names(field_cells)])
aggregate_cells <- do.call(rbind, lapply(split(cells, interaction(cells$arrangement, cells$arm, cells$group, drop = TRUE)), function(d)
  data.frame(arrangement = d$arrangement[1], arm = d$arm[1], group = d$group[1], n = nrow(d),
             centred_rmse = mean(d$centred_rmse), centred_rmse_min = min(d$centred_rmse), centred_rmse_max = max(d$centred_rmse),
             correlation = mean(d$centred_correlation), correlation_min = min(d$centred_correlation), correlation_max = max(d$centred_correlation),
             zero_field_rmse = mean(d$zero_field_rmse))))
write.csv(aggregate_cells, file.path(out, "aggregate.csv"), row.names = FALSE)
reading <- reading_labels(field_cells[field_cells$arm == "binary", ])
range_reading <- aggregate(mass_within_one_step ~ arrangement + arm, tables$range, function(v) all(v >= .5))
names(range_reading)[3] <- "range_recovered"
write.csv(reading, file.path(out, "reading.csv"), row.names = FALSE)
write.csv(range_reading, file.path(out, "range-reading.csv"), row.names = FALSE)
pairs <- merge(cells[cells$arm == "oracle", ], cells[cells$arm == "binary", ], by = c("community", "arrangement", "group"), suffixes = c("_oracle", "_binary"))
pairs <- merge(pairs, cells[cells$arm == "two_stage", ], by = c("community", "arrangement", "group"))
paired <- data.frame(pairs[c("community", "arrangement", "group")],
                     estimation_cost = pairs$centred_rmse_binary - pairs$centred_rmse_oracle,
                     detection_cost = pairs$centred_rmse - pairs$centred_rmse_binary)
write.csv(paired, file.path(out, "paired.csv"), row.names = FALSE)
write.csv(selection$manifest, file.path(out, "selected-fits.csv"), row.names = FALSE)
initial_groups <- collect(selection$initial, "groups")
sens <- merge(initial_groups, tables$groups, by = c("key", "metric", "group"), suffixes = c("_initial", "_selected"))
sens$bias_change <- sens$bias_selected - sens$bias_initial; sens$mae_change <- sens$mae_selected - sens$mae_initial
write.csv(sens[c("key", "metric", "group", "bias_change", "mae_change")], file.path(out, "long-run-sensitivity.csv"), row.names = FALSE)
saveRDS(selection$manifest, file.path(out, "selection.rds"))
saveRDS(list(scripts = tools::md5sum(file.path(repo, "dev/simstudy/spatial-design-sweep", c("summarise.R", "score.R", "generator.R", "oracle.R"))),
             session = sessionInfo(), time = Sys.time()), file.path(out, "provenance.rds"))
cat("Wrote final summary to", out, "\n")
