# Bridge from the pilot's community objects to the upstream TrainingData and PredictionData
# contracts, plus the oracle reference. Truth never enters TrainingData or PredictionData.
stage2_datasets <- function(pilot_study, config = stage2_config()) {
  man <- read.csv(file.path(pilot_study, "input-manifest.csv"), stringsAsFactors = FALSE)
  rows <- list()
  for (cell in config$cells) for (r in config$replicates) {
    file <- sprintf("%s-rep%d.rds", cell, r)
    m <- man[man$file == file, ]; stopifnot(nrow(m) == 1L)
    job <- sprintf("%s-rep%d-%s", cell, r, config$arrangement)
    resolved <- file.path(pilot_study, "jobs-resolution", paste0(job, ".rds"))
    oracle <- if (file.exists(resolved)) resolved else file.path(pilot_study, "jobs", paste0(job, ".rds"))
    stopifnot(file.exists(oracle))
    rows[[length(rows) + 1L]] <- data.frame(dataset_id = sprintf("%s-rep%d", cell, r), cell = cell, truth = substr(cell, 1, 2),
      n = m$n, replicate = r, input_file = file.path(pilot_study, "inputs", file), input_md5 = m$md5,
      oracle_file = oracle, oracle_md5 = unname(tools::md5sum(oracle)), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

site_ids <- function(n) sprintf("site%04d", seq_len(n))
lattice_ids <- function(n) sprintf("lattice%04d", seq_len(n))

load_training <- function(cm, dataset_id, input_md5) {
  idx <- cm$index$uniform
  xy <- unname(cm$points[idx, , drop = FALSE]); dimnames(xy) <- NULL
  y <- cm$z[idx, , drop = FALSE]; colnames(y) <- cm$species$species
  list(schema_version = 1L, dataset_id = dataset_id, site_id = site_ids(length(idx)), species_id = cm$species$species,
       y = y, environment = matrix(cm$environment[idx], ncol = 1L), xy = xy, input_hashes = c(input = input_md5))
}

load_prediction <- function(cm) {
  idx <- cm$index$lattice
  xy <- unname(cm$points[idx, , drop = FALSE]); dimnames(xy) <- NULL
  list(schema_version = 1L, site_id = lattice_ids(length(idx)), species_id = cm$species$species,
       environment = matrix(cm$environment[idx], ncol = 1L), xy = xy)
}

load_truth <- function(cm) {
  lat <- cm$index$lattice; idx <- cm$index$uniform
  list(psi_lattice = cm$psi[lat, , drop = FALSE], field_lattice = cm$spatial[lat, , drop = FALSE],
       psi_training = cm$psi[idx, , drop = FALSE], field_training = cm$spatial[idx, , drop = FALSE],
       eta_lattice = sweep(cm$environment[lat] %o% cm$B, 2, cm$B0, "+"),
       occupancy = cm$species$occupancy, species = cm$species$species, spatial_sd = cm$spatial_sd)
}

load_oracle_reference <- function(path) {
  r <- readRDS(path)
  list(probability_mean = r$prediction$probability_mean, field_mean = r$prediction$field_mean,
       screen = r$screen, status = if (all(r$screen$status == "pass")) "pass" else "unresolved",
       phases = r$phases, md5 = unname(tools::md5sum(path)), job = r$job)
}
