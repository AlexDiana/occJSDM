# Export the actual package table from the existing fit, without fitting again.
# Rscript dev/simstudy/vignette-lesson/summarise_latent_tables.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")

lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
lesson <- readRDS(lesson_path)
manifest <- lesson$manifests$default
fit_path <- file.path(archive, manifest$file)
stopifnot(identical(unname(tools::md5sum(fit_path)), manifest$md5),
          identical(lesson$input$source_hashes, lesson_source_hashes()))

saved <- readRDS(fit_path)
stopifnot(identical(saved$source_hashes, lesson$input$source_hashes),
          identical(saved$input_md5, lesson$input_md5))
fit <- saved$fit
validate_lesson_fit_identity(fit, lesson$input)

# Preserve every returned row and value. PCR identifies row order within each
# declared sample/primer combination; it is not a new biological replicate.
tables <- map_dfr(seq_along(fit$infos$speciesNames), function(index) {
  occJSDM::returnLatentPresences(fit, idx_species = index) |>
    as_tibble() |>
    mutate(species = fit$infos$speciesNames[index], .before = 1) |>
    group_by(Site, Sample, Primer) |>
    mutate(PCR = row_number()) |>
    ungroup()
})

keys <- c("species", "Site", "Sample", "Primer", "PCR")
stopifnot(!anyDuplicated(tables[keys]))
identity_check <- tables |>
  left_join(select(lesson$observations, all_of(keys), reads),
            by = keys, relationship = "one-to-one")
stopifnot(nrow(identity_check) == nrow(lesson$observations),
          identical(identity_check$OTU, identity_check$reads))

bundle <- list(
  schema = 1L,
  tables = tables,
  fit_manifest = manifest,
  lesson_md5 = unname(tools::md5sum(lesson_path)),
  source_hashes = lesson$input$source_hashes,
  exporter_md5 = unname(tools::md5sum(
    "dev/simstudy/vignette-lesson/summarise_latent_tables.R"))
)
saveRDS(bundle, "vignettes/teaching-data/latent-presence-lesson.rds", compress = "xz")
cat("Exported", nrow(tables), "native table rows for", n_distinct(tables$species),
    "species. No MCMC was run.\n")
