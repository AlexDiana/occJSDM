library(testthat)
repo <- normalizePath("../../../..")
if (file.exists("../config.R")) source("../config.R")
if (file.exists("../archive.R")) source("../archive.R")

test_that("both explicit input modes preserve the complete roster and seeds", {
  cfg <- study_config()
  expect_identical(names(cfg$arms), c("occJSDM", "sjSDM", "Hmsc", "gllvm"))
  for (mode in c("restored", "reconstructed")) {
    m <- comparison_manifest(cfg, mode)
    expect_equal(nrow(m), 192L)
    expect_false(anyDuplicated(m$key) > 0)
    expect_equal(sum(m$arm == "gp"), 12L)
    expect_true(all(m$new_fit))
    expect_equal(m$fit_seed, 27080000L + 100L * m$job_id + 1L)
    expect_equal(nrow(oracle_manifest(cfg, mode)), 96L)
  }
  expect_error(comparison_manifest(cfg, "automatic"), "input mode")
  expect_equal(cfg$remaining_wall_hours, 48)
})

test_that("missing archives cannot silently regenerate", {
  expect_error(verify_archive(tempfile(), repo), "archive prerequisite")
})

test_that("relative manifest paths survive moving roots and reject ambiguity or corruption", {
  root <- tempfile(); dir.create(root)
  dir.create(file.path(root, "inputs"))
  p <- file.path(root, "inputs/rep01.rds"); saveRDS(list(a=1), p)
  manifest <- data.frame(identity="rep01", path="/deleted/inputs/rep01.rds",
                         relative_path="inputs/rep01.rds", md5=unname(tools::md5sum(p)))
  expect_identical(resolve_manifest(manifest, root)$path, normalizePath(p))
  expect_error(resolve_manifest(rbind(manifest, manifest), root), "duplicate")
  saveRDS(list(a=2), p)
  expect_error(resolve_manifest(manifest, root), "hash mismatch")
  manifest$relative_path <- "../outside.rds"
  expect_error(resolve_manifest(manifest, root), "relative path")
})

test_that("audit rejects changed seeds, counts, values or generator hashes", {
  source(file.path(repo, "dev/simstudy/spatial-design-sweep/generator.R"), local=TRUE)
  inputs <- lapply(1:3, make_sweep_input)
  evidence <- audit_inputs(inputs, repo)
  expect_true(all(evidence$checks$passed))
  expect_equal(nrow(evidence$checks), 61L)
  expect_true(all(evidence$replay$minimum_uniform_threshold_distance > 1e-12))
  bad <- inputs; bad[[2]]$settings$seeds[1] <- 1L
  expect_error(audit_inputs(bad, repo), "seed")
  bad <- inputs; bad[[1]]$surveys$spread$binary$OTU[1,1] <- 1 - bad[[1]]$surveys$spread$binary$OTU[1,1]
  expect_error(audit_inputs(bad, repo), "detection_counts")
  bad <- inputs; bad[[2]]$surveys$spread$truth$field[1,2] <- bad[[2]]$surveys$spread$truth$field[1,2] + 1e-5
  expect_error(audit_inputs(bad, repo), "saved_two_species_fields")
  expect_error(audit_inputs(inputs, repo, generator_md5="wrong"), "generator hash")
  expect_error(audit_inputs(inputs, repo, tolerance=1e-5), "tolerance")
})

test_that("persisted reconstruction imports ordered fitter-only files without truth", {
  root <- tempfile(); dir.create(root)
  record <- reconstruct_inputs(repo, root)
  expect_identical(record$mode, "reconstructed")
  expect_equal(length(record$files$path), 3)
  m <- import_inputs(record, root)
  expect_equal(nrow(m), 36L)
  expect_false(anyDuplicated(m$relative_path) > 0)
  a <- readRDS(file.path(root, m$relative_path[m$identity == "rep01-spread-training"]))
  b <- readRDS(file.path(root, m$relative_path[m$identity == "rep01-spread-prediction"]))
  inp <- readRDS(record$files$path[1])
  expect_identical(a$y, inp$surveys$spread$binary$OTU)
  expect_identical(a$species_id, inp$landscape$species)
  expect_identical(a$site_id, sprintf("rep01-point%04d", 1:100))
  expect_setequal(names(a), c("schema_version", "dataset_id", "site_id", "species_id", "y", "environment", "xy", "input_hashes"))
  expect_setequal(names(b), c("schema_version", "dataset_id", "site_id", "species_id", "environment", "xy", "input_hashes"))
  expect_error(reconstruct_inputs(repo, root), "already exists")
  expect_identical(readRDS(file.path(root,"input-record.rds"))$mode,"reconstructed")
})
