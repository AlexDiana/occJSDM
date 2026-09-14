# These fixtures describe samples in the sampler's known site/sample order.
# Expected covariates come from that table, not from the grouping helper.
alignment_samples <- function(ids = "global") {
  samples <- data.frame(
    Site = c(1, 1, 2, 2, 10, 10),
    Sample = 1:6,
    effort = c(-3, 4, -1, 6, 2, 0),
    collection_type = c("dry", "wet", "wet", "dry", "dry", "wet")
  )
  if (ids == "local") samples$Sample <- rep(1:2, 3)
  if (ids == "character") {
    samples$Site <- rep(c("a", "a-b", "c"), each = 2)
    samples$Sample <- c("b-c", "d", "c", "d", "x", "y")
    # Samples 1 and 3 both paste to "a-b-c" but are distinct samples.
  }
  samples
}

alignment_data <- function(samples, model = "two_stage", shuffled = FALSE) {
  if (model == "two_stage") {
    # Unequal primer sets and PCR replication exercise P/K and idx_w_k too.
    block_sample <- c(1, 1, 2, 3, 3, 4, 5, 6, 6)
    block_primer <- c("A", "B", "B", "A", "B", "A", "A", "A", "B")
    block_k <- c(2, 1, 3, 1, 2, 2, 1, 2, 1)
    sample_row <- rep(block_sample, block_k)
    info <- samples[sample_row, , drop = FALSE]
    info$Primer <- rep(block_primer, block_k)
  } else {
    sample_row <- 1:6
    info <- samples
  }
  otu <- cbind(sp1 = c(5, 0, 8, 2, 0, 3)[sample_row],
               sp2 = c(0, 2, 0, 5, 3, 0)[sample_row])
  if (shuffled) {
    # Change group order and within-group order without drawing random values.
    perm <- rev(seq_len(nrow(info)))
    info <- info[perm, , drop = FALSE]
    otu <- otu[perm, , drop = FALSE]
  }
  list(info = info, OTU = otu)
}

fit_alignment_data <- function(data, covariates = c("effort", "collection_type")) {
  suppressMessages(suppressWarnings(runOccJSDM(
    data, collCovariates = covariates,
    MCMCparams = list(nchain = 2, nburn = 10, niter = 10, nthin = 1)
  )))
}

expect_sample_covariates <- function(fit, samples) {
  # Standardisation is calculated independently in this known sample order.
  expected_effort <- (samples$effort - mean(samples$effort)) / sd(samples$effort)
  expect_equal(unname(fit$X_theta[, "effort"]), expected_effort)
  expect_equal(unname(fit$X_theta[, "collection_typewet"]),
               as.numeric(samples$collection_type == "wet"))
  expect_equal(unname(fit$infos$M), c(2, 2, 2))
  expect_equal(fit$infos$list_idx$idx_z_w, c(1, 1, 2, 2, 3, 3))
}

test_that("two-stage covariates follow sample identities, not pasted-key order", {
  for (ids in c("global", "local", "character")) {
    samples <- alignment_samples(ids)
    for (shuffled in c(FALSE, TRUE)) {
      fit <- fit_alignment_data(alignment_data(samples, shuffled = shuffled))
      expect_sample_covariates(fit, samples)
      expect_equal(fit$infos$model, "two_stage")
      expect_equal(fit$infos$K, c(2, 1, 3, 1, 2, 2, 1, 2, 1))
      expected_sample <- rep(c(1, 1, 2, 3, 3, 4, 5, 6, 6),
                              c(2, 1, 3, 1, 2, 2, 1, 2, 1))
      expect_equal(fit$infos$list_idx$idx_w_k, expected_sample)
      # Tie each design row to the actual observations consumed by w.
      expect_equal(fit$infos$data_info$effort, samples$effort[expected_sample])
      expect_equal(unname(fit$infos$OTU[, "sp1"]),
                   c(5, 0, 8, 2, 0, 3)[expected_sample])
    }
  }
})

test_that("occupancy covariates follow numeric site and sample order", {
  samples <- alignment_samples()
  for (shuffled in c(FALSE, TRUE)) {
    fit <- fit_alignment_data(alignment_data(samples, "occupancy", shuffled))
    expect_equal(fit$infos$model, "occupancy")
    expect_sample_covariates(fit, samples)
    expect_equal(fit$infos$data_info$effort, samples$effort)
  }
})

test_that("intercept-only collection models keep one row per sample", {
  samples <- alignment_samples("character")
  fit <- fit_alignment_data(alignment_data(samples), covariates = character())
  expect_equal(unname(fit$X_theta), matrix(1, 6, 1))
  expect_equal(fit$infos$list_idx$idx_z_w, c(1, 1, 2, 2, 3, 3))
})

test_that("sample identifiers remain available as collection covariates", {
  samples <- alignment_samples("local")
  fit <- fit_alignment_data(alignment_data(samples, shuffled = TRUE),
                             covariates = c("Site", "Sample"))
  expect_equal(unname(fit$X_theta[, "Site"]),
               (samples$Site - mean(samples$Site)) / sd(samples$Site))
  expect_equal(unname(fit$X_theta[, "Sample"]),
               (samples$Sample - mean(samples$Sample)) / sd(samples$Sample))
})
