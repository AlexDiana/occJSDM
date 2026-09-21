# Extend the existing, verified lesson fits. This script never fits a model.
# From the repository root:
# Rscript dev/simstudy/vignette-lesson/summarise_outputs.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)

library(occJSDM)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")

lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
input <- readRDS(file.path(archive, "input.rds"))
stopifnot(identical(lesson$input, input))
stopifnot(identical(input$source_hashes, lesson_source_hashes()))

truth <- input$sim$true_params$jsdmParams_true
species <- colnames(input$sim$data_list$OTU)
traits <- colnames(input$sim$data_list$traits)
environments <- grep("^X_psi", names(input$sim$data_list$info), value = TRUE)
trait_sd <- apply(input$sim$data_list$traits, 2, sd)
trait_matrix <- scale(input$sim$data_list$traits)

# The simulator already standardizes X before applying B and G, whereas its
# traits remain on their original scale. Only the trait scale changes G here.
trait_truth <- sweep(truth$G, 1, trait_sd, "*")

summarise_draws <- function(draws) {
  interval <- quantile(draws, c(.025, .975), names = FALSE)
  tibble(estimate = mean(draws), lower = interval[1], upper = interval[2])
}

coefficient_table <- function(draws, values, block) {
  expand_grid(covariate = dimnames(draws)[[2]], term = dimnames(draws)[[3]]) |>
    mutate(
      block = block,
      truth = map2_dbl(covariate, term, ~ values[.x, .y]),
      summary = map2(covariate, term, ~ summarise_draws(draws[, .x, .y]))
    ) |>
    unnest(summary) |>
    mutate(excludes_zero = lower > 0 | upper < 0)
}

dimnames(truth$B) <- list(environments, species)
dimnames(trait_truth) <- list(traits, environments)

# An oracle diagnostic: regress the KNOWN species slopes on the measured
# traits. By linearity, the three contributions add to the realized slope.
# This is not a replacement for the generating trait coefficient G.
components <- list(
  "Generating measured-trait effect" = input$sim$data_list$traits %*% truth$G,
  "Unmeasured species trait" = truth$A %*% truth$C,
  "Remaining species differences" = truth$Bt
)
trait_components <- imap_dfr(components, function(component, label) {
  slopes <- coef(lm(component ~ trait_matrix))[-1, , drop = FALSE]
  expand_grid(trait = traits, covariate = environments) |>
    mutate(component = label,
           value = map2_dbl(trait, covariate,
                            ~ slopes[match(.x, traits), match(.y, environments)]))
})
realized_slopes <- coef(lm(t(truth$B) ~ trait_matrix))[-1, , drop = FALSE]

covariance <- crossprod(truth$L)
factor_sd <- sqrt(diag(covariance))
correlation_truth <- covariance / outer(factor_sd, factor_sd)
# A constant residual has no correlation, even with itself. Keep it undefined.
correlation_truth[outer(factor_sd == 0, factor_sd == 0, "|")] <- NA_real_
dimnames(correlation_truth) <- list(species, species)

all_outputs <- list()
for (arm in c("perfect", "default")) {
  manifest <- lesson$manifests[[arm]]
  fit_path <- file.path(archive, manifest$file)
  stopifnot(identical(unname(tools::md5sum(fit_path)), manifest$md5))
  saved <- readRDS(fit_path)
  fit <- saved$fit
  validate_lesson_fit_identity(fit, input)
  stopifnot(identical(saved$source_hashes, input$source_hashes), fit$infos$ps == 0)
  stopifnot(identical(fit$infos$speciesNames, species))
  stopifnot(max(abs(fit$Tr - trait_matrix)) < 1e-12)
  # Verifies the simulator/fitter environmental scale and the complete truth.
  reconstructed <- sweep(fit$X_psi %*% truth$B + truth$U %*% truth$L,
                         2, truth$B0, "+")
  stopifnot(max(abs(reconstructed - truth$eta)) < 1e-12)

  coefficients <- bind_rows(
    coefficient_table(returnOccupancyCovariates(fit), truth$B, "Environment"),
    coefficient_table(returnTraitsCoeff(fit), t(trait_truth), "Trait")
  )
  j <- fit$results_output$jsdm_output
  coefficient_diagnostics <- pmap_dfr(
    coefficients[c("block", "covariate", "term")],
    function(block, covariate, term) {
      e <- match(covariate, environments)
      draws <- if (block == "Environment") {
        j$B_output[e, match(term, species), , ]
      } else {
        j$G_output[match(term, traits), e, , ]
      }
      tibble(rhat = posterior::rhat(draws),
             ess_bulk = posterior::ess_bulk(draws),
             ess_tail = posterior::ess_tail(draws))
    }
  )
  coefficients <- bind_cols(coefficients, coefficient_diagnostics)

  gradients <- map_dfr(environments, function(covariate) {
    grid <- returnOccupancyGradient(fit, covName = covariate, n_grid = 40L)
    other <- setdiff(environments, covariate)
    other_median <- apply(fit$X_psi[, other, drop = FALSE], 2, median)
    # Calculate outside a data mask: the table also has a species column.
    generating_probability <- map2_dbl(grid$species, grid$x, function(sp, value) {
      s <- match(sp, species)
      plogis(truth$B0[s] + value * truth$B[covariate, s] +
               sum(other_median * truth$B[other, s]))
    })
    mutate(grid, covariate = covariate, truth = generating_probability)
  })

  correlation <- returnResidualCorrelationMatrix(fit)
  correlations <- expand_grid(species1 = species, species2 = species) |>
    mutate(
      truth = map2_dbl(species1, species2, ~ correlation_truth[.x, .y]),
      estimate = map2_dbl(species1, species2, ~ correlation[2, .x, .y]),
      lower = map2_dbl(species1, species2, ~ correlation[1, .x, .y]),
      upper = map2_dbl(species1, species2, ~ correlation[3, .x, .y])
    )

  variation <- returnVariancePartitioning(fit) |>
    as_tibble() |>
    select(species = Species, Environmental = Env, Spatial, Residual = Biotic) |>
    pivot_longer(-species, names_to = "component", values_to = "estimate")
  variation_truth <- as_tibble(truth$varPart) |>
    mutate(species = species) |>
    select(species, Environmental, Spatial, Residual = Biotic) |>
    pivot_longer(-species, names_to = "component", values_to = "truth")
  variation <- left_join(variation, variation_truth,
                         by = c("species", "component"), relationship = "one-to-one")

  j <- fit$results_output$jsdm_output
  ni <- dim(j$B0_output)[2]
  nc <- dim(j$B0_output)[3]
  residual_mean <- matrix(0, nrow(fit$X_psi), length(species))
  for (chain in seq_len(nc)) for (iteration in seq_len(ni)) {
    residual_mean <- residual_mean +
      j$U_output[, , iteration, chain] %*% j$L_output[, , iteration, chain]
  }
  residual_mean <- residual_mean / (ni * nc)
  residual <- expand_grid(species = species, Site = as.character(fit$infos$siteNames)) |>
    mutate(truth = c(truth$U %*% truth$L), estimate = c(residual_mean))

  baseline <- map_dfr(seq_along(species), function(s) {
    summarise_draws(plogis(j$B0_output[s, , ])) |>
      mutate(species = species[s], truth = plogis(truth$B0[s]))
  })

  detection_effort <- collection <- NULL
  if (arm == "default") {
    collection_truth <- input$sim$true_params$beta_theta_true
    scale_info <- fit$infos$list_X_theta_mat
    collection_truth[1, ] <- collection_truth[1, ] +
      scale_info$mean_df * collection_truth[2, ]
    collection_truth[2, ] <- scale_info$sd_df * collection_truth[2, ]
    dimnames(collection_truth) <- list(colnames(fit$X_theta), species)
    collection <- coefficient_table(returnCollectionCovariates(fit),
                                     collection_truth, "Collection")
    theta <- apply(fit$results_output$beta_theta_output[1, , , , drop = FALSE], 2, c) |>
      plogis()
    p <- apply(fit$results_output$p_output, c(1, 2), c)
    p_truth <- effective_detection_rate(input$params$p, input$params$mu1,
                                         input$params$sigma1)
    theta_truth <- plogis(collection_truth[1, ])
    detection_effort <- expand_grid(M = 1:4, K = 1:6) |>
      mutate(summary = map2(M, K, function(M, K) {
        missed_pcr <- apply((1 - p)^K, c(1, 3), prod)
        expected <- rowSums(1 - (1 - theta * (1 - missed_pcr))^M)
        missed_truth <- apply((1 - p_truth)^K, 2, prod)
        summarise_draws(expected) |>
          mutate(truth = sum(1 - (1 - theta_truth * (1 - missed_truth))^M))
      })) |>
      unnest(summary)
  }

  all_outputs[[arm]] <- list(coefficients = coefficients, gradients = gradients,
    correlations = correlations, variation = variation, residual = residual,
    baseline = baseline, collection = collection, detection_effort = detection_effort)
  rm(fit, saved, j)
  invisible(gc())
}

combine <- function(field) imap_dfr(all_outputs, ~ mutate(.x[[field]], arm = .y))
legacy <- new.env()
load("data/sampleresults.rda", envir = legacy)
legacy_counts <- imap_dfr(list(
  environment = returnOccupancyCovariates(legacy$sampleresults),
  trait = returnTraitsCoeff(legacy$sampleresults)
), function(draws, block) {
  interval <- apply(draws, 2:3, quantile, c(.025, .975))
  tibble(block = block, total = length(interval[1, , ]),
         excludes_zero = sum(interval[1, , ] * interval[2, , ] > 0))
})

outputs <- list(
  schema = 1L,
  coefficients = combine("coefficients"), gradients = combine("gradients"),
  correlations = combine("correlations"), variation = combine("variation"),
  residual = combine("residual"), baseline = combine("baseline"),
  collection = all_outputs$default$collection,
  detection_effort = all_outputs$default$detection_effort,
  trait_components = trait_components,
  realized_trait_slopes = realized_slopes,
  trait_hidden_correlation = cor(input$sim$data_list$traits, truth$A),
  legacy_counts = legacy_counts,
  legacy_fit_md5 = unname(tools::md5sum("data/sampleresults.rda")),
  source_hashes = input$source_hashes,
  lesson_md5 = unname(tools::md5sum("vignettes/teaching-data/nonspatial-lesson.rds")),
  fit_manifests = lesson$manifests[c("perfect", "default")],
  exporter_md5 = unname(tools::md5sum("dev/simstudy/vignette-lesson/summarise_outputs.R"))
)
saveRDS(outputs, "vignettes/teaching-data/output-lesson.rds", compress = "xz")
print(outputs$coefficients |> group_by(arm, block) |>
        summarise(total = n(), intervals_excluding_zero = sum(excludes_zero), .groups = "drop"))
print(trait_components |> filter(covariate == environments[1], trait == traits[1]))
