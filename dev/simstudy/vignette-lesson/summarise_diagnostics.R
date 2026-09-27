# Preserve selected complete chain traces for teaching; never run MCMC here.
# Rscript dev/simstudy/vignette-lesson/summarise_diagnostics.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")

lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
input <- readRDS(file.path(archive, "input.rds"))
stopifnot(identical(input, lesson$input),
          identical(input$source_hashes, lesson_source_hashes()))

read_checked_fit <- function(arm) {
  manifest <- lesson$manifests[[arm]]
  path <- file.path(archive, manifest$file)
  stopifnot(identical(unname(tools::md5sum(path)), manifest$md5))
  saved <- readRDS(path)
  stopifnot(identical(saved$source_hashes, input$source_hashes),
            identical(saved$input_md5, lesson$input_md5),
            identical(saved$mcmc, manifest$mcmc))
  validate_lesson_fit_identity(saved$fit, input)
  # The saved public diagnostics are used in the visible filtering example.
  stopifnot(isTRUE(all.equal(
    as.data.frame(occJSDM::returnConvergenceDiagnostics(saved$fit)),
    as.data.frame(lesson$diagnostics[[arm]]), check.attributes = FALSE
  )))
  saved$fit
}

fit <- read_checked_fit("default")
species <- c("OTU_1", "OTU_6")
species_index <- match(species, fit$infos$speciesNames)
covariate_index <- match("X_theta", colnames(fit$X_theta))
primer_index <- match("1", as.character(fit$infos$primerNames))
stopifnot(!anyNA(c(species_index, covariate_index, primer_index)))

collection_truth <- input$sim$true_params$beta_theta_true[2, species_index] *
  fit$infos$list_X_theta_mat$sd_df
p_truth <- effective_detection_rate(input$params$p, input$params$mu1,
                                     input$params$sigma1)
traces <- list(
  collection = list(
    arm = "default", param = "beta_theta", label1 = "X_theta", label2 = species,
    draws = fit$results_output$beta_theta_output[covariate_index, species_index, , , drop = FALSE],
    truth = data.frame(label1 = "X_theta", label2 = species, truth = collection_truth)
  ),
  detection = list(
    arm = "default", param = "p", label1 = "1", label2 = species,
    draws = fit$results_output$p_output[primer_index, species_index, , , drop = FALSE],
    truth = data.frame(label1 = "1", label2 = species, truth = p_truth[primer_index, species_index])
  )
)

# Choose the field-contamination parameter with the lowest public ESS in the
# existing longer alternative-prior fit, not the largest error against truth.
diagnostics <- lesson$diagnostics$alternative
candidates <- diagnostics[diagnostics$param == "theta0" & is.finite(diagnostics$ess), ]
candidates <- candidates[order(candidates$ess, candidates$label1), ]
flagged_species <- candidates$label1[1]
fit <- read_checked_fit("alternative")
s <- match(flagged_species, fit$infos$speciesNames)
stopifnot(!is.na(s))
traces$field_contamination <- list(
  arm = "alternative", param = "theta0", label1 = flagged_species, label2 = "1",
  draws = fit$results_output$theta0_output[s, , , drop = FALSE],
  truth = data.frame(label1 = flagged_species, label2 = "1", truth = input$params$theta0[s])
)

bundle <- list(
  schema = 1L, traces = traces,
  selection = "OTU_1 and OTU_6 for collection/PCR; lowest theta0 ESS in the saved longer alternative-prior fit",
  fit_manifests = lesson$manifests[c("default", "alternative")],
  source_hashes = input$source_hashes,
  lesson_md5 = unname(tools::md5sum("vignettes/teaching-data/nonspatial-lesson.rds")),
  exporter_md5 = unname(tools::md5sum("dev/simstudy/vignette-lesson/summarise_diagnostics.R"))
)
saveRDS(bundle, "vignettes/teaching-data/diagnostics-lesson.rds", compress = "xz")
print(candidates[1, c("param", "label1", "rhat", "ess")])
cat("Exported complete, unthinned post-burn-in traces; no model fitting.\n")
