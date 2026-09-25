# Reanalyse completed Lesson 4 fits. This script contains no model fitting.
# Rscript export-calibration.R RUN_ROOT REPOSITORY_ROOT OUTPUT_ROOT [JOB_NUMBERS]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
root <- normalizePath(args[1]); repo <- normalizePath(args[2]); out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
code <- file.path(repo, "dev/simstudy/jsdm-package-comparison")
source(file.path(code, "pilot-math.R"))
source(file.path(code, "extension/study.R"))
source(file.path(code, "extension/math.R"))
source(file.path(code, "extension/calibration-math.R"))
manifest <- read.csv(file.path(root, "status.csv"))
stopifnot(nrow(manifest) == 640L, all(manifest$completed))
selected <- if (length(args) >= 4L) as.integer(strsplit(args[4], ",", fixed = TRUE)[[1]]) else manifest$job_number
jobs <- manifest[manifest$scored & manifest$job_number %in% selected, , drop = FALSE]
grid <- fixed_grid()
rownames(grid) <- c("mean environment", "gradient 1 low", "gradient 1 high", "gradient 2 low", "gradient 2 high")
dir.create(file.path(out, "jobs"), showWarnings = FALSE)
source_files <- c(file.path(code, "pilot-math.R"), file.path(code, "extension",
  c("study.R", "math.R", "calibration-math.R", "export-calibration.R")))
source_hashes <- unname(tools::md5sum(source_files))
started <- Sys.time()
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i, ]; path <- file.path(root, "jobs", j$job)
  inputs <- c(file.path(path, "result.rds"), file.path(path, "score.rds"),
    file.path(root, "inputs", paste0(j$input, ".rds")),
    file.path(root, "truth", paste0(j$dataset, ".rds")))
  hashes <- tools::md5sum(inputs)
  cache <- file.path(out, "jobs", paste0(j$job, ".rds"))
  if (file.exists(cache)) {
    old <- readRDS(cache)
    stopifnot(identical(old$input_hashes, hashes), identical(old$source_hashes, source_hashes))
    next
  }
  r <- readRDS(inputs[1]); score <- readRDS(inputs[2]); input <- readRDS(inputs[3]); truth <- readRDS(inputs[4])
  stopifnot(isTRUE(r$ok), !r$truth_used, score$result_md5 == unname(hashes[1]),
    score$input_md5 == j$input_md5, unname(hashes[3]) == j$input_md5,
    identical(dim(score$estimate), dim(score$truth)))
  errors <- score$estimate - score$truth
  point <- cbind(j, elements = length(errors), bias = mean(errors), mae = mean(abs(errors)),
    rmse = sqrt(mean(errors^2)))
  stopifnot(abs(point$bias * 100 - score$overall$bias_pp) < 1e-10,
            abs(point$mae * 100 - score$overall$mae_pp) < 1e-10)
  # The five raw predictor locations are fixed for every community and package.
  # Integrate hidden conditions; intervals describe the marginal probability,
  # not a realised latent site state or a future zero/one observation.
  raw <- basis(grid, j$response)
  x <- as.data.frame(sweep(sweep(as.matrix(raw), 2, input$centre, "-"), 2, input$spread, "/"))
  true_probability <- marginal_truth(truth, grid)[, seq_len(j$n_species), drop = FALSE]
  bayesian <- j$package %in% c("occJSDM", "Hmsc")
  if (bayesian) {
    probability <- probability_intervals(r$parameters, x)
    integration <- attr(probability, "integration")
  } else {
    probability <- data.frame(estimate = as.vector(point_marginal(r$parameters, x)), lower = NA_real_, upper = NA_real_)
    integration <- list(method = "adaptive normal integration; no uncertainty saved", nodes = NA_integer_, error = NA_real_)
  }
  label <- function(d, target, term, element) {
    cbind(j[rep(1, nrow(d)), , drop = FALSE], target = target, term = term, element = element, d)
  }
  probability$truth <- as.vector(true_probability)
  grid_rows <- label(probability, "Marginal probability", "Five fixed environments",
    paste(rep(colnames(input$y), each = 5L), rep(rownames(grid), j$n_species), sep = ": "))
  coefficients <- NULL
  # Probit coefficients do not have the simulation's logit coefficient truth.
  # A deliberately misspecified straight response has no matching coefficient target.
  if (j$package != "Hmsc" && !(j$scenario == "curved" && j$response == "linear")) {
    b <- raw_coefficients(r$parameters$beta, input$centre, input$spread)
    true_beta <- truth$beta[, seq_len(j$n_species), drop = FALSE]
    if (j$response == "quadratic") true_beta <- rbind(true_beta, truth$curvature[seq_len(j$n_species)])
    terms <- c("Intercept", colnames(x))
    if (bayesian) {
      d <- data.frame(estimate = as.vector(apply(b, 1:2, mean)),
        lower = as.vector(apply(b, 1:2, quantile, probs = .025)),
        upper = as.vector(apply(b, 1:2, quantile, probs = .975)))
    } else d <- data.frame(estimate = as.vector(b), lower = NA_real_, upper = NA_real_)
    d$truth <- as.vector(true_beta)
    coefficients <- label(d, "Environmental coefficient", rep(terms, j$n_species),
      paste(rep(colnames(input$y), each = length(terms)), rep(terms, j$n_species), sep = ": "))
  }
  traits <- NULL
  if (!is.null(score$traits) && j$package != "Hmsc") {
    d <- score$traits[, c("estimate", "lower", "upper", "truth")]
    terms <- paste(score$traits$trait, score$traits$environment, sep = " / ")
    traits <- label(d, "Trait coefficient", terms, terms)
  }
  result <- list(prediction = point, rows = rbind(grid_rows, coefficients, traits),
    integration = integration, job = j, input_hashes = hashes, source_hashes = source_hashes,
    saved_draws_per_chain = if (bayesian) dim(r$parameters$beta)[3] else NA_integer_,
    saved_chains = if (bayesian) dim(r$parameters$beta)[4] else NA_integer_)
  saveRDS(result, paste0(cache, ".partial")); stopifnot(file.rename(paste0(cache, ".partial"), cache))
  cat(sprintf("%d/%d %s; %.1f seconds elapsed\n", i, nrow(jobs), j$job, as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
parts <- lapply(file.path(out, "jobs", paste0(jobs$job, ".rds")), readRDS)
rows <- do.call(rbind, lapply(parts, `[[`, "rows"))
point <- do.call(rbind, lapply(parts, `[[`, "prediction"))
summary <- calibration_group_summary(rows)
per_community <- do.call(rbind, lapply(split(rows, interaction(rows$job, rows$target, rows$term, drop = TRUE)), function(d) {
  cbind(d[1, c("package", "scenario", "n_sites", "n_species", "response", "use_traits", "target", "term", "replicate", "diagnostic_pass")], calibration_summary(d))
}))
integration <- do.call(rbind, lapply(parts, function(z) data.frame(job = z$job$job, package = z$job$package,
  method = z$integration$method, nodes = z$integration$nodes, endpoint_change = z$integration$error,
  draws_per_chain = z$saved_draws_per_chain, chains = z$saved_chains)))
bundle <- list(manifest = manifest, prediction = point, rows = rows, summary = summary,
  per_community = per_community, grid = grid, integration = integration,
  complete = nrow(jobs) == sum(manifest$scored), run_root = root, exported_at = Sys.time(),
  provenance = list(sources = tools::md5sum(source_files), archive_status = tools::md5sum(file.path(root, "status.csv")),
                    session = sessionInfo(), input_hashes = lapply(parts, `[[`, "input_hashes")))
saveRDS(bundle, file.path(out, "calibration.rds"), compress = "xz")
write.csv(summary, file.path(out, "interval-summary.csv"), row.names = FALSE)
write.csv(point, file.path(out, "prediction-by-community.csv"), row.names = FALSE)
write.csv(integration, file.path(out, "integration-checks.csv"), row.names = FALSE)
cat("Exported", nrow(jobs), "saved fits; complete =", bundle$complete, "; no models fitted.\n")
cat("Run verify-calibration.R before copying this bundle into teaching-data.\n")
