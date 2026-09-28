# Separate scorer: fit/selection is complete before this script opens truth.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1]); code <- normalizePath(args[2]); number <- as.integer(args[3])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "extension", "study.R"))
source(file.path(code, "extension", "math.R"))
manifest <- read.csv(file.path(root, "manifest.csv"))
job <- manifest[manifest$job_number == number, ]
out <- file.path(root, "jobs", job$job)
result_path <- file.path(out, "result.rds")
if (file.exists(file.path(out, "score.rds"))) quit(status = 0)
r <- readRDS(result_path)
stopifnot(!r$truth_used, r$input_md5 == job$input_md5)
if (!isTRUE(r$ok)) quit(status = 0)
input <- readRDS(file.path(root, "inputs", paste0(job$input, ".rds")))
if (job$package %in% c("occJSDM", "Hmsc")) {
  grid <- input$x[c(1, 21, 41, 61, 81), , drop = FALSE]
  reference_grid <- matrix(r$diagnostics$table$mean[r$diagnostics$table$block == "new_site_probability"], nrow = 5)
  numeric_check <- checked_posterior_mean(r$parameters, input$test_x, grid, reference_grid)
  estimate <- numeric_check$mean
  numeric_check$mean <- NULL
} else {
  estimate <- point_marginal(r$parameters, input$test_x)
  numeric_check <- list(method = "independent adaptive logit-normal integration")
}
stopifnot(identical(dim(estimate), c(300L, as.integer(job$n_species))),
          all(is.finite(estimate)), all(estimate >= 0 & estimate <= 1))
# This is the first read of generating truth in this job's workflow.
community <- readRDS(file.path(root, "truth", paste0(job$dataset, ".rds")))
truth_path <- file.path(root, "truth", paste0(job$dataset, "-marginal.rds"))
if (file.exists(truth_path)) truth_all <- readRDS(truth_path) else {
  truth_all <- marginal_truth(community, community$raw_x[301:600, ])
  # Per-job scoring can run independently; atomic rename prevents partial reads.
  temporary <- paste0(truth_path, ".", Sys.getpid())
  saveRDS(truth_all, temporary); file.rename(temporary, truth_path)
}
truth <- truth_all[, seq_len(job$n_species), drop = FALSE]
observations <- community$y[301:600, seq_len(job$n_species), drop = FALSE]
overall <- error_summary(truth, estimate)
overall$brier <- mean((estimate - observations)^2)
overall$negative_log_score <- -mean(dbinom(observations, 1, pmin(1 - 1e-12, pmax(1e-12, estimate)), log = TRUE))
species <- do.call(rbind, lapply(seq_len(job$n_species), function(s) {
  cbind(data.frame(species = colnames(input$y)[s], training_presences = sum(input$y[, s]),
    mean_true_probability = mean(truth[, s])), error_summary(truth[, s], estimate[, s]))
}))
band <- cut(truth, c(-Inf, .2, .8, Inf), right = FALSE, labels = c("below 20%", "20% to 80%", "80% or higher"))
bands <- do.call(rbind, lapply(levels(band), function(label) {
  keep <- band == label
  if (!any(keep)) return(data.frame(band = label, n = 0, bias_pp = NA_real_, mae_pp = NA_real_))
  cbind(data.frame(band = label), error_summary(truth[keep], estimate[keep]))
}))
trait_rows <- NULL
if (isTRUE(job$use_traits)) {
  trait_names <- names(input$traits)
  trait_spread <- apply(input$traits, 2, sd)
  for (t in 1:2) for (e in 1:2) {
    denominator <- trait_spread[t] * input$spread[e]
    if (job$package == "occJSDM") {
      draws <- as.vector(r$parameters$trait[t, e, , ]) / denominator
      estimate_g <- mean(draws); interval <- quantile(draws, c(.025, .975))
    } else if (job$package == "gllvm") {
      native <- readRDS(file.path(out, sprintf("start-%d.rds", r$selected_start)))
      term <- paste(names(input$x)[e], trait_names[t], sep = ":")
      estimate_g <- native$native_parameters$B[term] / denominator
      se <- native$standard_errors$B[term] / denominator
      interval <- estimate_g + c(-1, 1) * qnorm(.975) * se
    } else {
      # Hmsc's Gamma is environmental coefficient by trait, including intercepts.
      # Its probit coefficient has no exact logit generating counterpart.
      draws <- as.vector(r$parameters$trait[e + 1, t + 1, , ]) / input$spread[e]
      estimate_g <- mean(draws); interval <- quantile(draws, c(.025, .975))
    }
    trait_rows <- rbind(trait_rows, data.frame(trait = trait_names[t], environment = names(input$x)[e],
      estimate = unname(estimate_g), lower = unname(interval[1]), upper = unname(interval[2]),
      truth = if (job$package == "Hmsc") NA_real_ else community$trait_effect[t, e],
      true_direction = sign(community$trait_effect[t, e]), link = if (job$package == "Hmsc") "probit" else "logit",
      scale = if (job$package == "Hmsc") "raw trait and environmental units, probit; no exact logit coefficient truth" else "raw trait and environmental units"))
  }
}
score <- list(job = job, diagnostic_pass = r$diagnostic_pass, overall = overall,
  species = species, bands = bands, traits = trait_rows, estimate = estimate,
  truth = truth, observations = observations, numeric_check = numeric_check,
  result_md5 = unname(tools::md5sum(result_path)), input_md5 = job$input_md5,
  scored_at = Sys.time())
saveRDS(score, file.path(out, "score.partial.rds"))
stopifnot(file.rename(file.path(out, "score.partial.rds"), file.path(out, "score.rds")))
print(cbind(job[, c("job_number", "package", "scenario", "n_sites", "n_species")], overall))
