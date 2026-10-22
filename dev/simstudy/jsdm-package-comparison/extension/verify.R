args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1]); repo <- normalizePath(args[2])
bundle <- readRDS(file.path(repo, "vignettes/teaching-data/lesson-4-extension.rds"))
manifest <- read.csv(file.path(root, "manifest.csv"))
stopifnot(nrow(manifest) == 640L, !anyDuplicated(manifest$job),
          identical(bundle$manifest$job, manifest$job))
for (i in seq_len(nrow(manifest))) {
  job <- manifest[i, ]
  out <- file.path(root, "jobs", job$job)
  score_file <- file.path(out, "score.rds")
  if (!file.exists(score_file)) next
  s <- readRDS(score_file)
  input <- readRDS(file.path(root, "inputs", paste0(job$input, ".rds")))
  truth <- readRDS(file.path(root, "truth", paste0(job$dataset, ".rds")))
  stopifnot(unname(tools::md5sum(file.path(out, "result.rds"))) == s$result_md5,
            identical(input$y, truth$y[seq_len(job$n_sites), seq_len(job$n_species), drop = FALSE]),
            identical(s$observations, truth$y[301:600, seq_len(job$n_species), drop = FALSE]))
  # Recompute independently from prediction cells, rather than reusing the scorer.
  signed <- (as.vector(s$estimate) - as.vector(s$truth)) * 100
  stopifnot(abs(sum(signed) / length(signed) - s$overall$bias_pp) < 1e-10,
            abs(sum(abs(signed)) / length(signed) - s$overall$mae_pp) < 1e-10,
            sum(s$bands$n) == length(signed))
  # Independently check one true marginal against a direct normal integral.
  x <- truth$raw_x[301, ]
  mu <- sum(c(1, x) * truth$beta[, 1]) + truth$curvature[1] * x[1]^2
  sd <- sqrt(sum(truth$loading[, 1]^2))
  expected <- integrate(function(z) plogis(mu + sd * z) * dnorm(z), -Inf, Inf,
                         rel.tol = 1e-10)$value
  stopifnot(abs(expected - s$truth[1, 1]) < 1e-8)
  row <- bundle$overall[bundle$overall$job_number == job$job_number, ]
  stopifnot(nrow(row) == 1, abs(row$mae_pp - s$overall$mae_pp) < 1e-10)
}
cat("Verified input identity, held-out states, independent probability truth and error summaries for",
    bundle$scored, "scored results;", bundle$completed, "of 640 fits attempted.\n")
