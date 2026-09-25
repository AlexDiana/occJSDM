args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1]); repo <- normalizePath(args[2])
code <- file.path(repo, "dev/simstudy/jsdm-package-comparison")
source(file.path(code, "extension/study.R"))
source(file.path(code, "pilot-math.R"))
source(file.path(code, "extension/math.R"))
manifest <- read.csv(file.path(root, "manifest.csv"))
status <- manifest
status$completed <- FALSE; status$fit_ok <- FALSE; status$diagnostic_pass <- FALSE
status$scored <- FALSE; status$error <- NA_character_
overall <- species <- bands <- traits <- list()
for (i in seq_len(nrow(manifest))) {
  job <- manifest[i, ]; out <- file.path(root, "jobs", job$job)
  if (!file.exists(file.path(out, "status.json"))) next
  s <- jsonlite::read_json(file.path(out, "status.json"), simplifyVector = TRUE)
  stopifnot(s$input_md5 == job$input_md5,
    s$result_md5 == unname(tools::md5sum(file.path(out, "result.rds"))))
  status$completed[i] <- TRUE; status$fit_ok[i] <- isTRUE(s$ok)
  status$diagnostic_pass[i] <- isTRUE(s$diagnostic_pass)
  if (!is.null(s$error)) status$error[i] <- s$error
  score_path <- file.path(out, "score.rds")
  if (!file.exists(score_path)) next
  score <- readRDS(score_path)
  stopifnot(score$result_md5 == s$result_md5,
    identical(dim(score$truth), dim(score$estimate)),
    abs(score$overall$mae_pp - 100 * mean(abs(score$estimate - score$truth))) < 1e-10,
    abs(score$overall$bias_pp - 100 * mean(score$estimate - score$truth)) < 1e-10,
    abs(score$overall$brier - mean((score$estimate - score$observations)^2)) < 1e-10)
  status$scored[i] <- TRUE
  label <- cbind(job, diagnostic_pass = isTRUE(s$diagnostic_pass))
  overall[[length(overall) + 1L]] <- cbind(label, score$overall)
  species[[length(species) + 1L]] <- cbind(label[rep(1, nrow(score$species)), ], score$species)
  bands[[length(bands) + 1L]] <- cbind(label[rep(1, nrow(score$bands)), ], score$bands)
  if (!is.null(score$traits)) traits[[length(traits) + 1L]] <- cbind(label[rep(1, nrow(score$traits)), ], score$traits)
}
bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
# Generating curves illustrate the actual prespecified first community.
illustrations <- list()
for (scenario in c("baseline", "rare", "correlated", "curved")) {
  community <- readRDS(file.path(root, "truth", sprintf("%s-r01.rds", scenario)))
  x <- data.frame(environment_1 = seq(-2, 2, length.out = 81), environment_2 = 0)
  truth <- marginal_truth(community, x)
  illustrations[[scenario]] <- data.frame(scenario, environment = x$environment_1,
    probability = truth[, 1], species = "species_01")
}
correlated <- readRDS(file.path(root, "truth/correlated-r01.rds"))
curved <- readRDS(file.path(root, "truth/curved-r01.rds"))
# A prespecified species from community 1 illustrates fitted response shapes.
fitted_curves <- list()
curve_jobs <- status[status$scored & status$scenario == "curved" & status$replicate == 1, ]
for (i in seq_len(nrow(curve_jobs))) {
  job <- curve_jobs[i, ]
  out <- file.path(root, "jobs", job$job)
  r <- readRDS(file.path(out, "result.rds")); score <- readRDS(file.path(out, "score.rds"))
  input <- readRDS(file.path(root, "inputs", paste0(job$input, ".rds")))
  community <- readRDS(file.path(root, "truth", paste0(job$dataset, ".rds")))
  raw <- data.frame(environment_1 = seq(-2, 2, length.out = 81), environment_2 = 0)
  x <- as.data.frame(sweep(sweep(as.matrix(basis(raw, job$response)), 2, input$centre, "-"), 2, input$spread, "/"))
  truth <- marginal_truth(community, raw)[, 1]
  p <- r$parameters
  if (job$package %in% c("occJSDM", "Hmsc")) {
    p$beta <- p$beta[, 1, , , drop = FALSE]; p$loading <- p$loading[, 1, , , drop = FALSE]
    idx <- unique(round(seq(1, dim(p$beta)[3], length.out = score$numeric_check$draws_per_chain)))
    estimate <- as.vector(posterior_marginal(p, x, nodes = score$numeric_check$nodes, draw_indices = idx))
  } else {
    p$beta <- p$beta[, 1, drop = FALSE]; p$loading <- p$loading[, 1, drop = FALSE]
    estimate <- as.vector(point_marginal(p, x))
  }
  fitted_curves[[i]] <- data.frame(package = job$package, n_sites = job$n_sites,
    response = job$response, diagnostic_pass = job$diagnostic_pass,
    environment = raw$environment_1, truth, estimate)
}
bundle <- list(manifest = status, overall = bind(overall), species = bind(species),
  bands = bind(bands), traits = bind(traits), generating_curves = bind(illustrations),
  fitted_curves = bind(fitted_curves),
  curve_training_raw = as.data.frame(curved$raw_x[1:100, ]),
  curve_test_raw = as.data.frame(curved$raw_x[301:600, ]),
  correlated_environment = as.data.frame(correlated$raw_x[1:300, ]),
  exported_at = Sys.time(), planned_jobs = nrow(status), completed = sum(status$completed),
  scored = sum(status$scored), run_root = root,
  all_attempted = all(status$completed),
  sources = tools::md5sum(list.files(file.path(code, "extension"), full.names = TRUE)))
path <- file.path(repo, "vignettes/teaching-data/lesson-4-extension.rds")
saveRDS(bundle, path, compress = "xz")
write.csv(status, file.path(root, "status.csv"), row.names = FALSE)
write.csv(bundle$overall, file.path(root, "overall-errors.csv"), row.names = FALSE)
write.csv(bundle$traits, file.path(root, "trait-effects.csv"), row.names = FALSE)
cat("Exported", bundle$completed, "completed fits and", bundle$scored,
    "scored results from", bundle$planned_jobs, "planned combinations.\n")
