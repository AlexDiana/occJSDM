# Compare an sjSDM v0.1.0 rerun of one multistart start with the recorded
# v0.2.1 fit of the same start (TODO: "Check that sjSDM v0.1.0 gives the
# same fit as the recorded one"). Training data only; no truth is read.
# Verdict, stated before the rerun from check-sjsdm-multistart.R's within-
# basin criteria: "bit-identical" if both endpoints match exactly; else
# "same solution" if the continued endpoints refine to the same basin, their
# precise penalised scores differ by at most 0.1 and their fixed-grid
# predictions by at most 1 percentage point; else "different solution".
#   Rscript compare-sjsdm-v010.R ROOT CODE START
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
start <- as.integer(args[3])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
archive <- file.path(code, "stability-resolution-archive/stability-resolution")
old_prefix <- file.path(archive, "multistart", sprintf("sjSDM-multistart-start-%d", start))
new_prefix <- file.path(root, "stability-resolution/multistart", sprintf("sjSDM-v0.1.0-multistart-start-%d", start))
out_dir <- file.path(root, "comparison")
stopifnot(!file.exists(file.path(out_dir, "comparison.rds")))
dir.create(out_dir, showWarnings = FALSE)
read_run <- function(prefix) list(fit = readRDS(paste0(prefix, "-fit.rds")), final = readRDS(paste0(prefix, ".rds")),
  record = readRDS(paste0(prefix, "-record.rds")), trajectory = read.csv(paste0(prefix, "-trajectory.csv")))
old <- read_run(old_prefix)
new <- read_run(new_prefix)
stopifnot(isTRUE(old$record$ok), isTRUE(new$record$ok),
          identical(new$record$settings$data_md5, old$record$settings$data_md5),
          identical(unname(tools::md5sum(file.path(root, "inputs/training.rds"))), old$record$settings$data_md5),
          old$record$settings$backend == "torch_tmp", new$record$settings$backend == "tmp",
          new$record$settings$sjsdm_release == "v0.1.0")
common <- c("start", "fit_seed", "polish_seed", "epochs", "learning_rate", "polish_epochs",
            "polish_block", "polish_learning_rate", "sampling", "batch", "optimiser", "weight_decay")
stopifnot(identical(new$record$settings[common], old$record$settings[common]))

n_sites <- nrow(input$x); decay <- old$record$settings$weight_decay
penalty <- function(theta) n_sites * decay * sum(theta^2) / 2
precise <- make_joint_gradient(input$x, input$y, 121L)
coarse <- make_joint_gradient(input$x, input$y, 81L)
fine <- make_joint_gradient(input$x, input$y, 161L)
procrustes <- function(A, B) { s <- svd(A %*% t(B)); s$u %*% t(s$v) }
aligned_distance <- function(p, q) {
  Q <- procrustes(q$loading, p$loading)
  sqrt(sum((p$beta - q$beta)^2) + sum((Q %*% p$loading - q$loading)^2))
}
refs <- lapply(1:2, function(i) readRDS(file.path(archive, "fits",
  sprintf("penalised-reference-start-%d.rds", i)))$parameters)
names(refs) <- c("A", "B")

parameter_differences <- function(a, b) data.frame(
  identical = identical(a, b),
  max_abs_beta = max(abs(a$beta - b$beta)),
  max_abs_loading = max(abs(a$loading - b$loading)),
  max_abs_covariance = max(abs(crossprod(a$loading) - crossprod(b$loading))),
  aligned_distance = aligned_distance(a, b))
# As in check-sjsdm-multistart.R: precise native score, then a deterministic
# refinement to the basin's stationary point and its distance to the two
# penalised references.
assess <- function(p) {
  theta <- c(p$beta, p$loading)
  a <- precise(theta)
  fit <- optim(theta, function(t) coarse(t)$value + penalty(t),
               function(t) coarse(t)$gradient + n_sites * decay * t,
               method = "BFGS", control = list(maxit = 300, reltol = 1e-12))
  q <- list(beta = matrix(fit$par[1:30], 3, 10), loading = matrix(fit$par[-(1:30)], 2, 10))
  g <- fine(fit$par)
  dist <- vapply(refs, function(ref) aligned_distance(q, ref), numeric(1))
  data.frame(native_penalised_score = -a$value - penalty(theta),
             native_max_gradient = max(abs(a$gradient + n_sites * decay * theta)),
             refined_penalised_score = -g$value - penalty(fit$par),
             refined_convergence = fit$convergence,
             distance_to_A = dist[["A"]], distance_to_B = dist[["B"]],
             basin = if (min(dist) < 0.05) names(dist)[which.min(dist)] else "other")
}

stages <- rbind(
  cbind(stage = "after 3,000 epochs", parameter_differences(new$fit$parameters, old$fit$parameters),
        max_abs_history = max(abs(new$fit$history - old$fit$history))),
  cbind(stage = "after continuation", parameter_differences(new$final$parameters, old$final$parameters),
        max_abs_history = max(abs(new$final$polish_history - old$final$polish_history))))
assessment <- rbind(cbind(run = "recorded v0.2.1", assess(old$final$parameters)),
                    cbind(run = "v0.1.0 rerun", assess(new$final$parameters)))
grid_pp <- 100 * max(abs(point_marginal(new$final$parameters, fixed_grid()) -
                           point_marginal(old$final$parameters, fixed_grid())))
native_raw <- max(abs(new$final$native_raw - old$final$native_raw))
score_difference <- assessment$native_penalised_score[2] - assessment$native_penalised_score[1]
trajectory <- merge(old$trajectory, new$trajectory, by = c("start", "stage", "epoch"),
                    suffixes = c("_recorded", "_v010"))
trajectory$score_difference <- trajectory$penalised_score_v010 - trajectory$penalised_score_recorded
verdict <- if (all(stages$identical)) "bit-identical" else if (
  assessment$basin[1] == assessment$basin[2] && assessment$basin[1] != "other" &&
  abs(score_difference) <= 0.1 && grid_pp <= 1) "same solution" else "different solution"
summary <- data.frame(start = start, verdict = verdict, basin_recorded = assessment$basin[1],
  basin_v010 = assessment$basin[2], native_score_difference = score_difference,
  grid_prediction_max_difference_pp = grid_pp, native_raw_max_difference = native_raw,
  recorded_seconds = old$record$elapsed_seconds, v010_seconds = new$record$elapsed_seconds,
  recorded_python = old$record$settings$python, v010_python = new$record$settings$python,
  v010_torch = new$record$settings$torch, v010_warnings = length(new$record$warnings))
print(stages); print(assessment); print(trajectory[, c("stage", "epoch", "score_difference")]); print(summary)
write.csv(stages, file.path(out_dir, "stage-differences.csv"), row.names = FALSE)
write.csv(assessment, file.path(out_dir, "assessment.csv"), row.names = FALSE)
write.csv(trajectory, file.path(out_dir, "trajectory.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
saveRDS(list(stages = stages, assessment = assessment, trajectory = trajectory, summary = summary,
             session = sessionInfo()), file.path(out_dir, "comparison.rds"))
cat("Verdict:", verdict, "\n")
