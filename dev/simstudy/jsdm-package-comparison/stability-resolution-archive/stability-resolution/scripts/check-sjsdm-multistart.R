# Classify every native weak-penalty sjSDM endpoint by the stationary solution
# it belongs to, then apply the recorded stability assessment and selection
# rule. Training data only: no test outcomes or generating truth are read.
# Starts 1:3 are the original polished endpoints; 4:12 are the fresh starts.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
out_dir <- file.path(root, "stability-resolution/multistart")
stopifnot(!file.exists(file.path(out_dir, "multistart-summary.csv")))
n_sites <- nrow(input$x); decay <- 0.0001
penalty <- function(theta) n_sites * decay * sum(theta^2) / 2
precise <- make_joint_gradient(input$x, input$y, 121L)
coarse <- make_joint_gradient(input$x, input$y, 81L)
fine <- make_joint_gradient(input$x, input$y, 161L)
procrustes <- function(A, B) { s <- svd(A %*% t(B)); s$u %*% t(s$v) }
aligned_distance <- function(p, q) {
  Q <- procrustes(q$loading, p$loading)
  sqrt(sum((p$beta - q$beta)^2) + sum((Q %*% p$loading - q$loading)^2))
}
refs <- lapply(1:2, function(i) readRDS(file.path(root, "stability-resolution/fits",
  sprintf("penalised-reference-start-%d.rds", i)))$parameters)
names(refs) <- c("A", "B")   # A: score about -484.30 (start 1); B: about -484.43 (starts 2, 3)

endpoint <- function(start) {
  if (start <= 3L) {
    f <- readRDS(file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-polish-start-%d.rds", start)))
    r <- readRDS(file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-polish-start-%d-record.rds", start)))
  } else {
    f <- readRDS(file.path(out_dir, sprintf("sjSDM-multistart-start-%d.rds", start)))
    r <- readRDS(file.path(out_dir, sprintf("sjSDM-multistart-start-%d-record.rds", start)))
  }
  stopifnot(isTRUE(r$ok), r$settings$weight_decay == decay,
            max(abs(f$native_raw - ((cbind(1, as.matrix(input$x)) %*% f$parameters$beta) * .999999 + .0000005))) < 1e-10)
  f$parameters
}

rows <- list(); predictions <- list(); refined <- list()
for (start in 1:12) {
  p <- endpoint(start)
  theta <- c(p$beta, p$loading)
  a <- precise(theta)
  b <- joint_integration(p, input$x, input$y, 241L, clamp = TRUE)
  stopifnot(abs(a$value + b$loglik) < .001)
  native_score <- -a$value - penalty(theta)
  native_gradient <- max(abs(a$gradient + n_sites * decay * theta))
  # Deterministic refinement to the basin's stationary point, as in
  # sjsdm-penalised-reference.R, then checked at 161 nodes.
  objective <- function(t) coarse(t)$value + penalty(t)
  gradient <- function(t) coarse(t)$gradient + n_sites * decay * t
  fit <- optim(theta, objective, gradient, method = "BFGS",
               control = list(maxit = 300, reltol = 1e-12))
  q <- list(beta = matrix(fit$par[1:30], 3, 10), loading = matrix(fit$par[-(1:30)], 2, 10))
  g <- fine(fit$par)
  refined_score <- -g$value - penalty(fit$par)
  refined_gradient <- max(abs(g$gradient + n_sites * decay * fit$par))
  check_241 <- joint_integration(q, input$x, input$y, 241L, clamp = TRUE)
  dist <- vapply(refs, function(ref) aligned_distance(q, ref), numeric(1))
  basin <- if (min(dist) < 0.05) names(dist)[which.min(dist)] else "other"
  predictions[[start]] <- point_marginal(p, fixed_grid())
  refined[[start]] <- q
  rows[[start]] <- data.frame(start = start, origin = if (start <= 3L) "original" else "fresh",
    native_penalised_score = native_score, native_max_gradient = native_gradient,
    native_difference_121_241 = abs(a$value + b$loglik),
    refined_penalised_score = refined_score, refined_max_gradient = refined_gradient,
    refined_difference_161_241 = abs(g$value + check_241$loglik),
    refined_convergence = fit$convergence,
    distance_to_A = dist[["A"]], distance_to_B = dist[["B"]], basin = basin,
    native_to_refined_distance = aligned_distance(p, q),
    maximum_residual_sd = max(sqrt(colSums(p$loading^2))))
  print(rows[[start]])
}
summary <- do.call(rbind, rows)
grid_spread <- function(idx) if (length(idx) < 2L) NA_real_ else
  100 * max(apply(simplify2array(predictions[idx]), 1:2, function(v) diff(range(v))))
basins <- do.call(rbind, lapply(split(summary, summary$basin), function(d) data.frame(
  basin = d$basin[1], starts = paste(d$start, collapse = " "), count = nrow(d),
  fresh_count = sum(d$origin == "fresh"),
  best_native_score = max(d$native_penalised_score),
  native_score_spread = diff(range(d$native_penalised_score)),
  refined_score_spread = diff(range(d$refined_penalised_score)),
  grid_prediction_spread_pp = grid_spread(d$start))))
best_basin <- basins$basin[which.max(basins$best_native_score)]
in_best <- summary$start[summary$basin == best_basin]
across <- data.frame(all_native_score_spread = diff(range(summary$native_penalised_score)),
  all_grid_prediction_spread_pp = grid_spread(seq_len(nrow(summary))),
  best_basin = best_basin, best_basin_count = length(in_best),
  best_basin_native_spread = diff(range(summary$native_penalised_score[in_best])),
  best_basin_grid_spread_pp = grid_spread(in_best),
  passes_within_best_basin = diff(range(summary$native_penalised_score[in_best])) <= .1 &&
    (length(in_best) < 2L || grid_spread(in_best) <= 1),
  best_basin_reached_more_than_once = length(in_best) >= 2L)
# Declared selection rule: the native endpoint with the highest precise
# penalised training objective, provided its basin was reached by at least
# two independent starts and the within-basin checks pass.
selected <- summary$start[which.max(summary$native_penalised_score)]
selection <- data.frame(selected_start = selected, selected_basin = summary$basin[selected],
  selected_native_score = summary$native_penalised_score[selected],
  rule_satisfied = across$passes_within_best_basin && across$best_basin_reached_more_than_once &&
    summary$basin[selected] == best_basin)
print(basins); print(across); print(selection)
write.csv(summary, file.path(out_dir, "multistart-summary.csv"), row.names = FALSE)
write.csv(basins, file.path(out_dir, "multistart-basins.csv"), row.names = FALSE)
write.csv(cbind(across, selection), file.path(out_dir, "multistart-assessment.csv"), row.names = FALSE)
saveRDS(list(summary = summary, basins = basins, across = across, selection = selection,
             predictions = predictions, refined = refined, session = sessionInfo()),
        file.path(out_dir, "multistart-check.rds"))
