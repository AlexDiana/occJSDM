# Curvature check at the deterministic penalised-reference endpoints.
# Question: are the two distinct stationary solutions local maxima of the
# penalised training log likelihood, or is one a saddle? The two factors can
# be rotated without changing the model or the penalty, so exactly one flat
# direction is expected at any stationary point; it is identified and set aside
# before the remaining curvature is judged. Deterministic quadrature only; no
# test outcomes or generating truth are read. Nothing is fitted here.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
out_dir <- file.path(root, "stability-resolution/curvature")
dir.create(out_dir, showWarnings = FALSE)
stopifnot(!file.exists(file.path(out_dir, "curvature-summary.csv")))
n_sites <- nrow(input$x)
decay <- 0.0001
n_beta <- 3L * ncol(input$y)

penalised <- function(evaluate) function(theta) {
  v <- evaluate(theta)
  list(value = v$value + n_sites * decay * sum(theta^2) / 2,
       gradient = v$gradient + n_sites * decay * theta)
}

# Central differences of the analytic gradient, then symmetrised.
hessian <- function(f, theta, step) {
  k <- length(theta)
  H <- matrix(NA_real_, k, k)
  for (j in seq_len(k)) {
    d <- theta * 0; d[j] <- step
    H[, j] <- (f(theta + d)$gradient - f(theta - d)$gradient) / (2 * step)
  }
  list(H = (H + t(H)) / 2, asymmetry = max(abs(H - t(H))))
}

# Generator of the factor rotation: beta unchanged, each species' loading
# pair (l1, l2) moves along (-l2, l1).
rotation_generator <- function(theta) {
  loading <- matrix(theta[-seq_len(n_beta)], 2)
  g <- c(theta[seq_len(n_beta)] * 0, as.vector(rbind(-loading[2, ], loading[1, ])))
  g / sqrt(sum(g^2))
}

# Orthogonal 2x2 map aligning loading matrix B onto A in least squares.
procrustes <- function(A, B) {
  s <- svd(A %*% t(B))
  s$u %*% t(s$v)
}

theta_of <- function(p) c(p$beta, p$loading)
refs <- lapply(1:3, function(i) readRDS(file.path(root, "stability-resolution/fits",
  sprintf("penalised-reference-start-%d.rds", i))))

started <- Sys.time()
rows <- list(); eigen_out <- list()
for (i in 1:3) {
  theta <- theta_of(refs[[i]]$parameters)
  for (nodes in c(121L, 161L)) {
    f <- penalised(make_joint_gradient(input$x, input$y, nodes))
    g0 <- f(theta)
    step <- 1e-4
    h <- hessian(f, theta, step)
    e <- eigen(h$H, symmetric = TRUE)
    gen <- rotation_generator(theta)
    idx_flat <- which.min(abs(e$values))
    alignment <- abs(sum(e$vectors[, idx_flat] * gen))
    Hg <- h$H %*% gen
    others <- e$values[-idx_flat]
    rows[[length(rows) + 1L]] <- data.frame(start = i, nodes = nodes,
      penalised_score = -g0$value, max_gradient = max(abs(g0$gradient)),
      hessian_asymmetry = h$asymmetry, step = step,
      flat_eigenvalue = e$values[idx_flat], flat_alignment_with_rotation = alignment,
      max_abs_H_times_generator = max(abs(Hg)),
      smallest_other_eigenvalue = min(others), largest_eigenvalue = max(e$values),
      negative_eigenvalues = sum(others < -1e-6),
      near_zero_other_eigenvalues = sum(abs(others) < 1e-4),
      condition_excluding_flat = max(e$values) / min(others),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    eigen_out[[paste(i, nodes)]] <- list(values = e$values, flat_vector = e$vectors[, idx_flat],
      smallest_other_vector = e$vectors[, order(e$values)[if (idx_flat == order(e$values)[1]) 2 else 1]],
      generator = gen, H = h$H)
    print(rows[[length(rows)]])
  }
  # Step-size sensitivity at 121 nodes for the smallest non-flat eigenvalue.
  f <- penalised(make_joint_gradient(input$x, input$y, 121L))
  alt <- eigen(hessian(f, theta, 3e-4)$H, symmetric = TRUE)$values
  alt_flat <- which.min(abs(alt))
  rows[[length(rows)]]$smallest_other_eigenvalue_step_3e4 <- min(alt[-alt_flat])
  rows[[length(rows) - 1L]]$smallest_other_eigenvalue_step_3e4 <- NA_real_
}
summary <- do.call(rbind, rows)

# Are solutions 2 and 3 the same point up to factor rotation?
p2 <- refs[[2]]$parameters; p3 <- refs[[3]]$parameters
Q <- procrustes(p2$loading, p3$loading)
same <- data.frame(max_beta_difference_2_3 = max(abs(p2$beta - p3$beta)),
  max_rotated_loading_difference_2_3 = max(abs(p2$loading - Q %*% p3$loading)),
  max_implied_covariance_difference_2_3 = max(abs(crossprod(p2$loading) - crossprod(p3$loading))),
  rotation_determinant = det(Q))
print(same)

# Objective along the straight path between the two distinct solutions, after
# rotating solution 1's factors onto solution 2's. A rise between them means
# separate basins; a monotone path means one endpoint was not stationary.
p1 <- refs[[1]]$parameters
Q1 <- procrustes(p2$loading, p1$loading)
theta1 <- c(p1$beta, Q1 %*% p1$loading)
theta2 <- theta_of(p2)
f <- penalised(make_joint_gradient(input$x, input$y, 121L))
stopifnot(abs(f(theta1)$value - f(theta_of(p1))$value) < 1e-8)  # rotation invariance
path <- data.frame(t = seq(-0.25, 1.25, by = 0.05))
path$penalised_score <- vapply(path$t, function(t) -f(theta2 + t * (theta1 - theta2))$value, numeric(1))
path$distance_from_solution_2 <- path$t * sqrt(sum((theta1 - theta2)^2))
interior <- path[path$t > 0 & path$t < 1, ]
barrier <- data.frame(score_solution_2 = path$penalised_score[path$t == 0],
  score_solution_1 = path$penalised_score[path$t == 1],
  minimum_interior_score = min(interior$penalised_score),
  t_at_minimum = interior$t[which.min(interior$penalised_score)],
  barrier_below_lower_endpoint = min(path$penalised_score[path$t %in% c(0, 1)]) - min(interior$penalised_score),
  parameter_distance = sqrt(sum((theta1 - theta2)^2)))
print(path); print(barrier)

# Which species differ between the two solutions?
species_diff <- data.frame(species = colnames(input$y),
  beta_difference = sqrt(colSums((p1$beta - p2$beta)^2)),
  residual_sd_solution_1 = sqrt(colSums(p1$loading^2)),
  residual_sd_solution_2 = sqrt(colSums(p2$loading^2)))
print(species_diff)

write.csv(summary, file.path(out_dir, "curvature-summary.csv"), row.names = FALSE)
write.csv(path, file.path(out_dir, "path-profile.csv"), row.names = FALSE)
write.csv(species_diff, file.path(out_dir, "species-differences.csv"), row.names = FALSE)
saveRDS(list(summary = summary, eigen = eigen_out, same_solution_2_3 = same, path = path,
             barrier = barrier, species_diff = species_diff, theta1_aligned = theta1,
             theta2 = theta2, session = sessionInfo()),
        file.path(out_dir, "curvature.rds"))
cat("Curvature check finished in", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "seconds\n")
