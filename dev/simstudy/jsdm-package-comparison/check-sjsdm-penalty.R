args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
evaluate <- make_joint_gradient(input$x, input$y, 121L)
rows <- list()
predictions <- list()
for (i in 1:3) {
  label <- sprintf("sjSDM-default-penalty-start-%d", i)
  base <- file.path(root, "stability/fits", label)
  f <- readRDS(paste0(base, ".rds"))
  r <- readRDS(paste0(base, "-record.rds"))
  stopifnot(isTRUE(r$ok), identical(f$data$Y, input$y),
    max(abs(f$data$X - cbind(1, as.matrix(input$x)))) < 1e-12,
    length(f$history) == 3000L, f$seed == 26092340L + i,
    r$settings$weight_decay == .0001,
    f$settings$control$optimizer$params$weight_decay == .0001,
    f$settings$biotic$df == 2L, f$settings$biotic$l1_cov == 0,
    f$settings$biotic$l2_cov == 0, r$settings$backend == "torch_tmp")
  p <- point_parameters(f, "sjSDM")
  theta <- c(p$beta, p$loading)
  stopifnot(max(abs((cbind(1, as.matrix(input$x)) %*% p$beta)*.999999+.0000005-f$pilot_native_raw)) < 1e-10)
  a <- evaluate(theta)
  b <- joint_integration(p, input$x, input$y, 241L, clamp = TRUE)
  penalty <- nrow(input$x)*.0001*sum(theta^2)/2
  penalised_gradient <- a$gradient + nrow(input$x)*.0001*theta
  predictions[[i]] <- point_marginal(p, fixed_grid())
  rows[[i]] <- data.frame(start = i, loglik = -a$value,
    penalty = penalty, penalised_score = -a$value-penalty,
    difference_121_241 = abs(a$value+b$loglik),
    max_penalised_gradient = max(abs(penalised_gradient)),
    mean_loss_previous_100 = mean(tail(f$history, 200)[1:100]),
    mean_loss_final_100 = mean(tail(f$history, 100)),
    maximum_residual_sd = max(sqrt(colSums(p$loading^2))))
  stopifnot(abs(a$value+b$loglik) < .001)
  saveRDS(list(parameters = p, gradient = a, penalty = penalty,
               penalised_gradient = penalised_gradient, refined = b),
    file.path(root, "stability/checks", paste0(label, ".rds")))
}
rows <- do.call(rbind, rows)
spread <- max(apply(simplify2array(predictions), 1:2, function(p) diff(range(p))))
summary <- list(rows = rows, objective_spread = diff(range(rows$penalised_score)),
                prediction_spread_pp = 100*spread,
                passed = diff(range(rows$penalised_score)) <= .1 && spread <= .01)
saveRDS(summary, file.path(root, "stability/checks/default-penalty-summary.rds"))
write.csv(rows, file.path(root, "stability/checks/default-penalty-summary.csv"), row.names = FALSE)
print(summary)
