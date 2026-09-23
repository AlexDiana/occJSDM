args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
files <- list.files(file.path(root, "fits"), pattern = "([0-9]|res)[.]rds$", full.names = TRUE)
for (path in files) {
  label <- sub("[.]rds$", "", basename(path))
  package <- strsplit(label, "-", fixed = TRUE)[[1]][1]
  output <- file.path(root, "checks", label)
  if (file.exists(paste0(output, "-summary.rds"))) next
  fit <- readRDS(path)
  cat("Checking", label, "\n")
  if (package %in% c("occJSDM", "Hmsc")) {
    p <- bayesian_parameters(fit, package)
    beta_d <- mcmc_diagnostics(p$beta, paste0("beta_", 1:30), "coefficient", p$n, p$chains)
    covariance <- array(NA_real_, c(10, 10, p$n, p$chains))
    for (c in seq_len(p$chains)) for (i in seq_len(p$n)) {
      covariance[, , i, c] <- crossprod(p$loading[, , i, c])
    }
    mask <- lower.tri(matrix(0, 10, 10), diag = TRUE)
    covariance <- matrix(covariance, 100)[as.vector(mask), ]
    cov_d <- mcmc_diagnostics(covariance, paste0("covariance_", 1:55), "covariance", p$n, p$chains)
    probabilities <- bayesian_marginal(p, fixed_grid(), 31L)
    refined <- bayesian_marginal(p, fixed_grid(), 61L)
    integration_error <- max(abs(probabilities - refined))
    stopifnot(integration_error < 1e-4)
    pred_d <- mcmc_diagnostics(refined, paste0("grid_probability_", 1:50), "new_site_probability", p$n, p$chains)
    training <- bayesian_training(p, input$x)
    train_d <- mcmc_diagnostics(training, paste0("training_probability_", 1:1000), "training_probability", p$n, p$chains)
    diagnostics <- rbind(beta_d, cov_d, pred_d, train_d)
    diagnostics$pass <- with(diagnostics, is.finite(rhat) & rhat <= 1.01 &
                              bulk_ess >= 400 & tail_ess >= 400)
    write.csv(diagnostics, paste0(output, "-diagnostics.csv"), row.names = FALSE)
    summary <- list(package = package, passed = all(diagnostics$pass),
                    max_rhat = max(diagnostics$rhat), min_bulk_ess = min(diagnostics$bulk_ess),
                    min_tail_ess = min(diagnostics$tail_ess), failed = sum(!diagnostics$pass),
                    integration_error = integration_error,
                    grid_prediction = apply(refined, 1:2, mean))
    saveRDS(list(parameters = p, training_mean = apply(training, 1:2, mean)),
            paste0(output, "-parameters.rds"), compress = FALSE)
  } else {
    p <- point_parameters(fit, package)
    predictions <- point_marginal(p, fixed_grid())
    summary <- list(package = package, grid_prediction = predictions,
                    max_abs_beta = max(abs(p$beta)),
                    max_residual_sd = max(sqrt(colSums(p$loading^2))))
    if (package == "gllvm") {
      summary$convergence <- fit$convergence
      summary$native_loglik <- fit$logL
      summary$max_gradient <- max(abs(fit$pilot_gradient))
    } else {
      stopifnot(max(abs((cbind(1, as.matrix(input$x)) %*% p$beta) * 0.999999 + 0.0000005 - fit$pilot_native_raw)) < 1e-10)
      summary$first_loss <- fit$history[1]
      summary$mid_loss <- mean(fit$history[(length(fit$history)-199):(length(fit$history)-100)])
      summary$last_loss <- mean(tail(fit$history, 100))
      summary$native_loss <- fit$logLik[[1]]
      q1 <- joint_integration(p, input$x, input$y, 61L, clamp = TRUE)
      q2 <- joint_integration(p, input$x, input$y, 121L, clamp = TRUE)
      summary$joint_difference <- abs(q1$loglik - q2$loglik)
      summary$conditional_difference <- max(abs(q1$conditional - q2$conditional))
      summary$checked_loglik <- q2$loglik
      summary$quadrature_nodes <- 121L
      saveRDS(q2, paste0(output, "-joint.rds"))
    }
    saveRDS(p, paste0(output, "-parameters.rds"))
  }
  saveRDS(summary, paste0(output, "-summary.rds"))
  print(summary[setdiff(names(summary), "grid_prediction")])
}
