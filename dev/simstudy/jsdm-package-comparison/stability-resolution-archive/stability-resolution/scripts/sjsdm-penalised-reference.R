# Diagnostic deterministic optimisation of the fixed default-penalty arm.
# The loglik summary columns contain PENALISED log scores, for within-arm checks.
# This is not a native sjSDM fit and must never replace one without that label.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
start <- as.integer(args[3])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
parent <- readRDS(file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-polish-start-%d.rds", start)))
label <- sprintf("penalised-reference-start-%d", start)
out <- file.path(root, "stability-resolution/fits", label)
stopifnot(!file.exists(paste0(out,".rds")))
# 81 nodes per factor is a diagnostic approximation; independent 161/241-node
# checks decide whether the endpoint is numerically trustworthy.
f <- make_joint_gradient(input$x, input$y, 81L)
last <- NULL
last_value <- NULL
history <- list()
evaluations <- 0L
started <- Sys.time()
evaluate <- function(theta) {
  if (!identical(theta, last)) {
    value <- f(theta)
    value$unpenalised <- value$value
    value$value <- value$value + nrow(input$x)*.0001*sum(theta^2)/2
    value$gradient <- value$gradient + nrow(input$x)*.0001*theta
    last <<- theta
    last_value <<- value
    evaluations <<- evaluations + 1L
    history[[evaluations]] <<- c(evaluation = evaluations, loglik = -value$value,
      max_gradient = max(abs(value$gradient)),
      residual_sd = max(sqrt(colSums(matrix(theta[-(1:30)], 2, 10)^2))))
  }
  last_value
}
initial <- c(parent$parameters$beta, parent$parameters$loading)
fit <- optim(initial, function(t) evaluate(t)$value, function(t) evaluate(t)$gradient,
             method = "BFGS", control = list(maxit = 300, reltol = 1e-12, trace = 1, REPORT = 20))
p <- list(beta = matrix(fit$par[1:30], 3, 10), loading = matrix(fit$par[-(1:30)], 2, 10))
fine <- make_joint_gradient(input$x, input$y, 161L)(fit$par)
penalty <- nrow(input$x)*.0001*sum(fit$par^2)/2
fine$unpenalised <- fine$value
fine$value <- fine$value + penalty
fine$gradient <- fine$gradient + nrow(input$x)*.0001*fit$par
finer <- joint_integration(p, input$x, input$y, 241L, clamp = TRUE)
summary <- list(start = start, convergence = fit$convergence, counts = fit$counts,
                loglik_81 = -fit$value, loglik_161 = -fine$value,
                loglik_241 = finer$loglik-penalty, max_gradient_161 = max(abs(fine$gradient)),
                maximum_residual_sd = max(sqrt(colSums(p$loading^2))),
                grid_prediction = point_marginal(p, fixed_grid()),
                elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
saveRDS(list(parameters = p, optim = fit, fine = fine, finer = finer,
             summary = summary, history = do.call(rbind, history)), paste0(out, ".rds"))
write.csv(do.call(rbind, history), paste0(out, "-trajectory.csv"), row.names = FALSE)
print(summary[setdiff(names(summary), "grid_prediction")])
