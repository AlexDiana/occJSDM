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
  base <- file.path(root, "stability/fits", sprintf("sjSDM-low-rate-start-%d", i))
  f <- readRDS(paste0(base, ".rds"))
  r <- readRDS(paste0(base, "-record.rds"))
  stopifnot(isTRUE(r$ok), identical(f$input, input), length(f$history) == 1500,
            r$settings$data_md5 == unname(tools::md5sum(file.path(root, "inputs/training.rds"))))
  p <- f$parameters
  stopifnot(max(abs((cbind(1, as.matrix(input$x)) %*% p$beta)*.999999+.0000005-f$native_raw)) < 1e-10)
  a <- evaluate(c(p$beta, p$loading))
  b <- joint_integration(p, input$x, input$y, 241L, clamp = TRUE)
  predictions[[i]] <- point_marginal(p, fixed_grid())
  trajectory <- read.csv(paste0(base, "-trajectory.csv"))
  rows[[i]] <- data.frame(start = i, checked_loglik = -a$value,
    difference_121_241 = abs(a$value + b$loglik),
    max_gradient = max(abs(a$gradient)),
    improvement = (-a$value) - trajectory$loglik[1],
    late_improvement = tail(trajectory$loglik,1)-tail(trajectory$loglik,2)[1],
    maximum_residual_sd = max(sqrt(colSums(p$loading^2))))
  saveRDS(list(parameters = p, checked = a, finer = b),
          file.path(root, "stability/checks", sprintf("low-rate-%d.rds", i)))
}
rows <- do.call(rbind, rows)
spread <- max(apply(simplify2array(predictions), 1:2, function(p) diff(range(p))))
summary <- list(rows = rows, likelihood_spread = diff(range(rows$checked_loglik)),
                prediction_spread_pp = 100*spread,
                passed = diff(range(rows$checked_loglik)) <= .1 && spread <= .01)
saveRDS(summary, file.path(root, "stability/checks/low-rate-summary.rds"))
write.csv(rows, file.path(root, "stability/checks/low-rate-summary.csv"), row.names = FALSE)
print(summary)
