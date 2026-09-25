args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "study.R"))
input <- make_input(simulate_community(1, "traits"), 100, 10, "linear", TRUE)
library(gllvm)
fit <- gllvm(y = input$y, X = input$x, TR = input$traits,
  formula = y ~ environment_1 + environment_2 +
    (environment_1 + environment_2):(drought_tolerance + irrelevant_trait),
  randomX = ~ environment_1 + environment_2,
  family = binomial("logit"), link = "logit", num.lv = 2, method = "VA",
  scale.X = FALSE, sd.errors = TRUE, seed = 26100101,
  control.start = list(starting.val = "res", n.init = 1),
  control = list(maxit = 6000, max.iter = 6000))
saveRDS(fit, file.path(args[2], "gllvm-traits.rds"))
print(names(fit)); str(fit$params); print(fit$fourth.corner)
print(fit$convergence); print(fit$sd); print(fit$randomX)
print(args(gllvm:::predict.gllvm))
pr <- predict(fit, newX = input$x, type = "link", level = 0)
saveRDS(pr, file.path(args[2], "gllvm-traits-predict.rds"))
print(dim(pr)); print(head(pr[, 1:3])); print(names(fit))
