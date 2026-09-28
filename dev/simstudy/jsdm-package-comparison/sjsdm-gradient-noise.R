args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
suppressPackageStartupMessages(library(sjSDM))
input <- readRDS(file.path(root, "inputs/training.rds"))
evaluate <- make_joint_gradient(input$x, input$y, 121L)
rows <- list()
for (start in 1:3) {
  parent <- readRDS(file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-polish-start-%d.rds", start)))
  p <- parent$parameters
  theta <- c(p$beta, p$loading)
  precise <- evaluate(theta)
  precise$value <- precise$value + 100*.0001*sum(theta^2)/2
  precise$gradient <- precise$gradient + 100*.0001*theta
  direction <- precise$gradient / sqrt(sum(precise$gradient^2))
  fit <- sjSDM(Y = input$y, env = linear(input$x, lambda = 0),
    biotic = bioticStruct(df = 2, lambda = 0), family = binomial("logit"),
    device = "cpu", dtype = "float64", iter = 1L, sampling = 100L,
    step_size = 100L, learning_rate = .0002, parallel = 0L,
    control = sjSDMControl(optimizer = RMSprop(weight_decay = .0001),
                          scheduler = 0, early_stopping_training = 0),
    seed = 26092700L + start*100L, se = FALSE, verbose = FALSE)
  model <- fit$model
  model$set_env_weights(list(reticulate::r_to_py(t(p$beta))$copy()))
  model$set_sigma(reticulate::r_to_py(t(p$loading))$copy())
  stopifnot(max(abs(t(sjSDM:::force_r(model$env_weights)[[1]]) - p$beta)) == 0,
            max(abs(t(sjSDM:::force_r(model$get_sigma)) - p$loading)) == 0)
  python <- reticulate::py
  python$diagnostic_model <- model
  python$diagnostic_x <- as.matrix(cbind(1,input$x))
  python$diagnostic_y <- input$y
  python$diagnostic_direction <- direction
  python$diagnostic_seed <- 26092700L + start*100L
  path <- file.path(root, "stability-resolution/checks", sprintf("gradient-noise-%d.csv", start))
  stopifnot(!file.exists(path))
  python$diagnostic_output <- path
  reticulate::py_run_file(file.path(root,"stability-resolution/scripts/sjsdm-gradient-noise.py"))
  values <- read.csv(path)
  for (sampling in unique(values$sampling)) {
    d <- values[values$sampling == sampling,]
    rows[[paste(start,sampling)]] <- data.frame(start = start, sampling = sampling,
      evaluations = nrow(d), precise_objective = precise$value,
      mean_objective = mean(d$penalised_loss), objective_se = sd(d$penalised_loss)/sqrt(nrow(d)),
      precise_directional_gradient = sum(precise$gradient*direction),
      mean_directional_gradient = mean(d$directional_gradient),
      gradient_se = sd(d$directional_gradient)/sqrt(nrow(d)))
  }
  write.csv(do.call(rbind,rows),file.path(root,"stability-resolution/checks/gradient-noise-summary.csv"),row.names=FALSE)
  saveRDS(list(parameters=p,precise=precise,direction=direction,session=sessionInfo()),
          file.path(root,"stability-resolution/checks",sprintf("gradient-noise-%d.rds",start)))
  cat("Completed native gradient sampling at endpoint", start, "\n")
}
print(do.call(rbind, rows))
