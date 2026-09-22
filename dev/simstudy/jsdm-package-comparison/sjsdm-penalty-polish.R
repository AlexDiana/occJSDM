# Continue a saved native sjSDM fit without changing model, data or penalties.
# Always launch a frozen copy with source(), not an editable live script.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
start <- as.integer(args[3])
stopifnot(start %in% 1:3)
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
suppressPackageStartupMessages(library(sjSDM))
input <- readRDS(file.path(root, "inputs/training.rds"))
parent_path <- file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-start-%d.rds", start))
parent <- readRDS(parent_path)
label <- sprintf("sjSDM-default-penalty-polish-start-%d", start)
out <- file.path(root, "stability/fits", label)
stopifnot(!file.exists(paste0(out, "-record.rds")))
settings <- list(parent = basename(parent_path), parent_md5 = unname(tools::md5sum(parent_path)),
  seed = 26092640L + start, epochs = 1000L, block = 250L,
  learning_rate = .0002, sampling = 2000L, optimiser = "RMSprop", weight_decay = 0.0001,
  data_md5 = unname(tools::md5sum(file.path(root, "inputs/training.rds"))))
started <- Sys.time()
warnings <- character()
result <- tryCatch(withCallingHandlers({
  # One construction epoch creates the same package model, then is discarded.
  fit <- sjSDM(Y = input$y, env = linear(input$x, lambda = 0),
    biotic = bioticStruct(df = 2, lambda = 0), family = binomial("logit"),
    device = "cpu", dtype = "float64", iter = 1L, sampling = 2000L,
    step_size = 100L, learning_rate = settings$learning_rate, parallel = 0L,
    control = sjSDMControl(optimizer = RMSprop(weight_decay = 0.0001),
                          scheduler = 0, early_stopping_training = 0),
    seed = settings$seed, se = FALSE, verbose = FALSE)
  model <- fit$model
  model$set_env_weights(lapply(parent$weights, function(w) reticulate::r_to_py(w)$copy()))
  model$set_sigma(reticulate::r_to_py(parent$sigma)$copy())
  python <- reticulate::py
  python$continuation_model <- model
  # The weight setter replaces Parameter objects. Point a fresh optimiser at
  # those restored objects; the original optimiser still holds the old objects.
  reticulate::py_run_string(paste(
    "import torch",
    "m = continuation_model",
    "m.optimizer = torch.optim.RMSprop(list(m.env.parameters()) + [m.sigma], lr=0.0002, alpha=0.99, eps=1e-8, momentum=0.1, centered=False, weight_decay=0.0001)",
    "assert {id(p) for g in m.optimizer.param_groups for p in g['params']} == {id(p) for p in list(m.env.parameters()) + [m.sigma]}",
    "assert all(p.dtype == torch.float64 for p in list(m.env.parameters()) + [m.sigma])",
    "assert str(m.device) == 'cpu' and m.df == 2 and m._loss_function.__name__ == 'torch_tmp'",
    sep = "\n"))
  current_parameters <- function() list(
    beta = t(sjSDM:::force_r(model$env_weights)[[1]]),
    loading = t(sjSDM:::force_r(model$get_sigma)))
  p <- current_parameters()
  expected <- point_parameters(parent, "sjSDM")
  stopifnot(max(abs(p$beta-expected$beta)) == 0,
            max(abs(p$loading-expected$loading)) == 0)
  native <- predict(fit, newdata = input$x, type = "raw")
  stopifnot(max(abs(native - ((cbind(1, as.matrix(input$x)) %*% p$beta)*.999999+.0000005))) < 1e-10)
  evaluate <- make_joint_gradient(input$x, input$y, 61L)
  history <- numeric()
  rows <- list()
  for (block in 0:4) {
    if (block > 0) {
      model$fit(as.matrix(cbind(1, input$x)), input$y,
                batch_size = 100L, epochs = 250L, sampling = 2000L,
                parallel = 0L, early_stopping_training = -1L, verbose = FALSE)
      history <- c(history, as.vector(sjSDM:::force_r(model$history)))
    }
    p <- current_parameters()
    diagnostic <- evaluate(c(p$beta, p$loading))
    rows[[block+1L]] <- data.frame(start = start, epoch = block*250L,
      loglik = -diagnostic$value, max_gradient = max(abs(diagnostic$gradient)),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    print(rows[[block+1L]])
    saveRDS(list(parameters = p, history = history, diagnostic = diagnostic,
                 native_raw = predict(fit, newdata = input$x, type = "raw")),
             sprintf("%s-epoch-%04d.rds", out, block*250L))
    write.csv(do.call(rbind, rows), paste0(out, "-trajectory.csv"), row.names = FALSE)
  }
  saveRDS(list(parameters = p, history = history, settings = settings,
               input = input, native_raw = predict(fit, newdata = input$x, type = "raw")),
           paste0(out, ".rds"))
  list(ok = TRUE)
}, warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
}), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
saveRDS(c(result, list(settings = settings, warnings = unique(warnings),
                      started = started, ended = Sys.time(), session = sessionInfo())),
        paste0(out, "-record.rds"))
if (!isTRUE(result$ok)) stop(result$error)
cat("Saved continuation with clean script completion:", label, "\n")
