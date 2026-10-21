# sjsdm-multistart.R adapted to rerun one recorded start under Doug's sjSDM
# fork release v0.1.0, the version Lesson 5 tells readers to install. The
# recorded starts ran under v0.2.1 with SJSDM_MOJO_BACKEND=0. Differences
# from sjsdm-multistart.R, and nothing else:
# - v0.1.0 names the PyTorch Monte Carlo loss `tmp`. v0.2.1 renamed the same
#   function `torch_tmp` when it added the opt-in Mojo loss; on this path it
#   only adds int() casts to the noise dimensions. Both name checks accept
#   only `tmp`.
# - The settings, seeds and training-data MD5 must equal the archived record
#   of the same start, and sjSDM must load from the dedicated v0.1.0 library.
# - Outputs carry a v0.1.0 label, so nothing archived can be overwritten.
# SJSDM_MOJO_BACKEND does nothing in v0.1.0; it is still required to be "0",
# as in the recorded runs. Run with the v0.1.0 environment's run-r, one
# process per start, executing an immutable copy; never edit a running script.
#   run-r sjsdm-multistart-v010.R ROOT CODE START
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
start <- as.integer(args[3])
stopifnot(start %in% 4:12, Sys.getenv("SJSDM_MOJO_BACKEND") == "0")
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
archived <- readRDS(file.path(code, "stability-resolution-archive/stability-resolution/multistart",
                              sprintf("sjSDM-multistart-start-%d-record.rds", start)))
suppressPackageStartupMessages(library(sjSDM))
stopifnot(grepl("/sjsdm-v0\\.1\\.0-env-[0-9]+/library/sjSDM$", find.package("sjSDM")))
input <- readRDS(file.path(root, "inputs/training.rds"))
label <- sprintf("sjSDM-v0.1.0-multistart-start-%d", start)
out <- file.path(root, "stability-resolution/multistart", label)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
stopifnot(!file.exists(paste0(out, "-record.rds")))
settings <- list(start = start, fit_seed = 26092800L + start, polish_seed = 26092900L + start,
  epochs = 3000L, learning_rate = 0.002, polish_epochs = 1000L, polish_block = 250L,
  polish_learning_rate = 0.0002, sampling = 2000L, batch = 100L,
  optimiser = "RMSprop", weight_decay = 0.0001,
  data_md5 = unname(tools::md5sum(file.path(root, "inputs/training.rds"))))
stopifnot(identical(settings, archived$settings[names(settings)]))
settings$sjsdm_release <- "v0.1.0"
settings$sjsdm_library <- find.package("sjSDM")
settings$sjsdm_version <- as.character(utils::packageVersion("sjSDM"))
started <- Sys.time()
warnings <- character()
evaluate <- make_joint_gradient(input$x, input$y, 61L)
penalised <- function(theta) {
  v <- evaluate(theta)
  list(loglik = -v$value, penalised_score = -v$value - nrow(input$x) * settings$weight_decay * sum(theta^2) / 2,
       max_penalised_gradient = max(abs(v$gradient + nrow(input$x) * settings$weight_decay * theta)))
}
result <- tryCatch(withCallingHandlers({
  fit <- sjSDM(Y = input$y, env = linear(data = input$x, lambda = 0),
    biotic = bioticStruct(df = 2, lambda = 0),
    family = binomial("logit"), device = "cpu", dtype = "float64",
    iter = settings$epochs, sampling = settings$sampling, step_size = settings$batch,
    learning_rate = settings$learning_rate, parallel = 0L,
    control = sjSDMControl(optimizer = RMSprop(weight_decay = settings$weight_decay),
                          scheduler = 0, early_stopping_training = 0),
    seed = settings$fit_seed, se = FALSE, verbose = FALSE)
  model <- fit$model
  settings$backend <- reticulate::py_to_r(model$`_loss_function`$`__name__`)
  settings$python <- reticulate::py_config()$python
  settings$torch <- reticulate::py_to_r(reticulate::import("torch")$`__version__`)
  stopifnot(settings$backend == "tmp",
    as.character(reticulate::py_str(model$device)) == "cpu",
    as.integer(reticulate::py_to_r(model$df)) == 2L,
    fit$settings$control$optimizer$params$weight_decay == settings$weight_decay,
    identical(fit$data$Y, input$y),
    max(abs(fit$data$X - cbind(1, as.matrix(input$x)))) < 1e-12,
    length(fit$history) == settings$epochs)
  current_parameters <- function() list(
    beta = t(sjSDM:::force_r(model$env_weights)[[1]]),
    loading = t(sjSDM:::force_r(model$get_sigma)))
  p <- current_parameters()
  stopifnot(max(abs(p$beta - point_parameters(fit, "sjSDM")$beta)) == 0,
            max(abs(p$loading - point_parameters(fit, "sjSDM")$loading)) == 0)
  rows <- list(data.frame(start = start, stage = "fit", epoch = settings$epochs,
    penalised(c(p$beta, p$loading)),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))))
  print(rows[[1]])
  saveRDS(list(parameters = p, history = as.vector(fit$history), settings = settings),
          paste0(out, "-fit.rds"))
  # Low-step continuation with a fresh optimiser at the smaller learning rate,
  # exactly as in sjsdm-multistart.R.
  python <- reticulate::py
  python$multistart_model <- model
  reticulate::py_run_string(paste(
    "import torch",
    "m = multistart_model",
    sprintf("torch.manual_seed(%d)", settings$polish_seed),
    "m.optimizer = torch.optim.RMSprop(list(m.env.parameters()) + [m.sigma], lr=0.0002, alpha=0.99, eps=1e-8, momentum=0.1, centered=False, weight_decay=0.0001)",
    "assert {id(p) for g in m.optimizer.param_groups for p in g['params']} == {id(p) for p in list(m.env.parameters()) + [m.sigma]}",
    "assert all(p.dtype == torch.float64 for p in list(m.env.parameters()) + [m.sigma])",
    "assert str(m.device) == 'cpu' and m.df == 2 and m._loss_function.__name__ == 'tmp'",
    sep = "\n"))
  polish_history <- numeric()
  for (block in seq_len(settings$polish_epochs / settings$polish_block)) {
    model$fit(as.matrix(cbind(1, input$x)), input$y,
              batch_size = settings$batch, epochs = settings$polish_block,
              sampling = settings$sampling, parallel = 0L,
              early_stopping_training = -1L, verbose = FALSE)
    polish_history <- c(polish_history, as.vector(sjSDM:::force_r(model$history)))
    p <- current_parameters()
    rows[[length(rows) + 1L]] <- data.frame(start = start, stage = "polish",
      epoch = block * settings$polish_block, penalised(c(p$beta, p$loading)),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    print(rows[[length(rows)]])
    write.csv(do.call(rbind, rows), paste0(out, "-trajectory.csv"), row.names = FALSE)
  }
  saveRDS(list(parameters = p, polish_history = polish_history, settings = settings,
               native_raw = predict(fit, newdata = input$x, type = "raw")),
          paste0(out, ".rds"))
  write.csv(do.call(rbind, rows), paste0(out, "-trajectory.csv"), row.names = FALSE)
  list(ok = TRUE)
}, warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
}), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
saveRDS(c(result, list(settings = settings, warnings = unique(warnings), started = started,
                      ended = Sys.time(),
                      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
                      session = sessionInfo())),
        paste0(out, "-record.rds"))
if (!isTRUE(result$ok)) stop(result$error)
cat("Completed sjSDM v0.1.0 multistart fit and continuation:", label, "\n")
