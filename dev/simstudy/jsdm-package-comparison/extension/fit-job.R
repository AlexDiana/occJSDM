# One immutable worker per package/data combination. Never opens truth/.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1]); code <- normalizePath(args[2]); job_number <- as.integer(args[3])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "extension", "math.R"))
manifest <- read.csv(file.path(root, "manifest.csv"))
job <- manifest[manifest$job_number == job_number, ]
stopifnot(nrow(job) == 1)
input_path <- file.path(root, "inputs", paste0(job$input, ".rds"))
stopifnot(unname(tools::md5sum(input_path)) == job$input_md5)
input <- readRDS(input_path)
# Test predictors are not used for fitting, selection or diagnostics.
input$test_x <- NULL
out <- file.path(root, "jobs", job$job)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(!file.exists(file.path(out, "result.rds")))
suppressPackageStartupMessages(library(job$package, character.only = TRUE))
base_seed <- 26100000L + 100L * job_number
started <- Sys.time()
warning_messages <- character()
formula_env <- reformulate(names(input$x))
grid <- input$x[c(1, 21, 41, 61, 81), , drop = FALSE]
trait_model <- job$scenario == "traits"
records <- list(); candidates <- list()

run_bayesian <- function(attempt) {
  warmup <- if (attempt == 1) 2000L else 8000L
  draws <- if (attempt == 1) 4000L else 16000L
  set.seed(base_seed + attempt)
  if (job$package == "occJSDM") {
    fit <- occJSDM::runOccJSDM(list(info = input$x, OTU = input$y, traits = input$traits),
      occCovariates = names(input$x), listParams = list(n_factors = 2, n_lattrait = 0),
      MCMCparams = list(nchain = 4, nburn = warmup, niter = draws, nthin = 1))
    stopifnot(identical(fit$infos$OTU, input$y), fit$infos$model == "binary", fit$infos$ps == 0,
      max(abs(fit$X_psi - as.matrix(input$x))) < 1e-10)
  } else {
    design <- data.frame(site = factor(rownames(input$x), levels = rownames(input$x)))
    random <- Hmsc::setPriors(Hmsc::HmscRandomLevel(units = levels(design$site)), nfMin = 2, nfMax = 2)
    model <- Hmsc::Hmsc(Y = input$y, XData = input$x, XFormula = formula_env, XScale = FALSE,
      TrData = input$traits, TrFormula = if (is.null(input$traits)) NULL else ~ drought_tolerance + irrelevant_trait,
      distr = "probit", studyDesign = design, ranLevels = list(site = random))
    fit <- Hmsc::sampleMcmc(model, samples = draws, transient = warmup, thin = 1,
                          nChains = 4, nParallel = 1, verbose = 1000)
    stopifnot(identical(fit$Y, input$y), max(abs(fit$X - cbind(1, as.matrix(input$x)))) < 1e-10)
  }
  p <- extract_bayesian(fit, job$package)
  diagnostics <- parameter_diagnostics(p, input)
  # All unthinned global draws needed for the stated new-site target are saved.
  # Site-score draws are not needed for this extension's new-site question.
  saveRDS(p, file.path(out, sprintf("attempt-%d-parameters.rds", attempt)))
  write.csv(diagnostics$table, file.path(out, sprintf("attempt-%d-diagnostics.csv", attempt)), row.names = FALSE)
  list(parameters = p, diagnostics = diagnostics,
    settings = list(seed = base_seed + attempt, warmup = warmup, draws = draws, chains = 4,
      factors = 2, latent_traits = 0, link = p$link))
}

run_gllvm <- function(start) {
  options <- list(y = input$y, X = input$x, family = binomial("logit"), link = "logit",
    num.lv = 2, method = "VA", seed = base_seed + start,
    sd.errors = trait_model, control.start = list(starting.val = "res", n.init = 1, jitter.var = .1),
    control = list(maxit = 12000, max.iter = 12000))
  if (trait_model) {
    options$TR <- if (is.null(input$traits)) data.frame(unused_species_index = seq_len(ncol(input$y))) else input$traits
    options$formula <- if (is.null(input$traits)) y ~ environment_1 + environment_2 else
      y ~ environment_1 + environment_2 + (environment_1 + environment_2):(drought_tolerance + irrelevant_trait)
    options$randomX <- ~ environment_1 + environment_2
    options$scale.X <- FALSE
  }
  fit <- do.call(gllvm::gllvm, options)
  stopifnot(identical(fit$y, input$y), max(abs(fit$X - as.matrix(input$x))) < 1e-10)
  if (trait_model && is.null(input$traits)) {
    # gllvm requires a nonconstant TR column to initialise fourth-corner
    # machinery. The placeholder is absent from the formula and coefficients.
    stopifnot(setequal(names(fit$params$B), names(input$x)))
  }
  p <- extract_gllvm(fit, input, trait_model)
  gradient <- fit$TMBfn$gr(fit$TMBfn$env$last.par.best)
  admissible <- isTRUE(as.logical(fit$convergence)) && all(is.finite(gradient)) && max(abs(gradient)) < .01
  saveRDS(list(parameters = p, native_parameters = fit$params, standard_errors = fit$sd,
    traits = fit$TR, fourth_corner = fit$fourth.corner, gradient = gradient,
    native_loglik = fit$logL, convergence = fit$convergence, call = fit$call),
    file.path(out, sprintf("start-%d.rds", start)))
  list(parameters = p, score = fit$logL, acceptable = admissible,
       gradient = max(abs(gradient)), grid = point_marginal(p, grid))
}

run_sjsdm <- function(start) {
  stopifnot(Sys.getenv("SJSDM_MOJO_BACKEND") == "0")
  fit <- sjSDM::sjSDM(Y = input$y, env = sjSDM::linear(data = input$x, lambda = 0),
    biotic = sjSDM::bioticStruct(df = 2, lambda = 0), family = binomial("logit"),
    device = "cpu", dtype = "float64", iter = 3000L, sampling = 2000L,
    step_size = as.integer(nrow(input$x)), learning_rate = .002, parallel = 0L,
    control = sjSDM::sjSDMControl(optimizer = sjSDM::RMSprop(weight_decay = .0001),
      scheduler = 0, early_stopping_training = 0),
    seed = base_seed + start, se = FALSE, verbose = FALSE)
  model <- fit$model
  stopifnot(reticulate::py_to_r(model$`_loss_function`$`__name__`) == "torch_tmp",
            identical(fit$data$Y, input$y),
            max(abs(fit$data$X - cbind(1, as.matrix(input$x)))) < 1e-10)
  python <- reticulate::py
  python$extension_model <- model
  reticulate::py_run_string(paste(
    "import torch", "m = extension_model",
    sprintf("torch.manual_seed(%d)", base_seed + 50L + start),
    "m.optimizer = torch.optim.RMSprop(list(m.env.parameters()) + [m.sigma], lr=0.0002, alpha=0.99, eps=1e-8, momentum=0.1, centered=False, weight_decay=0.0001)", sep = "\n"))
  history <- fit$history
  model$fit(as.matrix(cbind(1, input$x)), input$y, batch_size = as.integer(nrow(input$x)),
    epochs = 1000L, sampling = 2000L, parallel = 0L, early_stopping_training = -1L, verbose = FALSE)
  p <- list(beta = t(sjSDM:::force_r(model$env_weights)[[1]]), loading = t(sjSDM:::force_r(model$get_sigma)))
  q <- joint_integration(p, input$x, input$y, 81L, clamp = TRUE)
  fine <- joint_integration(p, input$x, input$y, 161L, clamp = TRUE)
  difference <- abs(q$loglik - fine$loglik)
  if (difference > .001) {
    finer <- joint_integration(p, input$x, input$y, 241L, clamp = TRUE)
    difference <- abs(fine$loglik - finer$loglik); fine <- finer
  }
  score <- fine$loglik - nrow(input$x) * .0001 * sum(c(p$beta, p$loading)^2) / 2
  saveRDS(list(parameters = p, history = history, continuation_history = sjSDM:::force_r(model$history),
    seed = base_seed + start, continuation_seed = base_seed + 50 + start,
    score = score, integration_difference = difference, backend = "torch_tmp",
    epochs = c(3000L, 1000L), learning_rate = c(.002, .0002), sampling = 2000L,
    batch = nrow(input$x), weight_decay = .0001), file.path(out, sprintf("start-%d.rds", start)))
  list(parameters = p, score = score, acceptable = is.finite(score) && difference <= .001,
       integration_difference = difference, grid = point_marginal(p, grid))
}

stable_selection <- function(candidates) {
  eligible <- which(vapply(candidates, function(x) isTRUE(x$acceptable), logical(1)))
  if (!length(eligible)) return(list(selected = NA_integer_, stable = FALSE))
  best <- eligible[which.max(vapply(candidates[eligible], `[[`, numeric(1), "score"))]
  matches <- vapply(setdiff(eligible, best), function(i) {
    abs(candidates[[i]]$score - candidates[[best]]$score) <= .1 &&
      max(abs(candidates[[i]]$grid - candidates[[best]]$grid)) <= .01
  }, logical(1))
  list(selected = best, stable = any(matches))
}

result <- tryCatch(withCallingHandlers({
  if (job$package %in% c("occJSDM", "Hmsc")) {
    selected <- run_bayesian(1)
    if (!selected$diagnostics$passed) selected <- run_bayesian(2)
    p <- selected$parameters
    list(ok = TRUE, diagnostic_pass = selected$diagnostics$passed,
         parameters = p, diagnostics = selected$diagnostics, settings = selected$settings)
  } else {
    for (start in 1:6) {
      step_start <- Sys.time()
      candidate <- tryCatch(if (job$package == "gllvm") run_gllvm(start) else run_sjsdm(start),
        error = function(e) list(acceptable = FALSE, error = conditionMessage(e)))
      candidates[[start]] <- candidate
      record <- candidate[setdiff(names(candidate), c("parameters", "grid"))]
      record$seconds <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
      records[[start]] <- record
      saveRDS(records, file.path(out, "attempts.rds"))
      if (start >= 2 && !is.null(candidate$error) &&
          identical(candidate$error, candidates[[start - 1L]]$error)) {
        stop("Repeated identical fitting error: ", candidate$error)
      }
      selection <- stable_selection(candidates)
      if (start >= 3 && selection$stable) break
    }
    if (is.na(selection$selected)) stop("No fit satisfied the declared numerical checks; inspect attempts.rds")
    p <- candidates[[selection$selected]]$parameters
    list(ok = TRUE, diagnostic_pass = selection$stable, parameters = p,
         selected_start = selection$selected, attempts = records)
  }
}, warning = function(w) {
  warning_messages <<- c(warning_messages, conditionMessage(w))
  invokeRestart("muffleWarning")
}), error = function(e) list(ok = FALSE, diagnostic_pass = FALSE, error = conditionMessage(e)))
result <- c(result, list(job = job, warnings = unique(warning_messages), started = started,
  finished = Sys.time(), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  input_md5 = job$input_md5, truth_used = FALSE, session = sessionInfo()))
saveRDS(result, file.path(out, "result.partial.rds"))
stopifnot(file.rename(file.path(out, "result.partial.rds"), file.path(out, "result.rds")))
summary <- result[intersect(names(result), c("ok", "diagnostic_pass", "error", "selected_start", "elapsed_seconds", "warnings", "input_md5", "truth_used"))]
summary$job_number <- job_number
summary$result_md5 <- unname(tools::md5sum(file.path(out, "result.rds")))
jsonlite::write_json(summary, file.path(out, "status.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
print(result[intersect(names(result), c("ok", "diagnostic_pass", "error", "selected_start", "elapsed_seconds", "warnings"))])
if (!isTRUE(result$ok)) quit(status = 1L)
