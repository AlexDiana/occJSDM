# Execute in a fresh R process using the isolated run-r launcher.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(4L, 5L, 6L))
root <- normalizePath(args[1])
package <- args[2]
attempt <- as.integer(args[3])
start <- as.integer(args[4])
stopifnot(package %in% c("occJSDM", "Hmsc", "gllvm", "sjSDM"), attempt %in% 1:2)
label <- sprintf("%s-attempt-%d-start-%d", package, attempt, start)
method <- if (length(args) >= 5L) args[5] else "EVA"
if (package == "gllvm" && method != "EVA") label <- sub("gllvm-", paste0("gllvm-", method, "-"), label)
initialisation <- if (length(args) == 6L) args[6] else "random"
if (initialisation != "random") label <- paste0(label, "-", initialisation)
stopifnot(package == "sjSDM", attempt == 2L)
label <- sprintf("sjSDM-default-penalty-start-%d", start)
output <- file.path(root, "stability/fits", label)
stopifnot(!file.exists(paste0(output, ".rds")))
input <- readRDS(file.path(root, "inputs", "training.rds"))
input_hash <- unname(tools::md5sum(file.path(root, "inputs", "training.rds")))
warnings <- character()
started <- Sys.time()
settings <- list(package = package, attempt = attempt, start = start)

result <- tryCatch(withCallingHandlers({
  suppressPackageStartupMessages(library(package, character.only = TRUE))
  if (package %in% c("occJSDM", "Hmsc")) {
    stopifnot(start == 1L)
    settings$seed <- (if (attempt == 1) 26092200 else 26092300) +
      if (package == "occJSDM") 11 else 21
    settings$warmup <- if (attempt == 1) 2000L else 8000L
    settings$draws <- if (attempt == 1) 4000L else 16000L
    settings$chains <- 4L
    set.seed(settings$seed)
    if (package == "occJSDM") {
      fit <- runOccJSDM(list(info = input$x, OTU = input$y),
        occCovariates = names(input$x),
        listParams = list(n_factors = 2, n_lattrait = 0),
        MCMCparams = list(nchain = 4, nburn = settings$warmup,
                          niter = settings$draws, nthin = 1))
      stopifnot(fit$infos$model == "binary", fit$infos$ps == 0,
                fit$infos$n_factors == 2,
                dim(fit$results_output$jsdm_output$A_output)[2] == 0,
                identical(fit$infos$OTU, input$y),
                identical(fit$infos$speciesNames, colnames(input$y)),
                max(abs(fit$X_psi - as.matrix(input$x))) < 1e-12)
    } else {
      study_design <- data.frame(site = factor(rownames(input$x),
                                                levels = rownames(input$x)))
      random_level <- setPriors(HmscRandomLevel(units = levels(study_design$site)),
                               nfMin = 2, nfMax = 2)
      model <- Hmsc(Y = input$y, XData = input$x,
        XFormula = ~ environment_1 + environment_2, XScale = FALSE,
        distr = "probit", studyDesign = study_design, ranLevels = list(site = random_level))
      fit <- sampleMcmc(model, samples = settings$draws, transient = settings$warmup,
                       thin = 1, nChains = 4, nParallel = 1, verbose = 1000)
      stopifnot(identical(fit$Y, input$y),
                max(abs(fit$X - cbind(1, as.matrix(input$x)))) < 1e-12,
                max(abs(fit$X - fit$XScaled)) < 1e-12,
                all(vapply(fit$postList, function(chain) all(vapply(chain,
                  function(d) nrow(d$Lambda[[1]]) == 2, logical(1))), logical(1))))
    }
    saveRDS(fit, paste0(output, ".rds"), compress = FALSE)
  } else if (package == "gllvm") {
    settings$seed <- (if (attempt == 1) 26092230 else 26092330) + start
    settings$limit <- if (attempt == 1) 6000 else 12000
    settings$method <- method
    settings$initialisation <- initialisation
    fit <- gllvm(y = input$y, X = input$x, family = binomial("logit"),
      link = "logit", num.lv = 2, method = method, seed = settings$seed,
      sd.errors = FALSE,
      control.start = list(starting.val = initialisation, n.init = 1,
                           jitter.var = if (initialisation == "res") 0.1 else 0),
      control = list(maxit = settings$limit, max.iter = settings$limit))
    stopifnot(identical(fit$y, input$y), max(abs(fit$X - as.matrix(input$x))) < 1e-12,
              ncol(fit$lvs) == 2, fit$method == method, all(fit$link == "logit"))
    # Capture the gradient now, while the TMB external pointer is live.
    gradient <- tryCatch(fit$TMBfn$gr(fit$TMBfn$env$last.par.best), error = function(e) NA_real_)
    fit$pilot_gradient <- gradient
    saveRDS(fit, paste0(output, ".rds"), compress = FALSE)
  } else {
    stopifnot(Sys.getenv("SJSDM_MOJO_BACKEND") == "0")
    settings$seed <- (if (attempt == 1) 26092240 else 26092340) + start
    settings$epochs <- if (attempt == 1) 1000L else 3000L
    settings$sampling <- if (attempt == 1) 1000L else 2000L
    settings$learning_rate <- if (attempt == 1) 0.005 else 0.002
    fit <- sjSDM(Y = input$y, env = linear(data = input$x, lambda = 0),
      biotic = bioticStruct(df = 2, lambda = 0),
      family = binomial("logit"), device = "cpu", dtype = "float64",
      iter = settings$epochs, sampling = settings$sampling, step_size = 100L,
      learning_rate = settings$learning_rate, parallel = 0L,
      control = sjSDMControl(optimizer = RMSprop(weight_decay = 0.0001),
                            scheduler = 0, early_stopping_training = 0),
      seed = settings$seed, se = FALSE, verbose = FALSE)
    python <- reticulate::py
    python$pilot_model <- fit$model
    settings$weight_decay <- 0.0001
    settings$backend <- reticulate::py_eval("pilot_model._loss_function.__name__")
    settings$python <- reticulate::py_config()$python
    stopifnot(settings$backend == "torch_tmp",
      reticulate::py_eval("str(pilot_model.device)") == "cpu",
      reticulate::py_eval("int(pilot_model.df)") == 2L,
      identical(fit$data$Y, input$y),
      max(abs(fit$data$X - cbind(1, as.matrix(input$x)))) < 1e-12)
    # Keep a full numeric fit. Python external pointers cannot survive RDS reload.
    native_raw <- predict(fit, newdata = input$x, type = "raw")
    fit$model <- NULL
    fit$get_model <- NULL
    fit$pilot_native_raw <- native_raw
    fit$settings$control$optimizer$ff <- NULL
    saveRDS(fit, paste0(output, ".rds"), compress = FALSE)
  }
  list(ok = TRUE)
}, warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
}), error = function(e) list(ok = FALSE, error = conditionMessage(e)))

record <- c(result, list(settings = settings, input_hash = input_hash,
                        warnings = unique(warnings), started = started,
                        ended = Sys.time(), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
                        session = sessionInfo()))
saveRDS(record, paste0(output, "-record.rds"))
print(record[setdiff(names(record), "session")])
if (!isTRUE(record$ok)) quit(status = 1)
