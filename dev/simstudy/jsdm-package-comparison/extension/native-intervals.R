# Frozen-environment worker. Reads fits and training inputs, never truth/.
# run-r native-intervals.R ARCHIVE CODE OUTPUT [comma-separated job numbers]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
root <- normalizePath(args[1]); code <- normalizePath(args[2]); out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(dirname(code), "pilot-math.R"))
source(file.path(code, "calibration-math.R"))
source(file.path(code, "native-interval-math.R"))
stopifnot(as.character(packageVersion("gllvm")) == "2.0.15",
          as.character(packageVersion("sjSDM")) == "1.0.7",
          Sys.getenv("SJSDM_MOJO_BACKEND") == "0")
manifest <- read.csv(file.path(root, "status.csv"))
jobs <- subset(manifest, scored & package %in% c("gllvm", "sjSDM") & scenario != "traits" &
                 !(scenario == "curved" & response == "linear"))
if (length(args) >= 4L) jobs <- jobs[jobs$job_number %in% as.integer(strsplit(args[4], ",")[[1]]), ]
sources <- file.path(code, c("native-intervals.R", "native-interval-math.R", "native-sjsdm.py", "calibration-math.R"))
sources <- c(sources, file.path(dirname(code), "pilot-math.R"))
source_hashes <- tools::md5sum(sources)
dir.create(file.path(out, "jobs"), showWarnings = FALSE)
started <- Sys.time()
python_loaded <- FALSE
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i, ]
  result_path <- file.path(root, "jobs", j$job, "result.rds")
  input_path <- file.path(root, "inputs", paste0(j$input, ".rds"))
  r <- readRDS(result_path); input <- readRDS(input_path)
  input$test_x <- NULL
  start_path <- file.path(dirname(result_path), paste0("start-", r$selected_start, ".rds"))
  input_hashes <- tools::md5sum(c(result_path, input_path, start_path))
  stopifnot(isTRUE(r$ok), !r$truth_used, unname(input_hashes[2]) == j$input_md5)
  destination <- file.path(out, "jobs", paste0(j$job, ".rds"))
  if (file.exists(destination)) {
    old <- readRDS(destination)
    stopifnot(identical(old$source_hashes, source_hashes), identical(old$input_hashes, input_hashes))
    next
  }
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    k <- nrow(r$parameters$beta); species <- ncol(input$y)
    if (j$package == "gllvm") {
      # Restore the original selected start only, with the original seed/options.
      # The VA nuisance parameters/TMB tape were not archived, so replay is needed.
      fit <- gllvm::gllvm(y = input$y, X = input$x, family = binomial("logit"), link = "logit",
        num.lv = 2, method = "VA", seed = 26100000L + 100L*j$job_number + r$selected_start,
        sd.errors = TRUE, control.start = list(starting.val = "res", n.init = 1, jitter.var = .1),
        control = list(maxit = 12000, max.iter = 12000))
      p <- point_parameters(fit, "gllvm")
      parameter_change <- max(abs(c(p$beta - r$parameters$beta, p$loading - r$parameters$loading)))
      saved <- readRDS(start_path)
      loglik_change <- abs(fit$logL - saved$native_loglik)
      stopifnot(parameter_change <= 1e-8, loglik_change <= 1e-8)
      v <- vcov(fit)
      idx <- which(colnames(v) == "b")
      stopifnot(length(idx) == k*species,
        max(abs(fit$TMBfn$par[names(fit$TMBfn$par) == "b"] - as.vector(p$beta))) <= 1e-8)
      covariance <- lapply(seq_len(species), function(s) v[idx[(s-1L)*k + seq_len(k)], idx[(s-1L)*k + seq_len(k)], drop=FALSE])
      native_se <- rbind(fit$sd$beta0, t(fit$sd$Xcoef))
      for (s in seq_len(species)) stopifnot(max(abs(pmax(diag(covariance[[s]]), 0) - native_se[,s]^2)) < 1e-8)
      audit <- list(parameter_change = parameter_change, loglik_change = loglik_change,
        native_se = native_se, covariance = covariance,
        joint_min_eigenvalue = min(eigen(v, symmetric=TRUE, only.values=TRUE)$values),
        standard_error_method = deparse(get("se.gllvm", asNamespace("gllvm"))))
      method <- "gllvm native VA covariance; selected-start replay"
      numerical_pass <- audit$joint_min_eigenvalue > 0
    } else {
      if (!python_loaded) {
        suppressPackageStartupMessages(library(sjSDM))
        sjSDM:::check_module()
        reticulate::source_python(file.path(code, "native-sjsdm.py"))
        python_loaded <- TRUE
      }
      model <- sjSDM:::pkg.env$fa$Model_sjSDM(device="cpu", dtype="float64", seed=19L)
      model$add_env(as.integer(k), as.integer(species), intercept=TRUE, l1=0, l2=0)
      model$build(df=2L, link="logit", alpha=1, scheduler=FALSE)
      model$set_env_weights(list(t(r$parameters$beta)))
      model$set_sigma(t(r$parameters$loading))
      stopifnot(max(abs(sjSDM:::force_r(model$env_weights)[[1]] - t(r$parameters$beta))) == 0,
        max(abs(sjSDM:::force_r(model$get_sigma) - t(r$parameters$loading))) == 0)
      x <- as.matrix(cbind(1, input$x))
      seed <- 26200000L + 100L*j$job_number
      checks <- list()
      # The native routine defaults to only 100 integration draws. Increase its
      # supported sampling argument and audit with an independent random stream.
      # This rule uses numerical stability, never truth or achieved coverage.
      for (sampling in c(2000L, 8000L, 32000L)) {
        a <- native_covariance(model, x, input$y, sampling, seed)
        b <- native_covariance(model, x, input$y, sampling, seed+1L)
        covariance <- lapply(seq_len(species), function(s) a$covariance[s,,])
        other <- lapply(seq_len(species), function(s) b$covariance[s,,])
        primary <- native_coefficient_intervals(r$parameters$beta, covariance, input$centre, input$spread)
        secondary <- native_coefficient_intervals(r$parameters$beta, other, input$centre, input$spread)
        change <- if (anyNA(c(primary$se, secondary$se))) Inf else
          max(abs(primary$se-secondary$se) / pmax(primary$se, secondary$se))
        checks[[length(checks)+1L]] <- data.frame(sampling=sampling, max_relative_se_change=change)
        if (change <= .05) break
      }
      audit <- list(covariance=covariance, native_se=t(a$se),
        regularized_hessian=a$regularized_hessian, secondary_covariance=other,
        checks=do.call(rbind, checks), sampling=sampling, seeds=c(seed, seed+1L),
        parameter_change=0, loglik_change=NA_real_,
        standard_error_source=tools::md5sum(system.file("python/sjSDM_py/model_sjSDM.py", package="sjSDM")))
      method <- "sjSDM native conditional Hessian; frozen associations and other species"
      numerical_pass <- change <= .05
    }
    intervals <- native_coefficient_intervals(r$parameters$beta, covariance, input$centre, input$spread)
    list(ok=TRUE, method=method, intervals=intervals, audit=audit,
      numerical_pass=numerical_pass && all(intervals$covariance_valid))
  }, warning=function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error=function(e) list(ok=FALSE, numerical_pass=FALSE, error=conditionMessage(e)))
  result <- c(result, list(job=j, warnings=unique(warnings), input_hashes=input_hashes,
    source_hashes=source_hashes, truth_used=FALSE, session=sessionInfo(), finished=Sys.time()))
  saveRDS(result, paste0(destination, ".partial"))
  stopifnot(file.rename(paste0(destination, ".partial"), destination))
  cat(sprintf("%d/%d %s ok=%s numerical_pass=%s elapsed=%.1fs\n", i, nrow(jobs), j$job,
    result$ok, result$numerical_pass, as.numeric(difftime(Sys.time(), started, units="secs"))))
  if (!is.null(result$error)) cat(result$error, "\n")
}
