# Shared constants and helpers of the Hmsc longer-runs check (README.md in
# this directory). Sourced by select-jobs.R, refit.R, launch.R, compare.R and
# test-refit.R; never run on its own.

`%||%` <- function(x,y) if(is.null(x)) y else x

HERE <- local({
  here <- NULL
  for(i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if(!is.null(f)) {here <- dirname(normalizePath(f));break}
  }
  if(is.null(here)) {
    a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
    if(length(a)) here <- dirname(normalizePath(sub('^--file=','',a[1])))
  }
  here %||% getwd()
})
SCRIPT_REPO <- normalizePath(file.path(HERE,'../../..'))

STUDY_REL <- 'dev/simstudy/hmsc-longer-runs'
JOBS_FILE <- file.path(STUDY_REL,'results/jobs.csv')
# The Lesson 6 extension's analysis code, as committed, at the md5s the
# calibration bundle records for them (lesson-4-calibration.rds, provenance).
CODE_REL <- 'dev/simstudy/jsdm-package-comparison'
CODE_MD5 <- c('pilot-math.R'='89f36f29919ff0c51915fb6966c6d743',
  'extension/math.R'='793f7bad7ea7bb93df09eada0e0302ad',
  'extension/study.R'='5d3a5a6cdef9508bd29e3594d57052f1',
  'extension/calibration-math.R'='6a522e2450d5458f3152ac029d2e71f1')
# The frozen worker whose Hmsc statements fit_hmsc() copies verbatim
# (worker-final/extension/fit-job.R in the extension archive).
FIT_JOB_REL <- 'worker-final/extension/fit-job.R'
FIT_JOB_MD5 <- '67b62869071142a9426fe4467f342d53'
HMSC_VERSION <- '3.3.7'
# The original longer attempt (attempt 2) and this check's schedule (attempt 3):
# four times the iterations, the same number of retained draws per chain.
ORIGINAL_LONGER <- list(warmup=8000L,draws=16000L,thin=1L,chains=4L)
SCHEDULE <- list(warmup=32000L,draws=16000L,thin=4L,chains=4L)
ATTEMPT <- 3L
BUNDLE_REL <- 'vignettes/teaching-data/lesson-4-calibration.rds'
RESUMED_TEXT <- 'already complete with identical settings; resumed'

# The frozen worker's seed rule, base_seed + attempt, with attempt 3.
job_seed <- function(job_number) 26100000L+100L*as.integer(job_number)+ATTEMPT

iterations_per_chain <- function(s) s$warmup+s$draws*s$thin

parse_options <- function(args,known,required=character()) {
  bad <- !grepl('^--[A-Za-z][A-Za-z-]*=',args)
  if(any(bad)) stop('Malformed argument: ',args[bad][1])
  opt <- sub('^--([^=]+)=.*$','\\1',args);value <- sub('^--[^=]+=','',args)
  unknown <- setdiff(opt,known)
  if(length(unknown)) stop('Unknown option: --',unknown[1])
  if(anyDuplicated(opt)) stop('Repeated option: --',opt[duplicated(opt)][1])
  missing_opt <- setdiff(required,opt)
  if(length(missing_opt)) stop('Missing --',missing_opt[1])
  as.list(stats::setNames(value,opt))
}

hash_files <- function(files,root,label='') {
  root <- normalizePath(root,mustWork=TRUE);full <- normalizePath(files,mustWork=TRUE)
  prefix <- paste0(root,'/');inside <- startsWith(full,prefix)
  if(!all(inside)) stop('File outside ',root,': ',full[!inside][1])
  stats::setNames(unname(tools::md5sum(full)),paste0(label,substring(full,nchar(prefix)+1L)))
}

# The extension's analysis code, md5-checked, in one environment.
load_code <- function(repo) {
  files <- file.path(repo,CODE_REL,names(CODE_MD5))
  found <- unname(tools::md5sum(files))
  bad <- names(CODE_MD5)[found!=CODE_MD5]
  if(length(bad)) stop('Extension code changed since the calibration bundle recorded it: ',paste(bad,collapse=', '))
  code <- new.env(parent=globalenv())
  for(f in files) sys.source(f,envir=code)
  code
}

check_fit_job <- function(root) {
  f <- file.path(root,FIT_JOB_REL)
  found <- unname(tools::md5sum(f))
  if(!identical(found,FIT_JOB_MD5)) stop('Frozen worker ',f,' has md5 ',found,', not ',FIT_JOB_MD5)
  invisible(f)
}

load_hmsc <- function(study) {
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  .libPaths(c(file.path(study,'library'),.libPaths()))
  suppressPackageStartupMessages(library(Hmsc))
  pkg <- normalizePath(find.package('Hmsc'))
  if(!identical(pkg,normalizePath(file.path(study,'library/Hmsc')))) stop('Hmsc loaded from ',pkg,', not the study library')
  v <- as.character(utils::packageVersion('Hmsc'))
  if(!identical(v,HMSC_VERSION)) stop('Hmsc version ',v,' is not ',HMSC_VERSION)
  pkg
}

read_jobs <- function(repo) {
  f <- file.path(repo,JOBS_FILE)
  if(!file.exists(f)) stop('Missing job selection: ',f,'; run select-jobs.R first')
  utils::read.csv(f,stringsAsFactors=FALSE)
}

extension_manifest <- function(root) utils::read.csv(file.path(root,'manifest.csv'),stringsAsFactors=FALSE)

job_input <- function(root,job) {
  path <- file.path(root,'inputs',paste0(job$input,'.rds'))
  found <- unname(tools::md5sum(path))
  if(!identical(found,job$input_md5)) stop('Input md5 mismatch for ',path)
  readRDS(path)
}

# The Hmsc model and sampler exactly as the frozen worker built them (fit-job.R,
# run_bayesian(), the Hmsc branch), with this check's schedule and seed. Test
# predictors are dropped before fitting, as the worker drops them.
fit_hmsc <- function(input,schedule,seed) {
  input$test_x <- NULL
  formula_env <- reformulate(names(input$x))
  set.seed(seed)
  design <- data.frame(site = factor(rownames(input$x), levels = rownames(input$x)))
  random <- Hmsc::setPriors(Hmsc::HmscRandomLevel(units = levels(design$site)), nfMin = 2, nfMax = 2)
  model <- Hmsc::Hmsc(Y = input$y, XData = input$x, XFormula = formula_env, XScale = FALSE,
    TrData = input$traits, TrFormula = if (is.null(input$traits)) NULL else ~ drought_tolerance + irrelevant_trait,
    distr = "probit", studyDesign = design, ranLevels = list(site = random))
  fit <- Hmsc::sampleMcmc(model, samples = schedule$draws, transient = schedule$warmup, thin = schedule$thin,
                          nChains = schedule$chains, nParallel = 1, verbose = 1000)
  stopifnot(identical(fit$Y, input$y), max(abs(fit$X - cbind(1, as.matrix(input$x)))) < 1e-10)
  fit
}

study_job_dir <- function(study,job_name) file.path(study,'jobs',job_name)
