# Fit one study job and save it. Never overwrites a saved fit.
# Rscript OUTDIR/scripts/fit-job.R OUTDIR JOB   (JOB as in OUTDIR/jobs.txt)
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
outdir <- normalizePath(args[1], mustWork = TRUE)
job <- args[2]
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "cv-lib.R"))
stopifnot(identical(normalizePath(find.package("occJSDM")),
                    normalizePath(file.path(outdir, "library", "occJSDM"))),
          Sys.getenv("RCPP_PARALLEL_NUM_THREADS") == "1")
setup_path <- file.path(outdir, "setup.rds")
setup <- readRDS(setup_path)
stopifnot(job %in% setup$jobs$job)
spec <- cv_parse_job(job)
dir.create(file.path(outdir, "fits"), showWarnings = FALSE)
path <- file.path(outdir, "fits", paste0(job, ".rds"))
if (file.exists(path)) stop("Fit already saved: ", path)
sites <- sort(unique(setup$data$info$Site))
heldout <- if (is.na(spec$fold)) sites[0] else setup$folds$Site[setup$folds$fold == spec$fold]
result <- cv_fit(cv_subset(setup$data, setdiff(sites, heldout)), spec$model, spec$seed, setup$mcmc)
result$job <- job
result$fold <- spec$fold
result$heldout_sites <- heldout
result$session <- sessionInfo()
result$threads <- Sys.getenv(c("RCPP_PARALLEL_NUM_THREADS", "OMP_NUM_THREADS",
                               "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
result$setup_md5 <- unname(tools::md5sum(setup_path))
partial <- paste0(path, ".partial")
saveRDS(result, partial, compress = "gzip")
stopifnot(file.rename(partial, path))
cat("Saved", job, "in", round(as.numeric(difftime(result$finished, result$started, units = "secs")), 1),
    "seconds with", length(result$warnings), "warnings\n")
