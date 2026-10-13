# Freeze the study setup: data, folds, jobs, seeds and MCMC schedule.
# Run once from the repository root, with the study library first in R_LIBS:
# Rscript dev/simstudy/waic-cross-validation/prepare.R OUTDIR
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
outdir <- normalizePath(args[1], mustWork = TRUE)
setup_path <- file.path(outdir, "setup.rds")
if (file.exists(setup_path)) stop("setup.rds already exists; a study setup is never overwritten.")
source("dev/simstudy/waic-cross-validation/cv-lib.R")
stopifnot(identical(normalizePath(find.package("occJSDM")),
                    normalizePath(file.path(outdir, "library", "occJSDM"))))
md5 <- function(path) unname(tools::md5sum(path))
lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
lesson <- readRDS(lesson_path)
data <- lesson$input$sim$data_list
sites <- sort(unique(data$info$Site))
stopifnot(length(sites) == 100L, nrow(data$OTU) == 3600L, ncol(data$OTU) == 10L)
archive_input <- "/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-post-pr14-20261007/input.rds"
folds <- cv_make_folds(sites)
jobs <- do.call(rbind, lapply(cv_job_names(), function(job) as.data.frame(cv_parse_job(job))))
scripts <- list.files("dev/simstudy/waic-cross-validation", pattern = "\\.(R|sh)$", full.names = TRUE)
setup <- list(
  schema = 1L, created = Sys.time(), data = data, lesson_md5 = md5(lesson_path),
  archive_input_identical = if (file.exists(archive_input)) {
    identical(readRDS(archive_input), lesson$input)
  } else NA,
  folds = folds, jobs = jobs, mcmc = cv_mcmc, rng_kind = RNGkind(),
  package_commit = readLines(file.path(outdir, "source-commit.txt")),
  package_version = as.character(utils::packageVersion("occJSDM")),
  R = R.version.string, script_md5 = tools::md5sum(scripts))
saveRDS(setup, setup_path)
write.csv(folds, file.path(outdir, "folds.csv"), row.names = FALSE)
write.csv(jobs, file.path(outdir, "jobs.csv"), row.names = FALSE)
writeLines(jobs$job, file.path(outdir, "jobs.txt"))
cat("Saved setup for", nrow(jobs), "jobs; archived input identical:",
    setup$archive_input_identical, "\n")
print(table(folds$fold))
