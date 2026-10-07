# Cost and outcome record for the twelve native weak-penalty sjSDM starts, in
# the same columns as pilot-results/all-fit-attempts.csv so the lesson's
# runtime table can include them. Starts 1:3 combine the original 3,000-epoch
# fit and its 1,000-epoch continuation; starts 4:12 ran both stages in one
# process. Objectives are precise penalised training scores from the checker.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
check <- read.csv(file.path(root, "stability-resolution/multistart/multistart-summary.csv"))
rows <- lapply(1:12, function(start) {
  if (start <= 3L) {
    a <- readRDS(file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-start-%d-record.rds", start)))
    b <- readRDS(file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-polish-start-%d-record.rds", start)))
    seconds <- a$elapsed_seconds + as.numeric(difftime(b$ended, b$started, units = "secs"))
    warnings <- paste(unique(c(a$warnings, b$warnings)), collapse = "; ")
    fit <- sprintf("sjSDM-default-penalty-polish-start-%d", start)
    note <- "Original weak-penalty start plus separate low-step continuation; both completed"
  } else {
    r <- readRDS(file.path(root, "stability-resolution/multistart", sprintf("sjSDM-multistart-start-%d-record.rds", start)))
    stopifnot(isTRUE(r$ok))
    seconds <- r$elapsed_seconds
    warnings <- paste(r$warnings, collapse = "; ")
    fit <- sprintf("sjSDM-multistart-start-%d", start)
    note <- "Fresh weak-penalty start with low-step continuation in one process; completed"
  }
  s <- check[check$start == start, ]
  f <- if (start <= 3L) readRDS(file.path(root, "stability/fits", sprintf("sjSDM-default-penalty-polish-start-%d.rds", start))) else
    readRDS(file.path(root, "stability-resolution/multistart", sprintf("sjSDM-multistart-start-%d.rds", start)))
  data.frame(fit = fit, package = "sjSDM", seconds = seconds,
    process_note = paste0(note, "; local maximum ", s$basin),
    rhat = NA_real_, bulk_ess = NA_real_, tail_ess = NA_real_,
    objective = s$native_penalised_score, maximum_gradient = s$native_max_gradient,
    maximum_coefficient = max(abs(f$parameters$beta)),
    maximum_residual_sd = s$maximum_residual_sd, warnings = warnings)
})
attempts <- do.call(rbind, rows)
template <- read.csv(file.path(code, "pilot-results/all-fit-attempts.csv"))
stopifnot(identical(names(attempts), names(template)))
write.csv(attempts, file.path(code, "stability-resolution-results/sjsdm-weak-penalty-attempts.csv"), row.names = FALSE)
print(attempts[, c("fit", "seconds", "objective", "process_note")])
