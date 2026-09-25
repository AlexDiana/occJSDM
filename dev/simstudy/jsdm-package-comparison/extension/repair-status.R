# Add verified JSON summaries for completed implementation-check workers.
root <- normalizePath(commandArgs(trailingOnly = TRUE)[1])
manifest <- read.csv(file.path(root, "manifest.csv"))
for (i in seq_len(nrow(manifest))) {
  out <- file.path(root, "jobs", manifest$job[i])
  path <- file.path(out, "result.rds")
  if (!file.exists(path) || file.exists(file.path(out, "status.json"))) next
  r <- readRDS(path)
  stopifnot(r$job$job_number == manifest$job_number[i], r$input_md5 == manifest$input_md5[i], !r$truth_used)
  s <- r[intersect(names(r), c("ok", "diagnostic_pass", "error", "selected_start", "elapsed_seconds", "warnings", "input_md5", "truth_used"))]
  s$job_number <- manifest$job_number[i]
  s$result_md5 <- unname(tools::md5sum(path))
  jsonlite::write_json(s, file.path(out, "status.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
}
