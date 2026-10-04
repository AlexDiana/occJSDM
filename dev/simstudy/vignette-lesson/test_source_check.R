# Tests for lesson_source_check(); run from the repository root:
# Rscript dev/simstudy/vignette-lesson/test_source_check.R
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/source-check.R")

stopifnot(exists("lesson_source_check", mode = "function"),
          is.character(lesson_plot_only_changes),
          length(lesson_plot_only_changes) > 0L)

repo <- tempfile("source-check-")
dir.create(file.path(repo, "R"), recursive = TRUE)
dir.create(file.path(repo, "src"))
old_wd <- setwd(repo)

git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop("git ", paste(c(...), collapse = " "), " failed:\n",
                                             paste(out, collapse = "\n"))
  invisible(out)
}
put <- function(path, lines) writeLines(lines, path)
commit <- function(message) {
  git("add", "-A")
  git("-c", "user.name=test", "-c", "user.email=test@example.org",
      "commit", "-q", "-m", shQuote(message))
}
# Expect lesson_source_check() to fail with a message matching every pattern.
fails <- function(recorded, ..., allow = c("plotA")) {
  message <- tryCatch({ lesson_source_check(recorded, allow = allow); NULL },
                      error = function(e) conditionMessage(e))
  if (is.null(message)) stop("lesson_source_check() passed but should have failed")
  for (pattern in c(...)) if (!grepl(pattern, message, fixed = TRUE))
    stop("failure message lacks '", pattern, "':\n", message)
  invisible(message)
}
passes <- function(recorded, allow = c("plotA"))
  isTRUE(lesson_source_check(recorded, allow = allow))

git("init", "-q")
model <- c("#' Fit the model",
           "#' @export",
           "fitModel <- function(x, y) {",
           "  x + y",
           "}",
           "",
           "plotA = function(fit) {",
           "  scale <- c(\"red\", \"blue\")",
           "  scale",
           "}",
           "",
           "NULL")
put("R/model.R", model)
put("R/other.R", c("helper <- function(x) x * 2"))
put("DESCRIPTION", c("Package: toy", "Version: 0.0.1"))
put("NAMESPACE", c("export(fitModel)"))
put("src/core.cpp", c("// core", "int f() { return 1; }"))
commit("first")
files <- c("DESCRIPTION", "NAMESPACE", "R/model.R", "R/other.R", "src/core.cpp")
recorded <- tools::md5sum(files)

# A later commit, so the recorded version is not HEAD and must be found in history.
put("R/other.R", c("# A comment-only commit", "helper <- function(x) x * 2"))
commit("second")

# 1. Nothing changed in R/model.R since recording; R/other.R changed only in comments.
stopifnot(passes(recorded))

# 2. Comment-only edits: roxygen, plain comment lines, and a trailing comment.
put("R/model.R", c("#' Fit the model, now documented at length",
                   "#' @param x a number",
                   "#' @export",
                   "# A plain comment line",
                   "fitModel <- function(x, y) {",
                   "  x + y  # a trailing comment",
                   "}",
                   "",
                   "plotA = function(fit) {",
                   "  scale <- c(\"red\", \"blue\")",
                   "  scale",
                   "}",
                   "",
                   "NULL"))
stopifnot(passes(recorded))
commit("comments")
stopifnot(passes(recorded))

# 3. A code change in a function not on the allow list fails and names it.
put("R/model.R", sub("x + y", "x - y", model, fixed = TRUE))
fails(recorded, "R/model.R", "fitModel")
# ... even when that function is the only change and is committed.
commit("change fitModel")
fails(recorded, "R/model.R", "fitModel")
put("R/model.R", model); commit("restore fitModel")
stopifnot(passes(recorded))

# 4. A code change in an allow-listed function passes; the same change fails
#    when the allow list is empty.
put("R/model.R", sub("\"red\", \"blue\"", "\"#E69F00\", \"#56B4E9\"", model, fixed = TRUE))
stopifnot(passes(recorded))
fails(recorded, "R/model.R", "plotA", allow = character())

# 5. A changed non-function top-level expression fails.
put("R/model.R", c(model, "options(toy.flag = TRUE)"))
fails(recorded, "R/model.R", "top-level")
put("R/model.R", model)

# 6. A removed or added function fails unless allow-listed.
put("R/model.R", model[1:5])
fails(recorded, "R/model.R", "plotA", allow = character())
put("R/model.R", c(model, "extra <- function() 1"))
fails(recorded, "R/model.R", "extra")
put("R/model.R", model)
stopifnot(passes(recorded))

# 7. A deleted recorded file fails.
invisible(file.remove("R/other.R"))
fails(recorded, "R/other.R")
put("R/other.R", c("# A comment-only commit", "helper <- function(x) x * 2"))
stopifnot(passes(recorded))

# 8. New R/ files: only allow-listed functions may appear.
put("R/colours.R", c("#' Colours", "plotA <- function() NULL"))
stopifnot(passes(recorded))
put("R/colours.R", c("#' Colours", "plotA <- function() NULL", "newFit <- function() 2"))
fails(recorded, "R/colours.R", "newFit")
put("R/colours.R", c("palette_size <- 8L"))
fails(recorded, "R/colours.R", "top-level")
invisible(file.remove("R/colours.R"))
stopifnot(passes(recorded))

# 9. Non-R files are compared by raw md5, so even a comment there fails.
put("DESCRIPTION", c("Package: toy", "Version: 0.0.1", ""))
fails(recorded, "DESCRIPTION")
put("DESCRIPTION", c("Package: toy", "Version: 0.0.1"))
put("NAMESPACE", c("# comment", "export(fitModel)"))
fails(recorded, "NAMESPACE")
put("NAMESPACE", c("export(fitModel)"))
put("src/core.cpp", c("// core, edited", "int f() { return 1; }"))
fails(recorded, "src/core.cpp")
put("src/core.cpp", c("// core", "int f() { return 1; }"))
stopifnot(passes(recorded))

# 10. A recorded md5 that matches no version in history fails.
bogus <- recorded; bogus[["R/model.R"]] <- "00000000000000000000000000000000"
fails(bogus, "R/model.R", "no commit")

# 11. Several failures are all named in one message.
put("R/model.R", sub("x + y", "x - y", model, fixed = TRUE))
invisible(file.remove("R/other.R"))
fails(recorded, "R/model.R", "fitModel", "R/other.R")

setwd(old_wd); unlink(repo, recursive = TRUE)
cat("All lesson_source_check() tests passed.\n")
