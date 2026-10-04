# lesson_source_check() for the lesson verifiers. Source after helpers.R (it uses
# lesson_source_hashes()). Kept out of helpers.R because bundles record the md5
# of helpers.R and of the producer scripts, so editing those would break them.

# Package functions whose code may differ from the version a fit was made with.
# Added 4 October 2026 for PR #23 (Okabe-Ito colours): the changes are colour
# scales only. These functions draw plots from a finished fit and are never
# called while fitting, so a change to them cannot change any fit.
lesson_plot_only_changes <- c("plotFPTPStage2Rates", "plotDetectionRates",
                              "plotStage2FPRates", "plotTraceplot", "okabe_ito")

# Compare the source fingerprint recorded with a fit (lesson_source_hashes())
# against the current tree, ignoring comments. Run from the repository root.
# Files with equal raw md5 pass. A changed R/*.R file is looked up in git
# history by its recorded md5, both versions are parsed without source
# references (which drops comments) and compared function by function: only
# functions named in `allow` may be added, removed or changed, and every other
# top-level expression must deparse identically. New R/*.R files may define
# only allow-listed functions. Every other file (DESCRIPTION, NAMESPACE, src/)
# must match its raw md5. Returns TRUE invisibly or stops naming every
# difference.
lesson_source_check <- function(recorded, allow = lesson_plot_only_changes) {
  stopifnot(is.character(recorded), length(recorded) > 0L,
            !is.null(names(recorded)), is.character(allow))
  is_r <- function(path) grepl("^R/.*\\.[Rr]$", path)
  shown <- new.env(parent = emptyenv())
  # The file at a commit, written to a temporary file and cached within this
  # call; NULL if the commit does not have it.
  show <- function(commit, path) {
    key <- paste0(commit, ":", path)
    if (is.null(shown[[key]])) {
      tmp <- tempfile(fileext = ".R")
      rc <- suppressWarnings(system2("git", c("show", "--no-textconv",
                                              shQuote(paste0(commit, ":./", path))),
                                     stdout = tmp, stderr = FALSE))
      shown[[key]] <- if (identical(as.integer(rc), 0L)) tmp else NA_character_
    }
    if (is.na(shown[[key]])) NULL else shown[[key]]
  }
  recorded_version <- function(path, md5) {
    commits <- suppressWarnings(system2("git", c("log", "--format=%H", "--", shQuote(path)),
                                        stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(commits, "status"))) return(NULL)
    for (commit in commits) {
      file <- show(commit, path)
      if (!is.null(file) && identical(unname(tools::md5sum(file)), unname(md5)))
        return(file)
    }
    NULL
  }
  # Top-level functions (name -> deparsed definitions) and other expressions.
  definitions <- function(file) {
    exprs <- parse(file, keep.source = FALSE)
    text <- function(e) paste(deparse(e, width.cutoff = 500L), collapse = "\n")
    functions <- list(); others <- character()
    for (e in as.list(exprs)) {
      if (is.call(e) && length(e) == 3L &&
          (identical(e[[1L]], as.name("<-")) || identical(e[[1L]], as.name("="))) &&
          (is.name(e[[2L]]) || is.character(e[[2L]])) &&
          is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
        name <- as.character(e[[2L]])
        functions[[name]] <- c(functions[[name]], text(e[[3L]]))
      } else others <- c(others, text(e))
    }
    list(functions = functions, others = others)
  }
  problems <- character()
  note <- function(path, what) problems <<- c(problems, paste0(path, ": ", what))
  for (path in names(recorded)) {
    if (!file.exists(path)) { note(path, "recorded file is missing"); next }
    if (identical(unname(tools::md5sum(path)), unname(recorded[[path]]))) next
    if (!is_r(path)) { note(path, "md5 differs (compared byte for byte)"); next }
    old_file <- recorded_version(path, recorded[[path]])
    if (is.null(old_file)) {
      note(path, "no commit in git history has the recorded md5"); next
    }
    old <- definitions(old_file); new <- definitions(path)
    for (name in setdiff(union(names(old$functions), names(new$functions)), allow)) {
      if (!identical(old$functions[[name]], new$functions[[name]]))
        note(path, paste0("function ", name, " ",
                          if (is.null(old$functions[[name]])) "added"
                          else if (is.null(new$functions[[name]])) "removed"
                          else "changed"))
    }
    if (!identical(old$others, new$others))
      note(path, "top-level code other than function definitions changed")
  }
  for (path in setdiff(names(lesson_source_hashes()), names(recorded))) {
    if (!is_r(path)) { note(path, "new file not in the recorded fingerprint"); next }
    new <- definitions(path)
    for (name in setdiff(names(new$functions), allow))
      note(path, paste0("new file defines function ", name))
    if (length(new$others))
      note(path, "new file has top-level code other than function definitions")
  }
  if (length(problems))
    stop("Package source differs from the recorded fingerprint beyond comments ",
         "and allowed plot-only functions:\n  ", paste(problems, collapse = "\n  "),
         call. = FALSE)
  invisible(TRUE)
}
