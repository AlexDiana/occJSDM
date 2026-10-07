# lesson_source_check() for the lesson verifiers. Source after helpers.R (it uses
# lesson_source_hashes()). Kept out of helpers.R because bundles record the md5
# of helpers.R and of the producer scripts, so editing those would break them.

# Package functions whose code may differ from the version a fit was made with.
# Added 4 October 2026 for PR #23 (Okabe-Ito colours): the changes are colour
# scales only. These functions draw plots from a finished fit and are never
# called while fitting, so a change to them cannot change any fit.
# lesson_source_check() enforces "never called while fitting": it fails if any
# top-level code in R/ other than these definitions names one of them.
# Every entry needs a dated justification like this one. Entries do not expire:
# bundles keep the fit-time hashes, so a later change to a listed function also
# passes this check. What guards those later changes is the figure verifiers'
# comparison of the installed library against R/ (body(), formals() and
# package_files_md5).
lesson_plot_only_changes <- c("plotFPTPStage2Rates", "plotDetectionRates",
                              "plotStage2FPRates", "plotTraceplot", "okabe_ito")

# Compare the source fingerprint recorded with a fit (lesson_source_hashes())
# against the current tree, ignoring comments. Run from the repository root.
# Files with equal raw md5 pass. A changed R/*.R file is looked up in git
# history by its recorded md5, and both versions are parsed without source
# references, which drops comments. Top-level definitions of functions named in
# `allow` are removed from each version, and the remaining ordered sequences of
# top-level expressions must be identical() as language objects (so constants
# are compared exactly, 0.1 vs 0.1000000000000001 and 1L vs 1 included). New
# R/*.R files may contain only allow-listed function definitions. Every other
# file (DESCRIPTION, NAMESPACE, src/) must match its raw md5. In every current
# R/*.R file, top-level code other than the allow-listed definitions must not
# name an allow-listed function, as a symbol or as a string. Returns TRUE
# invisibly or stops naming every difference.
# The git lookup searches only history reachable from HEAD (git log). If the
# commit a fit was made at is not reachable (a rebased or squash-merged branch,
# a shallow clone), the check fails closed with "no commit in git history has
# the recorded md5", which is not evidence of a code change.
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
  # Ordered top-level expressions without allow-listed function definitions,
  # plus the remaining function definitions by name and the other expressions,
  # which are used only to name what differs.
  definitions <- function(file) {
    exprs <- as.list(parse(file, keep.source = FALSE))
    defined <- vapply(exprs, function(e) {
      if (is.call(e) && length(e) == 3L &&
          (identical(e[[1L]], as.name("<-")) || identical(e[[1L]], as.name("="))) &&
          (is.name(e[[2L]]) || is.character(e[[2L]])) &&
          is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function")))
        as.character(e[[2L]]) else NA_character_
    }, character(1))
    kept <- is.na(defined) | !(defined %in% allow)
    sequence <- exprs[kept]; names_kept <- defined[kept]
    is_fn <- !is.na(names_kept)
    list(sequence = sequence,
         functions = split(sequence[is_fn], factor(names_kept[is_fn])),
         others = sequence[!is_fn])
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
    if (identical(old$sequence, new$sequence)) next
    before <- length(problems)
    for (name in union(names(old$functions), names(new$functions))) {
      if (!identical(old$functions[[name]], new$functions[[name]]))
        note(path, paste0("function ", name, " ",
                          if (is.null(old$functions[[name]])) "added"
                          else if (is.null(new$functions[[name]])) "removed"
                          else "changed"))
    }
    if (!identical(old$others, new$others))
      note(path, "top-level code other than function definitions changed")
    if (length(problems) == before)
      note(path, "order of top-level code changed")
  }
  for (path in setdiff(names(lesson_source_hashes()), names(recorded))) {
    if (!is_r(path)) { note(path, "new file not in the recorded fingerprint"); next }
    new <- definitions(path)
    for (name in names(new$functions))
      note(path, paste0("new file defines function ", name))
    if (length(new$others))
      note(path, "new file has top-level code other than function definitions")
  }
  # Symbols and strings in an expression, so do.call("name") is caught too.
  referenced <- function(e) {
    if (is.name(e)) return(as.character(e))
    if (is.character(e)) return(e)
    if (is.call(e) || is.pairlist(e) || is.expression(e) || is.list(e))
      return(unique(unlist(lapply(as.list(e), referenced))))
    character()
  }
  for (path in names(lesson_source_hashes())) {
    if (!is_r(path) || !length(allow)) next
    current <- definitions(path)
    for (name in names(current$functions)) {
      used <- intersect(referenced(current$functions[[name]]), allow)
      if (length(used))
        note(path, paste0("function ", name, " refers to allow-listed ",
                          paste(used, collapse = ", ")))
    }
    used <- intersect(referenced(current$others), allow)
    if (length(used))
      note(path, paste0("top-level code refers to allow-listed ",
                        paste(used, collapse = ", ")))
  }
  if (length(problems))
    stop("Package source differs from the recorded fingerprint beyond comments ",
         "and allowed plot-only functions:\n  ", paste(problems, collapse = "\n  "),
         call. = FALSE)
  invisible(TRUE)
}
