# Namespace loading changes session-wide ggplot2 state. Use the actual R
# executable in a fresh process, including when another test has loaded ggtern.
.plotting_package_cache <- new.env(parent = emptyenv())

installed_plotting_package <- function() {
  path <- getNamespaceInfo(asNamespace("occJSDM"), "path")
  if (file.exists(file.path(path, "Meta", "package.rds"))) return(path)

  # pkgload checks every DESCRIPTION import with requireNamespace(), including
  # ggtern. Its development loader cannot reproduce ordinary namespace loading.
  # Install once per source-test session; R removes this temporary library when
  # the session exits. R CMD check already supplies an installed package.
  if (is.null(.plotting_package_cache$path)) {
    test_library <- tempfile("occJSDM-plot-tests-")
    dir.create(test_library)
    command <- file.path(R.home("bin"),
                         if (.Platform$OS.type == "windows") "Rcmd.exe" else "R")
    args <- c(if (.Platform$OS.type != "windows") "CMD", "INSTALL",
              "--no-multiarch", "--no-docs", paste0("--library=", shQuote(test_library)),
              shQuote(path))
    env <- c("R_TESTS=", paste0("R_LIBS=", shQuote(
      paste(.libPaths(), collapse = .Platform$path.sep)
    )))
    log <- system2(command, args, stdout = TRUE, stderr = TRUE, env = env)
    status <- attr(log, "status")
    if (!is.null(status) && status != 0L) {
      stop(paste(log, collapse = "\n"), call. = FALSE)
    }
    .plotting_package_cache$path <- file.path(test_library, "occJSDM")
  }
  .plotting_package_cache$path
}

fresh_plotting_session <- function(fun) {
  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(input, output, script)), add = TRUE)
  saveRDS(list(code = deparse(fun), libpath = .libPaths(),
               package_path = installed_plotting_package()),
          input)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    "input <- readRDS(args[1])",
    ".libPaths(input$libpath)",
    "load_occJSDM <- function() {",
    "  loadNamespace('occJSDM', lib.loc = dirname(input$package_path))",
    "}",
    "render_plot <- function(plot) {",
    "  tryCatch({ print(plot); '' }, error = function(e) conditionMessage(e))",
    "}",
    "grDevices::pdf(NULL)",
    "result <- eval(parse(text = input$code))()",
    "grDevices::dev.off()",
    "saveRDS(result, args[2])"
  ), script)
  executable <- if (.Platform$OS.type == "windows") {
    file.path(R.home("bin"), "Rterm.exe")
  } else {
    file.path(R.home("bin"), "exec", "R")
  }
  log <- system2(executable,
                 c("--vanilla", "--slave", paste0("--file=", shQuote(script)),
                   "--args", shQuote(input), shQuote(output)),
                 stdout = TRUE, stderr = TRUE, env = "R_TESTS=")
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(log, collapse = "\n"), call. = FALSE)
  }
  readRDS(output)
}

variance_partitioning_plot_fixture <- function() {
  list(
    infos = list(speciesNames = c("OTU_1", "OTU_2")),
    results_output = list(jsdm_output = list(
      varPart_output = array(rep(c(.2, .4, .5, .3, .3, .3, 1, 1), 6),
                             dim = c(2, 4, 3, 2))
    ))
  )
}
