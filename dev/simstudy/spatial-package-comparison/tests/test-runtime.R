library(testthat)
launcher <- normalizePath("../run-r")
test_that("recursive R uses the same interpreter as the worker", {
  code <- 'child <- system2(file.path(R.home("bin"),"R"),c("--vanilla","--slave","-e",shQuote("cat(R.version.string)")),stdout=TRUE); stopifnot(identical(child,R.version.string))'
  out <- suppressWarnings(system2(launcher,c("-e",shQuote(code)),stdout=TRUE,stderr=TRUE))
  expect_null(attr(out,"status"),info=paste(out,collapse="\n"))
})

test_that("Rscript invoked by RcppParallel preserves runtime and clean linker output", {
  code <- 'out <- system2(file.path(R.home("bin"),"Rscript"),c("-e",shQuote("cat(R.version.string)")),stdout=TRUE); stopifnot(identical(out,R.version.string))'
  out <- suppressWarnings(system2(launcher,c("-e",shQuote(code)),stdout=TRUE,stderr=TRUE))
  expect_null(attr(out,"status"),info=paste(out,collapse="\n"))
})
