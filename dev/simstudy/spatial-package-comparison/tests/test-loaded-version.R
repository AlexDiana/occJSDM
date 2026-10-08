library(testthat)
source("../archive.R")
source("../runner.R")
test_that("loaded hyphenated package versions match their numeric version", {
  loadNamespace("Hmsc")
  path <- find.package("Hmsc")
  files <- list.files(path,recursive=TRUE,full.names=TRUE)
  record <- list(package="Hmsc",path=path,version=as.character(packageVersion("Hmsc")),files=setNames(md5(files),files))
  expect_true(verify_package_identity(list(record)))
  record$version <- "0.0.0"
  expect_error(verify_package_identity(list(record)),"version mismatch")
})
