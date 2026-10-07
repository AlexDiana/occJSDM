library(testthat)
source('execution.R')
test_that('logged R jobs preserve stdout, stderr, arguments and numeric exit status', {
  folder<-tempfile('R logging ');dir.create(folder)
  script<-file.path(folder,'script with spaces.R');log<-file.path(folder,'output.log')
  writeLines(c('cat("ordinary output\\n")','cat("error output\\n",file=stderr())',
    'cat(commandArgs(TRUE)[1],"\\n")','quit(status=as.integer(commandArgs(TRUE)[2]))'),script)
  status<-run_logged_r(script,c('argument with spaces','0'),log)
  expect_identical(status,0L)
  expect_true(file.exists(log))
  expect_match(paste(readLines(log),collapse=' '),'ordinary output.*error output.*argument with spaces')
  expect_identical(run_logged_r(script,c('failure case','7'),log),7L)
})
