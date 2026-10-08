library(testthat)
test_that("the detached wrapper closes stdin and survives the starting shell", {
  root <- tempfile(); launcher <- normalizePath("../launch.sh")
  code <- "import os,time,pathlib; assert os.read(0,1)==b''; time.sleep(1); print('dummy finished',flush=True)"
  result <- system2(launcher,c(shQuote(root),"/usr/bin/python3","-c",shQuote(code)),stdout=TRUE,stderr=TRUE)
  expect_null(attr(result,"status"),info=paste(result,collapse="\n"))
  deadline <- Sys.time()+15
  while(!file.exists(file.path(root,"complete.json")) && Sys.time()<deadline) Sys.sleep(.1)
  expect_true(file.exists(file.path(root,"complete.json")))
  if(!file.exists(file.path(root,"complete.json"))) return(invisible(NULL))
  done <- jsonlite::read_json(file.path(root,"complete.json"),simplifyVector=TRUE)
  expect_equal(done$exit_status,0)
  expect_match(paste(readLines(file.path(root,"stdout.log")),collapse="\n"),"dummy finished")
  expect_match(paste(readLines(file.path(root,"launcher.log")),collapse="\n"),"SPATIAL COMMAND COMPLETE")
})
test_that("interrupted jobs leave status and deterministic job-level resume preserves completed work", {
  root <- tempfile(); launcher <- normalizePath("../launch.sh"); work <- tempfile(); dir.create(work)
  artifact <- file.path(work,"result.txt")
  code <- paste0("import pathlib,time; p=pathlib.Path(",encodeString(artifact,quote='"'),
    "); p.write_text('partial'); time.sleep(20); p.write_text('complete')")
  launch <- system2(launcher,c(shQuote(root),"/usr/bin/python3","-c",shQuote(code)),stdout=TRUE)
  deadline <- Sys.time()+10
  while(!file.exists(artifact) && Sys.time()<deadline) Sys.sleep(.1)
  expect_true(file.exists(artifact))
  owner <- jsonlite::read_json(file.path(root,"owner.json"),simplifyVector=TRUE)
  # This test signals only the supervisor it just created and recorded.
  tools::pskill(owner$pid,15L)
  deadline <- Sys.time()+15
  while(!file.exists(file.path(root,"complete.json")) && Sys.time()<deadline) Sys.sleep(.1)
  expect_true(file.exists(file.path(root,"complete.json")))
  expect_equal(readLines(artifact,warn=FALSE),"partial")
  resume <- paste0("import pathlib; p=pathlib.Path(",encodeString(artifact,quote='"'),
    "); p.write_text('complete') if p.read_text()!='complete' else None")
  for(i in 1:2) {
    log <- tempfile()
    system2(launcher,c(shQuote(log),"/usr/bin/python3","-c",shQuote(resume)),stdout=TRUE)
    deadline <- Sys.time()+10
    while(!file.exists(file.path(log,"complete.json")) && Sys.time()<deadline) Sys.sleep(.1)
    expect_equal(jsonlite::read_json(file.path(log,"complete.json"),simplifyVector=TRUE)$exit_status,0L)
    now <- file.info(artifact)$mtime
    if(i==1) first <- now else expect_equal(now,first)
  }
  expect_equal(readLines(artifact,warn=FALSE),"complete")
})
