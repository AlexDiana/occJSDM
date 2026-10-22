# With a file destination, stderr must also name that file. stderr=TRUE makes
# system2 capture text and return a character vector instead of an exit code.
run_logged_r <- function(script,args,log) {
  system2(file.path(R.home('bin'),'Rscript'),shQuote(c(script,args)),
    stdout=log,stderr=log)
}
