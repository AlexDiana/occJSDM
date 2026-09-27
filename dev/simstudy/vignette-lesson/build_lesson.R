# Run from repository root: Rscript dev/simstudy/vignette-lesson/build_lesson.R OUTDIR prepare|perfect|default|alternative [long]
args <- commandArgs(TRUE)
stopifnot(length(args)>=2)
outdir <- normalizePath(args[1],mustWork=TRUE)
arm <- args[2]
source("dev/simstudy/vignette-lesson/helpers.R")
input_path <- file.path(outdir,"input.rds")
if(arm=="prepare") {
  stopifnot(!file.exists(input_path))
  bundle <- make_lesson()
  bundle$source_hashes <- lesson_source_hashes()
  bundle$source_commit <- system2("git",c("rev-parse","HEAD"),stdout=TRUE)
  bundle$generator_hashes <- tools::md5sum(c("dev/simstudy/vignette-lesson/helpers.R",
                                           "dev/simstudy/vignette-lesson/build_lesson.R"))
  saveRDS(bundle,input_path,compress="xz")
  write.csv(bundle$cases,file.path(outdir,"cases-selected-before-fitting.csv"),row.names=FALSE)
  print(bundle$cases)
} else {
  bundle <- readRDS(input_path)
  stopifnot(identical(bundle$source_hashes,lesson_source_hashes()))
  long <- length(args)>=3 && args[3]=="long"
  path <- file.path(outdir,paste0(arm,if(long)"-long" else "","-fit.rds"))
  stopifnot(!file.exists(path))
  mcmc <- if(long) list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L) else
                   list(nchain=4L,nburn=3000L,niter=6000L,nthin=1L)
  result <- lesson_fit(bundle,arm,mcmc)
  result$input_md5 <- unname(tools::md5sum(input_path))
  saveRDS(result,path,compress="gzip")
  cat("Saved",path,"in",as.numeric(difftime(result$finished,result$started,units="secs")),"seconds\n")
}
