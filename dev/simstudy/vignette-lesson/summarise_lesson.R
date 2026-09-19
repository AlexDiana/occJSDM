# Run from repository root after all fits: Rscript .../summarise_lesson.R OUTDIR
args <- commandArgs(TRUE); stopifnot(length(args)==1)
outdir <- normalizePath(args[1],mustWork=TRUE)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
input_path <- file.path(outdir,"input.rds")
b <- readRDS(input_path)
stopifnot(identical(b$source_hashes,lesson_source_hashes()))
obs <- lesson_observations(b)
stopifnot(identical(b$cases,select_lesson_cases(obs)))
scores <- list(); manifests <- list()
for(arm in c("perfect","default","alternative")) {
  long_path <- file.path(outdir,paste0(arm,"-long-fit.rds"))
  path <- if(file.exists(long_path)) long_path else file.path(outdir,paste0(arm,"-fit.rds"))
  r <- readRDS(path)
  stopifnot(identical(r$source_hashes,b$source_hashes),
            identical(r$input_md5,unname(tools::md5sum(input_path))))
  message("Scoring ",arm)
  scores[[arm]] <- score_lesson_fit(r,b)
  manifests[[arm]] <- list(file=basename(path),md5=unname(tools::md5sum(path)),
                           seed=r$seed,mcmc=r$mcmc,priors=r$priors,
                           warnings=r$warnings,session=r$session,
                           seconds=as.numeric(difftime(r$finished,r$started,units="secs")))
  saveRDS(scores[[arm]],file.path(outdir,paste0(arm,"-summary.rds")))
  rm(r); invisible(gc())
}
bind <- function(field) do.call(rbind,lapply(scores,`[[`,field))
compact <- list(schema=1L,input=b,observations=obs,cases=b$cases,
                cells=bind("cells"),groups=bind("groups"),rates=bind("rates"),
                samples=bind("samples"),diagnostics=lapply(scores,`[[`,"diagnostics"),
                manifests=manifests,input_md5=unname(tools::md5sum(input_path)),
                summary_source_hashes=tools::md5sum(c("dev/simstudy/vignette-lesson/score_lesson.R",
                                                       "dev/simstudy/vignette-lesson/summarise_lesson.R")),
                reconstruction_difference=lapply(scores,`[[`,"reconstruction_difference"))
dir.create("vignettes/teaching-data",showWarnings=FALSE)
saveRDS(compact,"vignettes/teaching-data/nonspatial-lesson.rds",compress="xz")
print(compact$groups)
for(arm in names(scores)) {
  d <- scores[[arm]]$diagnostics
  cat(arm,": parameter max Rhat",max(d$rhat,na.rm=TRUE),
      "; cell max Rhat",max(scores[[arm]]$cells$rhat,na.rm=TRUE),"\n")
}
