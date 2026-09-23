# Evidence checks independent of the rendered tables. Run from repository root.
args <- commandArgs(TRUE); stopifnot(length(args)==1)
outdir <- normalizePath(args[1],mustWork=TRUE)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
x <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
b <- readRDS(file.path(outdir,"input.rds"))
stopifnot(identical(x$input,b),identical(b$source_hashes,lesson_source_hashes()),
          identical(b$generator_hashes,tools::md5sum(names(b$generator_hashes))),
          identical(x$summary_source_hashes,tools::md5sum(names(x$summary_source_hashes))),
          identical(x$input_md5,unname(tools::md5sum(file.path(outdir,"input.rds")))))
binary <- lesson_binary_data(b)
stopifnot(identical(unname(binary$OTU),unname(b$sim$true_params$z_true)),
          !any(c("Sample","Primer") %in% names(binary$info)),
          nrow(binary$OTU)==100,ncol(binary$OTU)==10,
          max(table(interaction(b$sim$data_list$info$Sample,b$sim$data_list$info$Primer)))==6)
stopifnot(identical(x$cases,select_lesson_cases(lesson_observations(b))))
o <- x$observations
stopifnot(all(o$w[o$source=="Laboratory false positive"]==0),
          all(o$z[o$source=="Field-stage false positive"]==0),
          all(o$w[o$source=="Field-stage false positive"]==1),
          all(o$z[o$source=="True detection"]==1),
          all(o$w[o$source=="True detection"]==1))
for(arm in names(x$manifests)) {
  manifest <- x$manifests[[arm]]
  path <- file.path(outdir,manifest$file)
  stopifnot(identical(manifest$md5,unname(tools::md5sum(path))))
  r <- readRDS(path); f <- r$fit
  validate_lesson_fit_identity(f,b)
  stopifnot(identical(r$source_hashes,b$source_hashes),
            identical(manifest$session,r$session),
            identical(r$input_md5,x$input_md5),f$infos$ps==0)
  c <- x$cells[x$cells$arm==arm,]
  j <- f$results_output$jsdm_output
  ni <- dim(j$B0_output)[2]; nc <- dim(j$B0_output)[3]
  # Independently reconstruct all draws for one complete species, using the
  # component sums rather than the scoring helper's matrix implementation.
  check <- matrix(0,length(f$infos$siteNames),ni*nc)
  column <- 0L
  for(ch in seq_len(nc))for(it in seq_len(ni)) {
    column <- column+1L
    linear <- rep(j$B0_output[1,it,ch],nrow(check))
    for(k in seq_len(ncol(f$X_psi))) linear <- linear+f$X_psi[,k]*j$B_output[k,1,it,ch]
    for(k in seq_len(dim(j$U_output)[2])) linear <- linear+j$U_output[,k,it,ch]*j$L_output[k,1,it,ch]
    check[,column] <- plogis(linear)
  }
  cc <- c[c$species==f$infos$speciesNames[1],]
  stopifnot(max(abs(rowMeans(check)-cc$estimate))<1e-10,
            max(abs(apply(check,1,quantile,.025)-cc$lower))<1e-10,
            max(abs(apply(check,1,quantile,.975)-cc$upper))<1e-10)
  if(arm!="perfect") {
    stopifnot(max(abs(c$estimate-c(f$results_output$psi_output)))<1e-10,
              max(abs(c$conditional-c(f$results_output$z_output)))<1e-10)
    samples <- x$samples[x$samples$arm==arm,]
    stopifnot(max(abs(samples$sample_probability-c(f$results_output$w_output)))<1e-10)
  }
  for(band in c("All","Low","Middle","High")) {
    keep <- switch(band,All=rep(TRUE,nrow(c)),Low=c$truth<.2,
                   Middle=c$truth>=.2 & c$truth<=.8,High=c$truth>.8)
    g <- x$groups[x$groups$arm==arm & x$groups$band==band,]
    stopifnot(abs(g$mae-mean(abs(c$estimate[keep]-c$truth[keep])))<1e-12,
              abs(g$signed_error-mean(c$estimate[keep]-c$truth[keep]))<1e-12)
  }
  cat(arm,": verified fit hashes, identities, independent draws, means, intervals and errors\n")
  rm(r,f,j,check);invisible(gc())
}
cat("All teaching-evidence checks passed.\n")
