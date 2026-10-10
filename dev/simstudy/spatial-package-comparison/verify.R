# Draw verification opens no truth and can run before scientific scoring.
verification_difference <- function(actual,expected,label,tolerance=1e-10) {
  if(!identical(dim(actual),dim(expected)) || length(actual)!=length(expected) ||
     any(!is.finite(actual)) || any(!is.finite(expected))) stop("Prediction summary mismatch: ",label)
  difference <- max(abs(actual-expected))
  if(difference>tolerance) stop("Prediction summary mismatch: ",label," (",difference,")")
  difference
}

# Independently reduce the archived unthinned draws, without the block summarizer.
verify_draw_summaries <- function(audit,prediction,tolerance=1e-10) {
  loc <- match(audit$site_id,prediction$site_id)
  if(anyNA(loc)||anyDuplicated(audit$site_id)||!identical(audit$species_id,prediction$species_id)) stop("Audit IDs mismatch")
  d <- audit$draws; shape <- dim(d$p_full)
  if(length(shape)!=4L || !identical(shape[1:2],c(length(loc),length(audit$species_id))) ||
     !all(vapply(d,function(x) identical(dim(x),shape),logical(1)))) stop("Audit draw dimensions mismatch")
  difference <- 0
  compare <- function(a,b,label) difference <<- max(difference,verification_difference(a,b,label,tolerance))
  for(nm in c("p_full","p_space_zero","p_env_zero","p_both_zero"))
    compare(apply(d[[nm]],c(1,2),mean),prediction[[nm]][loc,,drop=FALSE],nm)
  compare(apply(d$field,c(1,2),mean),prediction$field_mean[loc,,drop=FALSE],"field_mean")
  compare(apply(d$field,c(1,2),median),prediction$field_median[loc,,drop=FALSE],"field_median")
  compare(apply(d$p_full-d$p_space_zero,c(1,2),mean),prediction$delta_p[loc,,drop=FALSE],"delta_p")
  for(target in c("probability","field")) {
    values <- d[[if(target=="probability") "p_full" else "field"]]
    for(bound in c("lower","upper")) {
      nm <- paste(target,bound,sep="_")
      q <- apply(values,c(1,2),quantile,probs=if(bound=="lower") .025 else .975,type=7L)
      compare(q,prediction$intervals[[nm]][loc,,drop=FALSE],nm)
    }
  }
  for(nm in c("p_full","p_space_zero","p_env_zero","p_both_zero","delta_p","field")) {
    values <- if(nm=="delta_p") d$p_full-d$p_space_zero else d[[nm]]
    # Explicit site/species loops check iteration-by-chain shape and constant draws.
    errors <- matrix(0,length(loc),length(audit$species_id))
    for(i in seq_along(loc)) for(j in seq_along(audit$species_id)) {
      x <- matrix(values[i,j,,],shape[3],shape[4])
      errors[i,j] <- if(all(x==x[1])) 0 else suppressWarnings(posterior::mcse_mean(x))
    }
    compare(errors,prediction$mcse[[nm]][loc,,drop=FALSE],paste("MCSE",nm))
  }
  difference
}

compare_prediction_summaries <- function(actual,expected) {
  stopifnot(identical(actual$site_id,expected$site_id),identical(actual$species_id,expected$species_id),
    identical(actual$link,expected$link),identical(actual$uncertainty,expected$uncertainty))
  fields <- c("p_full","p_space_zero","p_env_zero","p_both_zero","delta_p","field_mean","field_median")
  differences <- vapply(fields,function(nm) verification_difference(actual[[nm]],expected[[nm]],nm),numeric(1))
  for(group in c("intervals","mcse")) {
    if(!identical(names(actual[[group]]),names(expected[[group]]))) stop("Prediction summary capability mismatch: ",group)
    if(length(actual[[group]])) differences <- c(differences,vapply(names(actual[[group]]),function(nm)
      verification_difference(actual[[group]][[nm]],expected[[group]][[nm]],paste(group,nm)),numeric(1)))
  }
  if(!is.null(expected$historical_field_median_centered)) differences <- c(differences,
    verification_difference(actual$historical_field_median_centered,expected$historical_field_median_centered,"historical field"))
  max(differences)
}

verify_prediction_archive <- function(run_root,job) {
  directory <- file.path(run_root,"jobs",job$key); state <- load_selected_predictions(run_root,job)
  if(is.null(state$training)) return(data.frame(key=job$key,check="prediction",status="unavailable",difference=NA_real_))
  selected <- file.path(directory,"attempts",sprintf("%02d",state$selection$selected_attempt))
  saved <- readRDS(file.path(selected,"fit.rds")); f <- saved$fit; b <- saved$basis
  tr <- readRDS(file.path(run_root,"training",paste0(job$dataset_id,".rds"))); tr$y <- NULL
  if(job$kind=="oracle") tr$species_id <- job$species
  old <- options(spatial.draw_sink=NULL,spatial.retain_raw_field=!is.null(state$training$historical_field_median_centered)); on.exit(options(old))
  p <- get(paste0("predict_",tolower(f$package)))(f,tr,b)
  difference <- compare_prediction_summaries(p,state$training)
  if(difference>1e-10) stop("Fresh training prediction mismatch: ",job$key)
  rows <- list(data.frame(key=job$key,check="fresh training reconstruction",status="passed",difference=difference))
  options(spatial.retain_raw_field=FALSE)
  if(!is.null(state$lattice)) {
    pred <- readRDS(file.path(run_root,"prediction",paste0(job$dataset_id,".rds")))
    if(job$kind=="oracle") pred$species_id <- job$species
    if(f$package %in% c("sjSDM","gllvm")) {
      new <- get(paste0("predict_",tolower(f$package)))(f,pred,b)
      difference <- compare_prediction_summaries(new,state$lattice)
      if(difference>1e-10) stop("Fresh point prediction mismatch: ",job$key)
      rows[[2]] <- data.frame(key=job$key,check="fresh lattice reconstruction",status="passed",difference=difference)
    } else {
      cfg <- readRDS(file.path(run_root,"config.rds")); expected <- pred$site_id[cfg$audit_lattice_cells[cfg$audit_lattice_cells<=length(pred$site_id)]]
      if(f$package %in% c("Hmsc","oracle")) pred <- subset_prediction(pred,order(pred$site_id,method="radix"))
      ex <- get(paste0(if(f$package=="occJSDM") "occ" else tolower(f$package),"_draw_extractor"))(f,pred,b)
      files <- list.files(file.path(directory,"prediction/audit"),pattern="rds$",full.names=TRUE)
      observed <- character()
      largest <- 0
      for(file in files) {
        audit <- readRDS(file); observed <- c(observed,audit$site_id)
        stopifnot(!is.null(audit$block_site_id),identical(audit$species_id,pred$species_id))
        stopifnot(identical(audit$draw_ids$iteration,seq_len(dim(audit$draws$p_full)[3])),
          identical(audit$draw_ids$chain,seq_len(dim(audit$draws$p_full)[4])))
        if(!is.null(audit$rng_before)) assign(".Random.seed",audit$rng_before,globalenv())
        idx <- match(audit$block_site_id,pred$site_id); stopifnot(!anyNA(idx),length(idx)<=25L)
        d <- ex(idx); keep <- match(audit$site_id,audit$block_site_id)
        difference <- max(vapply(names(audit$draws),function(nm) max(abs(d[[nm]][keep,,,,drop=FALSE]-audit$draws[[nm]])),numeric(1)))
        if(difference>1e-10) stop("Independent audit-draw reconstruction mismatch: ",file)
        largest <- max(largest,difference,verify_draw_summaries(audit,state$lattice))
      }
      if(!setequal(observed,expected)||anyDuplicated(observed)) stop("Fixed lattice audit IDs are incomplete or duplicated")
      rows[[2]] <- data.frame(key=job$key,check="native draws, fixed audit IDs, means, intervals and MCSE",status="passed",difference=largest)
    }
  }
  do.call(rbind,rows)
}
verify_results <- function(run_root,repo_root=NULL,phase="full") {
  verify_snapshot(run_root); jobs <- make_jobs(run_root)
  if(phase=="pilot") jobs <- subset(jobs,community=="rep01"&arrangement %in% c("spread","clustered"))
  out <- lapply(seq_len(nrow(jobs)),function(i) verify_prediction_archive(run_root,as.list(jobs[i,])))
  record <- list(schema_version=1L,checks=do.call(rbind,out),snapshot_hash=md5(file.path(run_root,"snapshot-record.rds")),phase=phase)
  write_once(record,file.path(run_root,paste0("verification-",phase,".rds"))); record$checks
}
