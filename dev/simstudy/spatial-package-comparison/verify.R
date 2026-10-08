# Draw verification opens no truth and can run before scientific scoring.
verify_prediction_archive <- function(run_root,job) {
  directory <- file.path(run_root,"jobs",job$key); state <- load_selected_predictions(run_root,job)
  if(is.null(state$training)) return(data.frame(key=job$key,check="prediction",status="unavailable",difference=NA_real_))
  selected <- file.path(directory,"attempts",sprintf("%02d",state$selection$selected_attempt))
  saved <- readRDS(file.path(selected,"fit.rds")); f <- saved$fit; b <- saved$basis
  tr <- readRDS(file.path(run_root,"training",paste0(job$dataset_id,".rds"))); tr$y <- NULL
  if(job$kind=="oracle") tr$species_id <- job$species
  old <- options(spatial.draw_sink=NULL,spatial.retain_raw_field=!is.null(state$training$historical_field_median_centered)); on.exit(options(old))
  p <- get(paste0("predict_",tolower(f$package)))(f,tr,b)
  compare <- function(a,b,fields) max(vapply(fields,function(nm) max(abs(a[[nm]]-b[[nm]])),numeric(1)))
  columns <- c("p_full","p_space_zero","p_env_zero","p_both_zero","field_mean")
  difference <- compare(p,state$training,columns)
  if(difference>1e-10) stop("Fresh training prediction mismatch: ",job$key)
  rows <- list(data.frame(key=job$key,check="fresh training reconstruction",status="passed",difference=difference))
  options(spatial.retain_raw_field=FALSE)
  if(!is.null(state$lattice)) {
    pred <- readRDS(file.path(run_root,"prediction",paste0(job$dataset_id,".rds")))
    if(job$kind=="oracle") pred$species_id <- job$species
    if(f$package %in% c("sjSDM","gllvm")) {
      new <- get(paste0("predict_",tolower(f$package)))(f,pred,b)
      difference <- compare(new,state$lattice,columns)
      if(difference>1e-10) stop("Fresh point prediction mismatch: ",job$key)
      rows[[2]] <- data.frame(key=job$key,check="fresh lattice reconstruction",status="passed",difference=difference)
    } else {
      cfg <- readRDS(file.path(run_root,"config.rds")); expected <- pred$site_id[cfg$audit_lattice_cells[cfg$audit_lattice_cells<=length(pred$site_id)]]
      if(f$package %in% c("Hmsc","oracle")) pred <- subset_prediction(pred,order(pred$site_id,method="radix"))
      ex <- get(paste0(if(f$package=="occJSDM") "occ" else tolower(f$package),"_draw_extractor"))(f,pred,b)
      files <- list.files(file.path(directory,"prediction/audit"),pattern="rds$",full.names=TRUE)
      observed <- character()
      for(file in files) {
        audit <- readRDS(file); observed <- c(observed,audit$site_id)
        stopifnot(!is.null(audit$block_site_id),identical(audit$species_id,pred$species_id))
        if(!is.null(audit$rng_before)) assign(".Random.seed",audit$rng_before,globalenv())
        idx <- match(audit$block_site_id,pred$site_id); stopifnot(!anyNA(idx),length(idx)<=25L)
        d <- ex(idx); keep <- match(audit$site_id,audit$block_site_id)
        difference <- max(vapply(names(audit$draws),function(nm) max(abs(d[[nm]][keep,,,,drop=FALSE]-audit$draws[[nm]])),numeric(1)))
        if(difference>1e-10) stop("Independent audit-draw reconstruction mismatch: ",file)
        loc <- match(audit$site_id,state$lattice$site_id)
        for(nm in columns[1:4]) {
          error <- max(abs(apply(audit$draws[[nm]],c(1,2),mean)-state$lattice[[nm]][loc,,drop=FALSE]))
          if(error>1e-10) stop("Audit draws disagree with lattice mean")
        }
      }
      if(!setequal(observed,expected)||anyDuplicated(observed)) stop("Fixed lattice audit IDs are incomplete or duplicated")
      rows[[2]] <- data.frame(key=job$key,check="native draw reconstruction and fixed audit means",status="passed",difference=0)
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
