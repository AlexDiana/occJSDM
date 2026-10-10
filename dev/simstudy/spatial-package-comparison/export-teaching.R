moran_reporting_summary <- function(basis,prediction) {
  stopifnot(basis$arm=="mem")
  distance <- xy_distances(prediction$xy,basis$training_xy)
  nearest <- apply(distance,1,min)
  train <- basis$training_matrix
  future <- if(basis$status=="available") predict_basis(basis,prediction) else matrix(numeric(),length(prediction$site_id),0L)
  amplitude <- data.frame(column=colnames(train),training_sd=apply(train,2,sd),
    training_rms=sqrt(colMeans(train^2)),lattice_sd=apply(future,2,sd),lattice_rms=sqrt(colMeans(future^2)),
    training_maximum_absolute=apply(abs(train),2,max),lattice_maximum_absolute=apply(abs(future),2,max))
  list(graph=list(k=basis$graph$k,components=basis$graph$components,
    neighbor_radius_quantiles=quantile(basis$graph$radius,c(0,.25,.5,.75,1))),
    distance=list(units="raw unit-square Euclidean distance",nearest_training_quantiles=quantile(nearest,c(0,.25,.5,.75,1)),maximum_nearest_training_distance=max(nearest)),
    basis=amplitude,fitted_range_available=FALSE,
    interpretation="A fixed neighbor graph and training-scaled eigenvectors define this basis. Graph distances and basis amplitudes are not fitted kernel ranges or probability uncertainty.")
}

# Run only after full diagnostics selection, fresh native verification and truth scoring.
# Requires study input/runner/score helpers; the resulting bundle needs base R only.
build_teaching_bundle <- function(run_root,repo_root) {
  verification_path <- file.path(run_root,"verification-full.rds")
  if(!file.exists(verification_path)) stop("Export requires independent full-run verification")
  audit <- readRDS(verification_path); snapshot <- verify_snapshot(run_root)
  stopifnot(identical(audit$snapshot_hash,md5(file.path(run_root,"snapshot-record.rds"))),all(audit$checks$status %in% c("passed","unavailable")))
  cfg <- readRDS(file.path(run_root,"config.rds")); if(isTRUE(cfg$smoke_only)) stop("Synthetic runs cannot produce the lesson bundle")
  scores <- readRDS(file.path(run_root,"scores-full.rds"))
  stopifnot(identical(scores$snapshot_hash,audit$snapshot_hash),sum(scores$accounting$kind=="comparison")==192L,sum(scores$accounting$kind=="oracle")==96L)
  manifest <- resolve_manifest(readRDS(file.path(run_root,"import-manifest.rds")),run_root)
  training <- prediction_data <- truth_index <- landscapes <- list()
  for(ds in unique(manifest$dataset_id)) {
    read_kind <- function(kind) readRDS(manifest$path[manifest$dataset_id==ds & manifest$kind==kind])
    training[[ds]] <- read_kind("training"); prediction_data[[ds]] <- read_kind("prediction")
    truth <- read_kind("truth"); comm <- sub("-.*$","",ds)
    if(is.null(landscapes[[comm]])) landscapes[[comm]] <- truth$landscape else stopifnot(identical(landscapes[[comm]],truth$landscape))
    truth$landscape <- NULL; truth$community <- comm; truth_index[[ds]] <- truth
  }
  landscapes <- landscapes[cfg$communities]
  basis_summaries <- setNames(lapply(names(training),function(ds) {
    path <- file.path(run_root,"bases",paste0(ds,"-mem.rds"))
    if(file.exists(path)) moran_reporting_summary(readRDS(path),prediction_data[[ds]]) else NULL
  }),names(training))
  jobs <- make_jobs(run_root); predictions <- selections <- attempts <- hyperparameters <- site_attribution <- list()
  for(i in seq_len(nrow(jobs))) {
    job <- as.list(jobs[i,]); key <- job$key; directory <- file.path(run_root,"jobs",key)
    state <- load_selected_predictions(run_root,job)
    predictions[[key]] <- state[c("training","lattice")]; selections[[key]] <- state$selection
    site_attribution[[key]] <- list()
    for(scope in c("training","lattice")) if(!is.null(state[[scope]])) {
      score_path <- file.path(run_root,"scores",key,paste0(scope,".rds"))
      if(!file.exists(score_path)) stop("Missing attribution score archive: ",score_path)
      site_attribution[[key]][[scope]] <- readRDS(score_path)$sites
    }
    paths <- list.dirs(file.path(directory,"attempts"),recursive=FALSE,full.names=TRUE)
    attempts[[key]] <- lapply(paths,function(path) {
      result <- readRDS(file.path(path,"summary.rds")); result$training_probability <- NULL
      result$files <- readRDS(file.path(path,"complete.rds"))$files
      result$relative_path <- file.path("jobs",key,"attempts",basename(path)); result
    })
    hp <- file.path(run_root,"scores",key,"hyperparameters.rds")
    hyperparameters[[key]] <- if(file.exists(hp)) readRDS(hp) else NULL
  }
  input <- readRDS(file.path(run_root,"input-record.rds"))
  historical_path <- file.path(repo_root,"vignettes/teaching-data/spatial-lesson.rds")
  if(md5(historical_path)!=cfg$bundle_md5) stop("Historical teaching reference identity changed")
  historical <- readRDS(historical_path)
  old <- historical[c("oracle","oracle_lattice","fits","field_convergence")]
  old$label <- "Historical reference only; different fitted schedules and pointwise-median-then-centre field convention"
  bundle <- list(schema_version=1L,input_mode=input$mode,config=cfg,landscapes=landscapes,truth_index=truth_index,
    training=training,prediction_data=prediction_data,predictions=predictions,scores=scores,selections=selections,
    attempts=attempts,hyperparameters=hyperparameters,basis_summaries=basis_summaries,site_attribution=site_attribution,historical=old,
    capabilities=read.csv(file.path(run_root,"snapshot/study/capabilities.csv")),sources=read.csv(file.path(run_root,"snapshot/study/sources.csv")),
    verification=list(archive_required=audit,compact=NULL),
    provenance=list(input_limitations=input$limitations,reconstruction_evidence=input$evidence,
      relative_input_manifest=manifest[setdiff(names(manifest),"path")],source_snapshot=snapshot,
      original_input_files=input$files[c("identity","relative_path","md5")],runtime_hash=md5(file.path(run_root,"runtime/package-records.rds")),
      runtime_checkpoint=readRDS(file.path(run_root,"runtime-checkpoint.rds")),historical_bundle_md5=cfg$bundle_md5))
  repair <- file.path(run_root,"runtime-repair-v1/manifest.json")
  if(file.exists(repair)) {
    recovery <- readRDS(file.path(run_root,"runtime-recovery.rds"))
    bundle$provenance$runtime_repair <- list(amendment="E12",manifest=jsonlite::read_json(repair),manifest_md5=md5(repair),
      recovery_manifest_md5=md5(file.path(run_root,"runtime-recovery.rds")),
      recovered_attempts=substring(recovery$attempts,nchar(run_root)+2L),completion_hashes=unname(recovery$completion_hashes),
      administrative_extensions=sum(basename(recovery$attempts)=="02"),
      interpretation="These eight pilot extensions followed an administrative loaded-version validation failure; they are not an observed poor-mixing extension rate")
  }
  policy_paths <- c(E11="execution-policy.json",E12="execution-policy-runtime-repair.json",E14="window-20261008-2224/policy.json",E15="resume-20261009-completion/policy.json",E16="resume-20261009-eight-workers/policy.json")
  policy_paths <- policy_paths[file.exists(file.path(run_root,policy_paths))]
  bundle$provenance$execution_policies <- data.frame(amendment=names(policy_paths),relative_path=unname(policy_paths),md5=md5(file.path(run_root,policy_paths)),row.names=NULL)
  bundle$verification$compact <- verify_teaching_bundle(bundle)
  bundle
}

save_teaching_bundle <- function(bundle,path) {
  verify_teaching_bundle(bundle)
  if(file.exists(path)) stop("Existing teaching bundle preserved; choose a new verified export path")
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- tempfile(tmpdir=dirname(path)); on.exit(unlink(temporary),add=TRUE)
  saveRDS(bundle,temporary,compress="xz",version=3L)
  digest <- unname(tools::md5sum(temporary))
  stopifnot(file.rename(temporary,path))
  writeLines(digest,paste0(path,".md5"))
  invisible(path)
}
