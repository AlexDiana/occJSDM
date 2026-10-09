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
  jobs <- make_jobs(run_root); predictions <- selections <- attempts <- hyperparameters <- list()
  for(i in seq_len(nrow(jobs))) {
    job <- as.list(jobs[i,]); key <- job$key; directory <- file.path(run_root,"jobs",key)
    state <- load_selected_predictions(run_root,job)
    predictions[[key]] <- state[c("training","lattice")]; selections[[key]] <- state$selection
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
    attempts=attempts,hyperparameters=hyperparameters,historical=old,
    capabilities=read.csv(file.path(run_root,"snapshot/study/capabilities.csv")),sources=read.csv(file.path(run_root,"snapshot/study/sources.csv")),
    verification=list(archive_required=audit,compact=NULL),
    provenance=list(input_limitations=input$limitations,reconstruction_evidence=input$evidence,
      relative_input_manifest=manifest[setdiff(names(manifest),"path")],source_snapshot=snapshot,
      original_input_files=input$files[c("identity","relative_path","md5")],runtime_hash=md5(file.path(run_root,"runtime/package-records.rds")),
      runtime_checkpoint=readRDS(file.path(run_root,"runtime-checkpoint.rds")),historical_bundle_md5=cfg$bundle_md5))
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
