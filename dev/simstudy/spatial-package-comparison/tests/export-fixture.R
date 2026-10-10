# Synthetic archive of 288 accounted-for cells with one tiny saved prediction.
# No models are fitted here; every file is written under a new temporary directory.
fixture_export_archive <- function(bundle,code,repo) {
  root <- tempfile(); dir.create(file.path(root,"sources/acceptance"),recursive=TRUE)
  writeLines("synthetic fixture",file.path(root,"sources/acceptance/generator.R"))
  cfg <- bundle$config; write_once(cfg,file.path(root,"config.rds"))
  write_once(comparison_manifest(cfg,"reconstructed"),file.path(root,"comparison-manifest.rds"))
  write_once(oracle_manifest(cfg,"reconstructed"),file.path(root,"oracle-manifest.rds"))
  manifest <- list()
  for(ds in names(bundle$training)) {
    truth <- bundle$truth_index[[ds]]; truth$landscape <- bundle$landscapes[[truth$community]]
    values <- list(training=bundle$training[[ds]],prediction=bundle$prediction_data[[ds]],truth=truth)
    for(kind in names(values)) {
      relative <- paste0(kind,"/",ds,".rds"); path <- file.path(root,relative); write_once(values[[kind]],path)
      manifest[[length(manifest)+1L]] <- data.frame(identity=paste(ds,kind,sep="-"),dataset_id=ds,kind=kind,relative_path=relative,md5=md5(path))
    }
  }
  write_once(do.call(rbind,manifest),file.path(root,"import-manifest.rds"))
  input <- list(mode="reconstructed",limitations="Synthetic exporter test only",evidence=list(synthetic=TRUE),
    files=data.frame(identity=cfg$communities,relative_path=paste0("inputs/",cfg$communities,".rds"),md5="synthetic fixture"))
  write_once(input,file.path(root,"input-record.rds"))
  write_once(list(synthetic=TRUE),file.path(root,"runtime/package-records.rds"))
  write_once(list(synthetic=TRUE,passed=TRUE),file.path(root,"runtime-checkpoint.rds"))
  jobs <- make_jobs(root); accounting <- transform(jobs,status="failed",training_supported=FALSE,lattice_supported=FALSE)
  selected_key <- names(bundle$predictions)[1]
  for(i in seq_len(nrow(jobs))) {
    job <- jobs[i,]; dir <- file.path(root,"jobs",job$key)
    identity <- c(synthetic="export unit fixture")
    selection <- list(identity=identity,status="failed",selected_attempt=NA_integer_,attempts=integer(),reason="synthetic failed cell")
    if(job$key==selected_key) {
      ap <- file.path(dir,"attempts/01")
      write_once(bundle$predictions[[selected_key]]$training,file.path(ap,"training-prediction.rds"))
      complete_attempt(ap,identity,list(status="complete",attempt_id=1L,acceptable=TRUE))
      selection$status <- "selected"; selection$selected_attempt <- 1L; selection$attempts <- 1L
      selection$attempt_hashes <- setNames(md5(file.path(ap,"complete.rds")),"1")
      accounting$status[i] <- "selected"; accounting$training_supported[i] <- TRUE
      tv <- bundle$truth_index[[job$dataset_id]]; tv$landscape <- bundle$landscapes[[tv$community]]
      write_once(score_predictions(bundle$predictions[[selected_key]]$training,truth_view(tv,"training")),file.path(root,"scores",job$key,"training.rds"))
    }
    write_once(selection,file.path(dir,"selection.rds"))
  }
  freeze_sources(code,repo,root)
  hash <- md5(file.path(root,"snapshot-record.rds"))
  write_once(list(snapshot_hash=hash,checks=data.frame(key=selected_key,check="fresh training reconstruction",status="passed",difference=0)),file.path(root,"verification-full.rds"))
  scores <- bundle$scores; scores$accounting <- accounting; scores$snapshot_hash <- hash
  write_once(scores,file.path(root,"scores-full.rds")); root
}
