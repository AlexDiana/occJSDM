# Input operations run only in preprocessing, never in competing fit workers.
md5 <- function(path) unname(tools::md5sum(path))
write_once <- function(object,path) {
  if (file.exists(path)) {
    if (!identical(readRDS(path),object)) stop("Existing different result preserved: ",path)
    return(invisible(path))
  }
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(pattern=".pending-",tmpdir=dirname(path))
  on.exit(unlink(tmp),add=TRUE); saveRDS(object,tmp,version=3L)
  if (!file.rename(tmp,path)) stop("Atomic write failed: ",path)
  invisible(path)
}
resolve_manifest <- function(manifest,root) {
  required <- c("identity","relative_path","md5")
  if (!all(required %in% names(manifest))) stop("Incomplete input manifest")
  if (anyDuplicated(manifest$identity) || anyDuplicated(manifest$relative_path)) stop("duplicate manifest identity/path")
  rel <- manifest$relative_path
  if (anyNA(rel) || any(!nzchar(rel)) || any(grepl("^/|(^|/)\\.\\.(/|$)",rel))) stop("Unsafe relative path")
  root <- normalizePath(root,mustWork=TRUE)
  manifest$path <- file.path(root,rel)
  if (any(!file.exists(manifest$path))) stop("archive prerequisite: missing manifest files")
  real <- normalizePath(manifest$path,mustWork=TRUE)
  if (any(!startsWith(real,paste0(root,"/")))) stop("Unsafe relative path escapes archive")
  if (any(md5(real)!=manifest$md5)) stop("Input hash mismatch")
  manifest
}
verify_archive <- function(archive_root,repo_root) {
  path <- file.path(archive_root,"input-manifest.rds")
  if (!file.exists(path)) stop("archive prerequisite: original input-manifest.rds and inputs are required")
  old <- readRDS(path)
  if (!all(c("community","arrangement","arm","input_md5") %in% names(old)) ||
      nrow(old)!=24L || anyDuplicated(old[c("community","arrangement","arm")]))
    stop("archive prerequisite: incomplete or duplicate historical manifest")
  f <- unique(old[c("community","input_md5")])
  files <- resolve_manifest(data.frame(identity=f$community,
    relative_path=paste0("inputs/",f$community,".rds"),md5=f$input_md5),archive_root)
  if (!setequal(files$identity,sprintf("rep%02d",1:3))) stop("archive prerequisite: all three inputs required")
  files <- files[match(sprintf("rep%02d",1:3),files$identity),]
  catalog <- read.csv(file.path(repo_root,"dev/maintenance/archive-review-20261007/lesson-7-backup-paths.csv"))
  checks <- list()
  for (col in c("expected_selected_result_md5","expected_selected_oracle_md5")) {
    expected <- catalog[!is.na(catalog[[col]]) & nzchar(catalog[[col]]),]
    present <- file.exists(file.path(archive_root,expected$relative_path))
    if (any(md5(file.path(archive_root,expected$relative_path[present])) != expected[[col]][present])) stop("Historical reference hash mismatch")
    checks[[col]] <- data.frame(path=expected$relative_path,present=present,expected_md5=expected[[col]])
  }
  evidence <- audit_inputs(lapply(files$path,readRDS),repo_root)
  list(schema_version=1L,mode="restored",files=files,evidence=evidence,
    historical_checks=checks,limitations="Input hashes verified against restored manifest; independent historical checksums exist only for listed reference fits.")
}
reconstruct_inputs <- function(repo_root,output_root) {
  if (file.exists(file.path(output_root,"input-record.rds")) || dir.exists(file.path(output_root,"inputs")))
    stop("Reconstruction destination already exists; verify and resume without regeneration")
  cfg <- study_config(); generator <- file.path(repo_root,"dev/simstudy/spatial-design-sweep/generator.R")
  if (!identical(md5(generator),cfg$generator_md5)) stop("generator hash mismatch")
  gen <- new.env(parent=globalenv()); sys.source(generator,gen)
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  inputs <- lapply(1:3,gen$make_sweep_input)
  evidence <- audit_inputs(inputs,repo_root)
  files <- data.frame(identity=sprintf("rep%02d",1:3),relative_path=sprintf("inputs/rep%02d.rds",1:3))
  files$path <- file.path(normalizePath(output_root,mustWork=TRUE),files$relative_path)
  for (i in 1:3) write_once(inputs[[i]],files$path[i])
  files$md5 <- md5(files$path)
  record <- list(schema_version=1L,mode="reconstructed",files=files,evidence=evidence,
    runtime=list(R=R.version.string,R_home=R.home(),RNGkind=RNGkind(),rng_state=.Random.seed),
    limitations="Numerically audited seeded reconstruction; missing original matrices and bytes prevent full historical identity certification.")
  write_once(record,file.path(output_root,"input-record.rds")); record
}
import_inputs <- function(input_record,output_root) {
  check_input_mode(input_record$mode)
  files <- input_record$files
  if (any(md5(files$path)!=files$md5)) stop("Input hash mismatch before import")
  if (!all(input_record$evidence$checks$passed) || length(files$path)!=3L) stop("Unverified inputs")
  rows <- list()
  for (i in seq_len(nrow(files))) {
    inp <- readRDS(files$path[i]); land <- inp$landscape; comm <- files$identity[i]
    for (arr in c("spread","pairs","clustered","grid")) {
      dataset <- paste(comm,arr,sep="-"); tr <- inp$surveys[[arr]]; idx <- land$index[[arr]]; lat <- land$index$lattice
      pack <- function(ids) list(schema_version=1L,dataset_id=dataset,
        site_id=sprintf("%s-point%04d",comm,ids),species_id=land$species,
        environment=matrix(land$environment[ids],ncol=1L,dimnames=list(NULL,"environment")),
        xy=land$points[ids,,drop=FALSE],input_hashes=setNames(files$md5[i],comm))
      training <- pack(idx); training$y <- tr$binary$OTU
      prediction <- pack(lat)
      truth <- list(schema_version=1L,dataset_id=dataset,species_id=land$species,
                    training_site_id=training$site_id,prediction_site_id=prediction$site_id,
                    training_index=idx,lattice_index=lat,landscape=land,input_hashes=training$input_hashes)
      for (kind in c("training","prediction","truth")) {
        rel <- paste0(kind,"/",dataset,".rds"); path <- file.path(output_root,rel)
        write_once(get(kind),path)
        rows[[length(rows)+1L]] <- data.frame(identity=paste(dataset,kind,sep="-"),
          dataset_id=dataset,kind=kind,relative_path=rel,md5=md5(path),input_mode=input_record$mode)
      }
    }
  }
  m <- do.call(rbind,rows); write_once(m,file.path(output_root,"import-manifest.rds")); m
}
verify_runtime <- function(package,launcher) {
  output <- tempfile(fileext=".rds"); on.exit(unlink(output),add=TRUE)
  # Values are R literals passed as a single quoted command argument, not shell code.
  q <- function(x) encodeString(x,quote='"')
  code <- paste0('p <- ',q(package),'; library(p,character.only=TRUE); path <- find.package(p); ',
    'f <- list.files(path,recursive=TRUE,full.names=TRUE); ',
    'saveRDS(list(schema_version=1L,package=p,version=as.character(packageVersion(p)),',
    'R=R.version.string,R_home=R.home(),path=path,files=tools::md5sum(f)),',q(output),')')
  logs <- system2(launcher,c("-e",shQuote(code)),stdout=TRUE,stderr=TRUE)
  status <- attr(logs,"status"); if (is.null(status)) status <- 0L
  if (status!=0L || !file.exists(output)) stop("Runtime probe failed for ",package,": ",paste(logs,collapse="\n"))
  ans <- readRDS(output); ans$log <- logs; ans
}

# Reused checks from the committed planning audit; accepting inputs lets tests corrupt evidence.
audit_inputs <- function(inputs, repo_root, generator_md5=study_config()$generator_md5, tolerance=1e-12) {
  cfg <- study_config()
  if (!identical(tolerance,1e-12)) stop("Frozen audit tolerance must be 1e-12")
  generator_path <- file.path(repo_root,"dev/simstudy/spatial-design-sweep/generator.R")
  if (!identical(generator_md5,cfg$generator_md5) || !identical(md5(generator_path),cfg$generator_md5)) stop("generator hash mismatch")
  bundle_path <- file.path(repo_root,"vignettes/teaching-data/spatial-lesson.rds")
  if (!identical(md5(bundle_path),cfg$bundle_md5)) stop("Teaching evidence hash mismatch")
  b <- readRDS(bundle_path)
  source(generator_path,local=TRUE)
  stopifnot(length(inputs)==3L)
checks <- list()
replay <- list()
record <- function(community, check, n, difference = 0, exact = FALSE) {
  if (!is.finite(difference) || difference > if (exact) 0 else tolerance) stop("Audit mismatch: ", community, " ", check)
  checks[[length(checks) + 1L]] <<- data.frame(community, check, n,
    max_absolute_difference = difference, passed = TRUE, criterion = if (exact) "exact" else "absolute <= 1e-12")
}
for (r in 1:3) {
  inp <- inputs[[r]]
  expected_seeds <- c(sites=sweep_seed("spatial-design-sites",r), data=sweep_seed("spatial-design-data",r), detection=sweep_seed("spatial-design-detection",r))
  if (!identical(inp$settings$seeds,expected_seeds)) stop("seed mismatch")
  id <- sprintf("rep%02d", r)
  if (r == 1) for (nm in names(b$landscape)) {
    old <- b$landscape[[nm]]; new <- inp$landscape[[nm]]
    if (is.numeric(old)) {
      stopifnot(identical(dim(old), dim(new)))
      record(id, paste0("landscape_", nm), length(old), max(abs(old - new)))
    } else {
      stopifnot(identical(old, new))
      record(id, paste0("landscape_", nm), length(old), exact = TRUE)
    }
  }
  ar <- b$arrangements[b$arrangements$community == id, ]
  stopifnot(nrow(ar) == 400L)
  xy <- t(mapply(function(a, i) inp$surveys[[a]]$truth$xy[i, ], ar$arrangement, ar$site))
  stopifnot(identical(unname(as.matrix(ar[c("x", "y")])), unname(xy)))
  record(id, "survey_coordinates", length(xy), exact = TRUE)
  fm <- b$field_maps[b$field_maps$community == id & b$field_maps$source == "truth", ]
  stopifnot(nrow(fm) == 800L)
  values <- mapply(function(a, i, s) inp$surveys[[a]]$truth$field[i, s], fm$arrangement, fm$site, fm$species)
  record(id, "saved_two_species_fields", length(values), max(abs(values - fm$value)))
  oracle <- b$oracle[b$oracle$community == id, ]
  stopifnot(nrow(oracle) == 32L)
  null_rmse <- mapply(function(a, s) {
    f <- inp$surveys[[a]]$truth$field[, s]; sqrt(mean((f - mean(f))^2))
  }, oracle$arrangement, oracle$species)
  record(id, "all_species_zero_field_rmse", length(null_rmse), max(abs(null_rmse - oracle$zero_field_rmse)))
  stats <- b$statistics[b$statistics$replicate == r, ]
  stats <- stats[match(inp$statistics$arrangement, stats$arrangement), ]
  stopifnot(nrow(stats) == 4L)
  for (nm in names(stats)) {
    if (is.numeric(stats[[nm]])) {
      record(id, paste0("design_", nm), nrow(stats), max(abs(inp$statistics[[nm]] - stats[[nm]])))
    } else stopifnot(identical(inp$statistics[[nm]], stats[[nm]]))
  }
  ss <- b$fits$species[b$fits$species$community == id, ]
  stopifnot(nrow(ss) == 64L)
  occupied <- mapply(function(a, s) sum(inp$surveys[[a]]$truth$z[, s]), ss$arrangement, ss$species)
  detections <- mapply(function(a, arm, s) sum(inp$surveys[[a]][[arm]]$OTU[, s]), ss$arrangement, ss$arm, ss$species)
  mean_psi <- mapply(function(a, s) mean(inp$surveys[[a]]$truth$psi[, s]), ss$arrangement, ss$species)
  record(id, "occupied_counts_both_arms", length(occupied), max(abs(occupied - ss$occupied)), TRUE)
  record(id, "detection_counts_both_arms", length(detections), max(abs(detections - ss$detections)), TRUE)
  record(id, "mean_true_probability", length(mean_psi), max(abs(mean_psi - ss$truth)))

  # Replay the generator's data stream through the two environmental normal
  # vectors and species-field normal matrix, then recover its occupancy uniforms.
  N <- nrow(inp$landscape$points); S <- length(inp$landscape$species)
  set.seed(sweep_seed("spatial-design-data", r))
  invisible(rnorm(N)); invisible(rnorm(N)); invisible(rnorm(N * S))
  u <- matrix(runif(N * S), N, S)
  z <- matrix(as.integer(u < inp$landscape$psi), N, S)
  stopifnot(identical(z, unname(inp$landscape$z)))
  record(id, "replayed_uniforms_match_regenerated_z", length(z), exact = TRUE)
  gap <- min(abs(u - inp$landscape$psi))
  stopifnot(is.finite(gap), gap > tolerance)
  replay[[r]] <- data.frame(community = id, n = length(u),
    minimum_uniform_threshold_distance = gap, comparison_tolerance = tolerance)
  if (r == 1L) {
    z_saved_psi <- matrix(as.integer(u < b$landscape$psi), N, S)
    stopifnot(identical(z_saved_psi, unname(inp$landscape$z)))
    record(id, "saved_original_psi_and_replayed_uniforms_match_z", length(z), exact = TRUE)
  }
}
  list(checks=do.call(rbind,checks),replay=do.call(rbind,replay),generator_md5=md5(generator_path),bundle_md5=md5(bundle_path))
}
