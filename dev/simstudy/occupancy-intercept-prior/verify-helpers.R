# Verification-only helpers for the occupancy-intercept prior study: control
# selection, control provenance, equivalence and pilot checks. run.R does not
# source or hash this file, so it can change without touching the fit identity.
# Source jobs.R first.

# Relocate a recorded control-fit path into the given archives root.
relocate_archive_path <- function(path,archives,default) {
  for(name in names(ARCHIVES)) {
    tag <- paste0('/',ARCHIVES[[name]],'/');at <- regexpr(tag,path,fixed=TRUE)
    if(at>0) return(file.path(archive_dir(archives,name),substring(path,at+nchar(tag))))
  }
  if(startsWith(path,'/')) stop('Control fit is outside the known archives: ',path)
  file.path(default,path)
}

# For each community, the schedule of the control arm's selected fit.
control_schedules <- function(phase,archives) {
  keys <- phase_keys(phase)
  if(phase %in% c('A1','B')) {
    pr11 <- archive_dir(archives,'pr11')
    selection <- readRDS(file.path(pr11,'selection.rds'))
    plan <- read.csv(file.path(pr11,'long-selection.csv'),stringsAsFactors=FALSE)
    i <- match(keys,selection$key);j <- match(keys,plan$key)
    if(anyNA(i) || anyNA(j)) stop('Recorded pr11 selection is missing phase ',phase,' keys')
    schedule <- selection$schedule[i]
    if(!all(schedule %in% c('initial','long')) || !identical(schedule=='long',as.logical(plan$run_long[j])))
      stop('pr11 selection.rds and long-selection.csv disagree for phase ',phase)
    result <- selection$result_file[i]
    if(!identical(basename(result),paste0(keys,'-result.rds')) || !identical(basename(dirname(result)),schedule))
      stop('pr11 selected result files disagree with the selected schedules')
    fit <- file.path(pr11,schedule,paste0(keys,'-fit.rds'));md5 <- rep(NA_character_,length(keys))
  } else if(phase=='A2') {
    amplitude <- archive_dir(archives,'amplitude')
    final <- file.path(amplitude,'robust-v1/summary-binary-final')
    fits <- read.csv(file.path(final,'fits.csv'),stringsAsFactors=FALSE)
    fits <- fits[fits$prior=='inverse_gamma',]
    selection <- read.csv(file.path(final,'selection.csv'),stringsAsFactors=FALSE)
    i <- match(keys,fits$key);j <- match(keys,selection$key)
    if(anyNA(i) || anyNA(j) || anyDuplicated(fits$key)) stop('Recorded spatial-amplitude selection is missing A2 keys')
    schedule <- fits$phase[i]
    if(!all(schedule %in% c('initial','long')) || !identical(schedule=='long',as.logical(selection$needs_long[j])))
      stop('spatial-amplitude fits.csv and selection.csv disagree')
    fit <- vapply(fits$fit_file[i],relocate_archive_path,character(1),archives=archives,default=amplitude,USE.NAMES=FALSE)
    if(!identical(basename(fit),paste0(keys,'-fit.rds')) || !identical(basename(dirname(fit)),schedule))
      stop('spatial-amplitude selected fit files disagree with the selected schedules')
    md5 <- fits$fit_md5[i]
  } else stop('Unknown phase: ',phase)
  data.frame(phase=phase,key=keys,schedule=schedule,control_fit=fit,control_fit_md5=md5,stringsAsFactors=FALSE)
}

installed_library_files <- function(pkg)
  c(file.path(pkg,'libs/occJSDM.so'),file.path(pkg,c('DESCRIPTION','NAMESPACE')),list.files(file.path(pkg,'R'),full.names=TRUE))

# Names and dimensions of every posterior output array, for structure checks.
posterior_layout <- function(results) {
  walk <- function(x,prefix) {
    if(is.list(x) && !is.data.frame(x)) return(do.call(rbind,c(list(NULL),lapply(names(x),function(n) walk(x[[n]],paste0(prefix,n,'/'))))))
    data.frame(path=sub('/$','',prefix),class=class(x)[1],dim=paste(dim(x) %||% length(x),collapse='x'),stringsAsFactors=FALSE)
  }
  walk(results,'')
}

# Largest absolute difference over all numeric leaves; Inf when structures differ.
max_abs_difference <- function(a,b) {
  if(is.list(a) || is.list(b)) {
    if(!is.list(a) || !is.list(b) || !identical(names(a),names(b)) || length(a)!=length(b)) return(Inf)
    if(!length(a)) return(0)
    return(max(vapply(seq_along(a),function(i) max_abs_difference(a[[i]],b[[i]]),numeric(1))))
  }
  if(!identical(dim(a),dim(b)) || length(a)!=length(b)) return(Inf)
  if(is.numeric(a) && is.numeric(b)) {
    if(!identical(is.na(a),is.na(b))) return(Inf)
    ok <- !is.na(a);if(!any(ok)) return(0)
    inf <- is.infinite(a[ok]) | is.infinite(b[ok])
    if(any(inf) && !identical(a[ok][inf],b[ok][inf])) return(Inf)
    d <- abs(a[ok][!inf]-b[ok][!inf]);return(if(length(d)) max(d) else 0)
  }
  if(identical(a,b)) 0 else Inf
}
count_numeric <- function(x) if(is.list(x)) sum(vapply(x,count_numeric,numeric(1))) else if(is.numeric(x)) length(x) else 0

# Record the installed build that run.R will insist on (see check_library_fingerprint).
write_library_fingerprint <- function(study,csv) {
  revision <- readLines(file.path(study,'source-revision.txt'))
  stopifnot(length(revision)==1L,grepl('^[0-9a-f]{40}$',revision))
  h <- hash_files(library_fingerprint_files(study),study)
  utils::write.csv(data.frame(file=names(h),md5=unname(h),revision=revision),csv,row.names=FALSE)
  invisible(csv)
}

# Rebuild a saved fit's resume metadata through the same constructor run.R uses.
saved_fit_metadata <- function(saved) {
  fields <- names(formals(fit_metadata))
  absent <- setdiff(fields,names(saved))
  if(length(absent)) stop('Saved fit lacks metadata fields: ',paste(absent,collapse=', '))
  do.call(fit_metadata,saved[fields])
}

# End a verification script: exit status 1 unless every check is TRUE.
finish_gate <- function(pass,what) {
  if(!length(pass) || anyNA(pass) || !all(pass)) {
    bad <- if(is.null(names(pass))) which(is.na(pass) | !pass) else names(pass)[is.na(pass) | !pass]
    cat(what,'FAILED:',paste(bad,collapse=', '),'\n')
    quit(save='no',status=1)
  }
  cat('All',what,'checks pass.\n')
  invisible(TRUE)
}
