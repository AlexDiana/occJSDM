# Restore into a new project, never reusing the development package library.
# Usage: Rscript --vanilla restore-runtime.R TARGET [--use-local-binaries]
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Supply a new absolute target project directory")
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
source_project <- dirname(normalizePath(sub("^--file=","",script_arg[[1]]),mustWork=TRUE))
target <- args[[1]]
if (!grepl("^/",target)) stop("Target must be absolute")
if (dir.exists(target) && length(list.files(target,all.files=TRUE,no..=TRUE))) stop("Target must be new or empty")
dir.create(target,recursive=TRUE,showWarnings=FALSE)
target <- normalizePath(target,mustWork=TRUE)
use_binaries <- "--use-local-binaries" %in% args
if (use_binaries && !(Sys.info()[["sysname"]]=="Darwin" && R.version$arch=="aarch64" && startsWith(as.character(getRversion()),"4.5."))) stop("The recorded binary cache is only for macOS arm64 R 4.5")
for (f in c("activate.R","renv.lock","source-provenance.json")) stopifnot(file.copy(file.path(source_project,f),file.path(target,f)))
dir.create(file.path(target,"renv"),showWarnings=FALSE)
for (f in c("activate.R","settings.json")) stopifnot(file.copy(file.path(source_project,"renv",f),file.path(target,"renv",f)))
writeLines('source("activate.R")',file.path(target,".Rprofile"))
source_cellar <- file.path(source_project,"renv","cellar")
target_cellar <- file.path(target,"renv","cellar")
dir.create(target_cellar,recursive=TRUE,showWarnings=FALSE)
if (use_binaries) {
  for (p in list.files(source_cellar,full.names=TRUE)) stopifnot(file.copy(p,target_cellar,recursive=TRUE))
} else {
  stopifnot(file.copy(file.path(source_cellar,"occJSDM"),target_cellar,recursive=TRUE))
}
boot <- file.path(target,"bootstrap-library")
dir.create(boot,recursive=TRUE,showWarnings=FALSE)
if (use_binaries) {
  stopifnot(file.copy(file.path(source_project,"bootstrap-library","renv"),boot,recursive=TRUE))
} else {
  install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.2.3.tar.gz",lib=boot,repos=NULL,type="source")
}
Sys.setenv(P2A_R_PROJECT=target,RENV_PROJECT=target,RENV_PATHS_ROOT=file.path(target,"state"),RENV_PATHS_LIBRARY=file.path(target,"library"),RENV_PATHS_CACHE=file.path(target,"cache"),RENV_PATHS_SANDBOX=file.path(target,"sandbox"),RENV_CONFIG_CACHE_ENABLED="FALSE")
.libPaths(c(boot,.Library))
source_lock <- renv::lockfile_read(file.path(target,"renv.lock"))
pinned_name <- paste0("occJSDM_",source_lock$Packages$occJSDM$RemoteSha,".tar.gz")
versioned_source <- file.path(target_cellar,"occJSDM","occJSDM_0.1.0.tar.gz")
if (!file.exists(file.path(target_cellar,"occJSDM",pinned_name))) stopifnot(file.copy(versioned_source,file.path(target_cellar,"occJSDM",pinned_name)))
Sys.setenv(RENV_PATHS_CELLAR=target_cellar)
stopifnot(identical(as.character(packageVersion("renv")),"1.2.3"),startsWith(find.package("renv"),target))
renv::load(project=target)
renv::settings$external.libraries(character(),project=target)
renv::settings$use.cache(FALSE,project=target)
if (use_binaries) {
  # Metadata comes from the tested lock; all install payloads still require cellar files.
  local_repo <- file.path(target,"offline-repository")
  contrib <- file.path(local_repo,"src","contrib")
  dir.create(contrib,recursive=TRUE,showWarnings=FALSE)
  lock_records <- renv::lockfile_read(file.path(target,"renv.lock"))$Packages
  fields <- c("Package","Version","Depends","Imports","LinkingTo","Suggests","NeedsCompilation","Priority")
  index <- matrix(NA_character_,nrow=length(lock_records),ncol=length(fields),dimnames=list(names(lock_records),fields))
  for (p in names(lock_records)) for (f in fields) if (!is.null(lock_records[[p]][[f]])) index[p,f] <- paste(unlist(lock_records[[p]][[f]]),collapse=", ")
  write.dcf(index,file=file.path(contrib,"PACKAGES"))
  saveRDS(index,file=file.path(contrib,"PACKAGES.rds"))
  Sys.setenv(RENV_CONFIG_REPOS_OVERRIDE=paste0("file://",local_repo))
}
renv::restore(project=target,library=renv::paths$library(project=target),lockfile=file.path(target,"renv.lock"),prompt=FALSE)
source(file.path(target,"activate.R"))
library(occJSDM)
library(jsonlite)
library(digest)
lock <- renv::lockfile_read(file.path(target,"renv.lock"))
ip <- installed.packages(lib.loc=renv::paths$library(project=target))
stopifnot(all(vapply(rownames(ip),function(p)startsWith(find.package(p),target),logical(1))))
stopifnot(all(vapply(rownames(ip),function(p)packageVersion(p)==package_version(lock$Packages[[p]]$Version),logical(1))))
provenance <- read_json(file.path(target,"source-provenance.json"),simplifyVector=TRUE)
stopifnot(identical(digest(file.path(target,provenance$cellar_path),algo="sha256",file=TRUE),provenance$cellar_sha256))
write_json(list(status="passed",r_version=R.version.string,project=target,library_paths=.libPaths(),occjsdm_origin=find.package("occJSDM"),packages=as.list(setNames(ip[,"Version"],ip[,"Package"])),restoration=if(use_binaries)"documented macOS dependency binary cellar with occJSDM rebuilt from source" else "repository source dependencies with retained occJSDM source cellar"),file.path(target,"restore-verification.json"),auto_unbox=TRUE,pretty=TRUE)
cat("Restore and isolated origin/version checks passed:",target,"\n")
