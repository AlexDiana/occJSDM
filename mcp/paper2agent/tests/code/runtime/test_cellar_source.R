.libPaths(c("r-runtime/bootstrap-library", .Library))
lock <- renv::lockfile_read("r-runtime/renv.lock")
record <- lock$Packages$occJSDM
archive <- file.path("r-runtime/renv/cellar/occJSDM", paste0("occJSDM_", record$RemoteSha, ".tar.gz"))
stopifnot(file.exists(archive))
cat("Pinned source archive is discoverable using renv's RemoteSha naming.\n")
