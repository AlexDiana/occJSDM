# Collect evidence about the BLAS this R uses: where it is, whether the file
# is a symlink, what it links to, whether it defines the BLAS symbols itself,
# its md5, and, on macOS, the installer receipt that lists it. This script
# makes no verdict and sets no variable. A person reads the evidence and
# declares OCCJSDM_REFERENCE_BLAS in the launching environment, with the basis
# of the declaration in OCCJSDM_REFERENCE_BLAS_BASIS; the launcher records both.
# Spec: "Tests", the declaration.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
out <- args[[1]]
darwin <- Sys.info()[["sysname"]] == "Darwin"
blas <- unname(extSoftVersion()[["BLAS"]])
lapack <- La_library()
info <- Sys.readlink(blas)
linkage <- if (darwin) system2("otool", c("-L", shQuote(blas)), stdout = TRUE) else system2("ldd", shQuote(blas), stdout = TRUE)
symbols <- if (darwin) system2("nm", c("-gU", shQuote(blas)), stdout = TRUE) else system2("nm", c("-D", shQuote(blas)), stdout = TRUE)
receipt <- NULL
if (darwin) {
  pkgs <- system2("pkgutil", "--pkgs", stdout = TRUE)
  fw <- grep("^org[.]R-project[.]R[.]fw[.]pkg$", pkgs, value = TRUE)
  if (length(fw) == 1L) {
    receipt <- list(package = fw,
                    info = system2("pkgutil", c("--pkg-info", fw), stdout = TRUE),
                    lists_blas = any(grepl(basename(blas), system2("pkgutil", c("--files", fw), stdout = TRUE), fixed = TRUE)))
  }
}
facts <- list(blas = blas, lapack = lapack, blas_is_symlink = nzchar(info), symlink_target = info,
              file_info = unclass(file.info(blas))[c("size", "mtime", "uid", "gid")],
              linkage = linkage,
              names_other_blas = any(grepl("accelerate|veclib|openblas|mkl|blis|atlas", linkage, ignore.case = TRUE)),
              defines_dgemm = any(grepl(if (darwin) " T _dgemm_$" else " T dgemm_$", symbols)),
              blas_md5 = unname(tools::md5sum(blas)), installer_receipt = receipt,
              r_version = R.version.string, checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(facts, out)
cat("BLAS evidence written to", out, "\n")
cat("  library:", blas, "\n  symlink:", facts$blas_is_symlink, "\n  links to another BLAS by name:", facts$names_other_blas,
    "\n  defines dgemm_ itself:", facts$defines_dgemm, "\n  md5:", facts$blas_md5,
    "\n  installer receipt lists it:", if (is.null(receipt)) "not available" else receipt$lists_blas, "\n")
cat("No verdict is made here. Declare OCCJSDM_REFERENCE_BLAS yourself, with the basis in OCCJSDM_REFERENCE_BLAS_BASIS, if the evidence and the installation's provenance warrant it.\n")
