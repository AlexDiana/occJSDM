# Write ROOT/inputs/training.rds for the sjSDM v0.1.0 rerun. The original
# run root is gone; the committed teaching bundle holds its training object
# verbatim (export-teaching.R copied it, verify-teaching.R checked identity).
# Written by R 4.5.0, as the original was, the file reproduces the recorded
# MD5; R 4.6.1 writes the same object with different bytes.
# From the repository root, with the pinned R 4.5.0 launcher:
#   sh dev/simstudy/spatial-recoverability-pilot/run-r \
#     -f dev/simstudy/jsdm-package-comparison/prepare-sjsdm-v010-input.R --args ROOT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, R.version$major == "4", R.version$minor == "5.0")
root <- args[1]
path <- file.path(root, "inputs", "training.rds")
stopifnot(!file.exists(path))
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
bundle <- readRDS("vignettes/teaching-data/jsdm-comparison.rds")
recorded <- "c37f293de37bd932bfb24d426537641a"
stopifnot(identical(unname(bundle$provenance$source_hashes[1]), recorded),
          grepl("/inputs/training.rds$", names(bundle$provenance$source_hashes)[1]))
saveRDS(bundle$training, path)
stopifnot(identical(unname(tools::md5sum(path)), recorded))
cat("Wrote", path, "with the recorded MD5", recorded, "\n")
