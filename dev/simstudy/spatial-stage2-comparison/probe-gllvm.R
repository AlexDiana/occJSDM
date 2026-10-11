# Fits gllvm/spatial2 on the first n training sites of a pilot dataset with a small iteration cap, to measure how time and memory scale
# with the site count. Not a study attempt. Run through probe-gllvm.sh, which watches the process footprint and enforces a cap.
#   run-r -f probe-gllvm.R --args STUDY DATASET_ID N MAXIT OUT_CSV [VARIANT]
# VARIANT is "unstructured" (the upstream adapter, the default) or a gllvm variational covariance structure such as "bdNN", which is a
# probe of an alternative the study has NOT adopted: it calls gllvm directly, because the adapter asserts the unstructured form.
args <- commandArgs(trailingOnly = TRUE); study <- args[1]; dataset <- args[2]; n <- as.integer(args[3]); maxit <- as.integer(args[4]); out <- args[5]; variant <- if (length(args) >= 6) args[6] else "unstructured"
source("dev/simstudy/spatial-stage2-comparison/jobs.R")
job <- job_table(study); job <- job[job$dataset_id == dataset & job$key == "gllvm/spatial2", ]
cfg <- readRDS(file.path(study, "config.rds"))
tr <- slice_training(load_training(readRDS(job$input_file), job$dataset_id, job$input_md5), n, seq_along(readRDS(job$input_file)$species$species))
basis <- fit_basis(tr, "spatial2"); settings <- attempt_settings(job, 1L, cfg, cfg$frozen); settings$gllvm$maxit <- maxit
if (variant == "unstructured") {
  t <- system.time(fit <- fit_gllvm(tr, basis, settings))[["elapsed"]]
} else {
  X <- occ_design(tr, basis)
  t <- system.time(fit <- gllvm::gllvm(y = tr$y, X = X, family = binomial("logit"), link = "logit", num.lv = 2L, method = "VA", row.eff = FALSE, scale.X = FALSE, sd.errors = FALSE,
    seed = settings$seed, control = list(maxit = maxit, max.iter = maxit), control.va = list(Lambda.struc = variant, NN = 10), control.start = list(starting.val = "random", n.init = 1L, jitter.var = .1),
    lvCor = ~corExp(1 | site), studyDesign = data.frame(site = factor(tr$site_id, levels = tr$site_id)), distLV = tr$xy))[["elapsed"]]
}
conv <- if (variant == "unstructured") isTRUE(fit$native$convergence) else isTRUE(fit$convergence)
write.csv(data.frame(dataset = dataset, n = n, maxit = maxit, variant = variant, fit_seconds = t, converged = conv, stringsAsFactors = FALSE), out, row.names = FALSE)
