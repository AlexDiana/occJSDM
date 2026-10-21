# Cost probes cited in SHELVED.md. Not study attempts; nothing here feeds a decision rule. Run through run-r:
#   run-r -f probe-costs.R --args STUDY gpp OUT_TXT        Hmsc with the GPP method at several knot counts on a 900-site community
#   run-r -f probe-costs.R --args STUDY gp-slope OUT_TXT   occJSDM gp at 300 support points on two short schedules, to separate fixed from per-iteration cost
args <- commandArgs(trailingOnly = TRUE); study <- args[1]; mode <- args[2]; out <- args[3]
source("dev/simstudy/spatial-stage2-comparison/jobs.R")
log <- function(...) { line <- paste0(...); cat(line, "\n"); cat(line, "\n", file = out, append = TRUE) }
cat("", file = out)
cfg <- readRDS(file.path(study, "config.rds"))
if (mode == "gpp") {
  job <- job_table(study); job <- job[job$dataset_id == "t2-n900-rep1" & job$key == "Hmsc/spatial2", ]
  tr <- load_training(readRDS(job$input_file), job$dataset_id, job$input_md5); xy <- tr$xy; rownames(xy) <- tr$site_id
  log("Hmsc GPP probe on ", job$dataset_id, " (900 sites, 12 species), k-means knots, one chain, 5 warm-up and 5 retained draws, nf = 2")
  for (k in c(50L, 100L, 200L, 300L)) {
    set.seed(1); knots <- stats::kmeans(xy, centers = k, nstart = 5L, iter.max = 100L)$centers
    d <- as.matrix(dist(knots)); diag(d) <- Inf
    lvl <- Hmsc::setPriors(Hmsc::HmscRandomLevel(sData = xy, sMethod = "GPP", sKnot = knots), nfMin = 2L, nfMax = 2L)
    r <- tryCatch({
      m <- Hmsc::Hmsc(Y = tr$y, XData = data.frame(environment = as.numeric(scale(tr$environment))), XScale = FALSE, XFormula = ~., distr = "probit",
                      studyDesign = data.frame(site = factor(tr$site_id, levels = tr$site_id)), ranLevels = list(site = lvl))
      t0 <- Sys.time(); Hmsc::sampleMcmc(m, samples = 5L, transient = 5L, thin = 1L, nChains = 1L, nParallel = 1L, verbose = 0)
      sprintf("ran, %.1f seconds for 10 iterations", as.numeric(Sys.time() - t0, units = "secs"))
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    log("knots ", k, ", minimum knot distance ", signif(min(d), 3), ": ", r)
  }
} else if (mode == "gp-slope") {
  frozen <- list(support_points_900 = 300L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L)
  job <- job_table(study); job <- job[job$dataset_id == "t1-r0.08-n900-rep1" & job$key == "occJSDM/gp", ]
  tr <- load_training(readRDS(job$input_file), job$dataset_id, job$input_md5); basis <- fit_basis(tr, "gp")
  log("occJSDM gp probe on ", job$dataset_id, " (900 sites), 300 support points, one chain")
  secs <- c(); its <- c()
  for (sched in list(c(nburn = 5L, niter = 10L), c(nburn = 10L, niter = 20L))) {
    st <- attempt_settings(job, 1L, cfg, frozen, list(initial = c(nchain = 1L, sched, nthin = 1L), extended = NULL))
    t <- system.time(suppressMessages(capture.output(f <- fit_occjsdm(tr, basis, st))))[["elapsed"]]
    secs <- c(secs, t); its <- c(its, sum(sched)); log("iterations ", sum(sched), ": ", round(t, 1), " seconds")
  }
  slope <- diff(secs) / diff(its); log(sprintf("per-iteration cost %.1f seconds, fixed cost about %.1f seconds", slope, secs[1] - slope * its[1]))
}
