# Independent verification for stage 2. The paired intervals and the development reading are recomputed from the saved
# per-species scores with code written without reference to score.R. Run through run-r:
#   run-r -f verify.R --args --study=DIR
# Exits non-zero on any failed check. A failure is a defect to diagnose, not a result.
D2 <- "dev/simstudy/spatial-stage2-comparison"
source(file.path(D2, "upstream", "adapter-contract.R")); source(file.path(D2, "data.R"))

check <- function(name, passed, detail = "") data.frame(check = name, passed = isTRUE(passed), detail = as.character(detail), stringsAsFactors = FALSE)

mean_finite <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }

# One headroom share per dataset and arm: species averaged within each gate occupancy level, then the levels averaged.
dataset_shares <- function(species, gate) {
  s <- species[species$occupancy %in% gate, ]
  level <- aggregate(headroom_share ~ dataset_id + key + occupancy, s, mean_finite, na.action = na.pass)
  share <- aggregate(headroom_share ~ dataset_id + key, level, mean_finite, na.action = na.pass)
  status <- aggregate(status ~ dataset_id + key, s, function(x) x[1])
  out <- merge(share, status, by = c("dataset_id", "key")); out$cell <- sub("-rep[0-9]+$", "", out$dataset_id); out
}

verify_pairs <- function(species, pairs, cfg) {
  d <- dataset_shares(species, cfg$gate_occupancy); out <- list()
  for (cell in unique(pairs$cell)) for (cmp in cfg$comparisons) {
    a <- d[d$cell == cell & d$key == cmp[1], c("dataset_id", "headroom_share", "status")]; b <- d[d$cell == cell & d$key == cmp[2], c("dataset_id", "headroom_share", "status")]
    ids <- intersect(a$dataset_id, b$dataset_id); a <- a[match(ids, a$dataset_id), ]; b <- b[match(ids, b$dataset_id), ]
    keep <- a$status == "pass" & b$status == "pass" & is.finite(a$headroom_share) & is.finite(b$headroom_share)
    diffs <- b$headroom_share[keep] - a$headroom_share[keep]; k <- length(diffs)
    est <- if (k) mean(diffs) else NA_real_; half <- if (k > 1) qt(0.975, k - 1) * sd(diffs) / sqrt(k) else NA_real_
    verdict <- if (k < cfg$min_paired) "inconclusive" else if (est >= cfg$delta && est - half > 0) "b clearly better" else
      if (est <= -cfg$delta && est + half < 0) "a clearly better" else "no clear difference"
    row <- pairs[pairs$cell == cell & pairs$a == cmp[1] & pairs$b == cmp[2], ]
    same <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-9)) || (is.na(x) && is.na(y))
    ok <- nrow(row) == 1L && row$paired == k && row$excluded == length(ids) - k && same(row$estimate, est) &&
      same(row$lower, if (is.na(half)) NA else est - half) && same(row$upper, if (is.na(half)) NA else est + half) && identical(row$verdict, verdict)
    out[[length(out) + 1L]] <- check(sprintf("pair %s %s vs %s", cell, cmp[1], cmp[2]), ok,
      sprintf("recomputed n=%d excluded=%d estimate=%.4f verdict=%s", k, length(ids) - k, est, verdict))
  }
  do.call(rbind, out)
}

# The development reading, written out from the spec's decision rule at the 900-site cells.
verify_reading <- function(species, pairs, stored, cfg) {
  d <- dataset_shares(species, cfg$gate_occupancy)
  medians <- aggregate(headroom_share ~ cell + key, d[is.finite(d$headroom_share), ], median)
  max_share <- max(medians$headroom_share)
  big <- pairs[grepl("n900$", pairs$cell), ]
  versus_gp <- big[big$a == "occJSDM/gp", ]
  truth1 <- versus_gp[startsWith(versus_gp$cell, "t1"), "verdict"]; truth2 <- versus_gp[startsWith(versus_gp$cell, "t2"), "verdict"]
  factors_win_1 <- "b clearly better" %in% truth1; gp_wins_2 <- "a clearly better" %in% truth2
  quiet <- all(pairs$verdict %in% c("no clear difference", "inconclusive"))
  reading <- if (factors_win_1 && !gp_wins_2) "shared factors" else if (!factors_win_1 && gp_wins_2) "per-species ranges" else
    if (factors_win_1 && gp_wins_2) "trade-off" else if (quiet && max_share < 0.5) "representation not the bottleneck" else
    if (quiet) "no clear difference" else "no pre-specified branch"
  dimension <- any(big$a == "Hmsc/spatial6" & big$b == "Hmsc/spatial2" & startsWith(big$cell, "t2") & big$verdict == "a clearly better")
  method <- any(big$a == "Hmsc/spatial2" & big$b == "gllvm/spatial2" & big$verdict != "no clear difference" & big$verdict != "inconclusive")
  rbind(check("development reading", identical(stored$reading, reading), sprintf("recomputed '%s', stored '%s'", reading, stored$reading)),
        check("largest median headroom share", isTRUE(all.equal(as.numeric(stored$max_share), max_share, tolerance = 1e-9)), sprintf("recomputed %.4f", max_share)),
        check("factor-dimension qualifier", identical(isTRUE(stored$factor_dimension), dimension), dimension),
        check("fitting-method qualifier", identical(isTRUE(stored$fitting_method), method), method))
}

# Status recomputed from the saved attempt summaries alone, following the spec's finite-budget rules.
recompute_status <- function(summaries, kind, cfg) {
  ids <- vapply(summaries, function(x) x$attempt_id, numeric(1)); done <- vapply(summaries, function(x) identical(x$status, "complete"), logical(1))
  good <- done & vapply(summaries, function(x) isTRUE(x$acceptable), logical(1))
  if (kind == "bayesian") return(if (any(good)) "pass" else if (length(ids) >= 2L && any(done)) "unresolved" else "no fit")
  if (!any(good)) return("no fit")
  obj <- vapply(summaries[good], function(x) x$objective, numeric(1)); best <- summaries[good][[which.min(obj)]]
  agree <- vapply(summaries[good][-which.min(obj)], function(x) abs(x$objective - best$objective) <= cfg$diagnostics$objective_agreement &&
                    max(abs(x$training_probability - best$training_probability)) <= cfg$diagnostics$probability_agreement, logical(1))
  if (any(agree)) "pass" else "unstable"
}

verify_study <- function(study, results_dir = file.path(D2, "results")) {
  cfg <- readRDS(file.path(study, "config.rds")); out <- list(); add <- function(x) out[[length(out) + 1L]] <<- x
  ds <- read.csv(file.path(study, "datasets.csv"), stringsAsFactors = FALSE)
  man <- read.csv(file.path(cfg$pilot_study, "input-manifest.csv"), stringsAsFactors = FALSE)
  add(check("40 datasets, ten per cell", nrow(ds) == 40L && all(table(ds$cell) == 10L), paste(table(ds$cell), collapse = " ")))
  add(check("recorded input hashes equal the pilot manifest", all(ds$input_md5 == man$md5[match(basename(ds$input_file), man$file)])))
  add(check("input files on disk hash to the manifest", all(unname(tools::md5sum(ds$input_file)) == ds$input_md5)))
  add(check("oracle files on disk hash to the recorded values", all(unname(tools::md5sum(ds$oracle_file)) == ds$oracle_md5)))
  resolved <- file.exists(file.path(cfg$pilot_study, "jobs-resolution", paste0(ds$cell, "-rep", ds$replicate, "-uniform.rds")))
  add(check("oracle reference is the resolution refit where one exists", all(grepl("jobs-resolution", ds$oracle_file) == resolved), sprintf("%d resolution refits", sum(resolved))))
  add(check("config thresholds equal the frozen protocol", { p <- source(file.path(D2, "config.R"), local = TRUE); f <- stage2_config()
    identical(cfg$delta, f$delta) && identical(cfg$mcmc, f$mcmc) && identical(cfg$diagnostics, f$diagnostics) && identical(cfg$limits, f$limits) && identical(cfg$arms, f$arms) }))
  # Counts
  status_all <- list.files(file.path(study, "status"), pattern = "\\.json$"); status_main <- status_all[!grepl("--timing", status_all)]
  add(check("200 status files after excluding timing jobs", length(status_main) == 200L, length(status_main)))
  add(check("no error files", length(list.files(file.path(study, "status"), pattern = "\\.error$")) == 0L))
  att <- list.files(file.path(study, "jobs", "attempts"), pattern = "\\.rds$"); att_main <- att[!grepl("--timing", att)]
  per_job <- table(sub("-[0-9]+\\.rds$", "", att_main)); kind <- ifelse(grepl("--gllvm-", names(per_job)), "point", "bayesian")
  add(check("at most 400 study attempts", length(att_main) <= cfg$limits$attempts, length(att_main)))
  add(check("at most two attempts per Bayesian job and six per point job", all(per_job[kind == "bayesian"] <= 2L) && all(per_job[kind == "point"] <= 6L)))
  # Frozen settings between the timing step and the first study attempt
  if (!is.null(cfg$frozen_from)) {
    step <- cfg$frozen_from$step; timing_file <- file.path(study, sprintf("timing-step%d.json", step))
    t_json <- file.mtime(timing_file); first_attempt <- min(file.mtime(file.path(study, "jobs", "attempts", att_main)))
    add(check("frozen settings set after the timing step and before the first study attempt", cfg$frozen_at > t_json && cfg$frozen_at < first_attempt,
              sprintf("timing %s, frozen %s, first attempt %s", format(t_json), format(cfg$frozen_at), format(first_attempt))))
    add(check("config written before every attempt", file.mtime(file.path(study, "config.rds")) < min(file.mtime(file.path(study, "jobs", "attempts", att)))))
    add(check("frozen step passed the 48-hour rule", isTRUE(jsonlite::fromJSON(timing_file)$passes)))
  } else add(check("frozen settings recorded", FALSE, "config has no frozen_from"))
  # Statuses recomputed from summaries, and every prediction validated
  jobs <- read.csv(file.path(study, "scores", "jobs.csv"), stringsAsFactors = FALSE)
  bad_status <- character(); bad_pred <- character(); species_of <- list()
  for (i in seq_len(nrow(jobs))) {
    id <- jobs$job_id[i]; k <- if (grepl("--gllvm-", id)) "point" else "bayesian"
    files <- list.files(file.path(study, "jobs", "summaries"), pattern = paste0("^", id, "-[0-9]+\\.rds$"), full.names = TRUE)
    st <- recompute_status(lapply(files, readRDS), k, cfg)
    stored <- jsonlite::fromJSON(file.path(study, "status", paste0(id, ".json")))$status
    if (!(identical(st, stored) || (stored == "unsupported" && st != "no fit"))) bad_status <- c(bad_status, id)
    pf <- file.path(study, "jobs", "predictions", paste0(id, ".rds"))
    if (file.exists(pf)) {
      ds_id <- jobs$dataset_id[i]
      if (is.null(species_of[[ds_id]])) species_of[[ds_id]] <- readRDS(ds$input_file[ds$dataset_id == ds_id])$species$species
      ok <- tryCatch(validate_prediction(readRDS(pf)$prediction, list(site_id = lattice_ids(1600L), species_id = species_of[[ds_id]])), error = function(e) FALSE)
      if (!isTRUE(ok)) bad_pred <- c(bad_pred, id)
    } else if (jobs$status[i] %in% c("pass", "unresolved", "unstable")) bad_pred <- c(bad_pred, id)
  }
  add(check("statuses agree with a recomputation from the attempt summaries", length(bad_status) == 0L, paste(head(bad_status, 5), collapse = ",")))
  add(check("every selected fit has a prediction that passes validate_prediction", length(bad_pred) == 0L, paste(head(bad_pred, 5), collapse = ",")))
  # Scores
  species <- read.csv(file.path(study, "scores", "species.csv"), stringsAsFactors = FALSE); pairs <- read.csv(file.path(study, "scores", "pairs.csv"), stringsAsFactors = FALSE)
  add(check("species.csv has 200 x 12 rows", nrow(species) == 2400L, nrow(species)))
  add(check("aggregate.csv has 200 rows", nrow(read.csv(file.path(study, "scores", "aggregate.csv"))) == 200L))
  add(check("pairs.csv has 4 cells x 5 comparisons", nrow(pairs) == 20L))
  add(check("every status is one of the five allowed", all(species$status %in% c("pass", "unresolved", "unstable", "unsupported", "no fit"))))
  add(check("oracle status carried into every score row", all(species$oracle_status %in% c("pass", "unresolved"))))
  p <- verify_pairs(species, pairs, cfg); for (i in seq_len(nrow(p))) add(p[i, ])
  r <- verify_reading(species, pairs, jsonlite::fromJSON(file.path(study, "scores", "reading.json")), cfg); for (i in seq_len(nrow(r))) add(r[i, ])
  res <- do.call(rbind, out)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE); write.csv(res, file.path(results_dir, "verification.csv"), row.names = FALSE)
  res
}

if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  study <- sub("^--study=", "", grep("^--study=", args, value = TRUE))
  res <- verify_study(study); print(res[, c("check", "passed")], right = FALSE)
  if (all(res$passed)) cat("Verification passed.\n") else { cat("Verification FAILED:", sum(!res$passed), "check(s)\n"); quit(status = 1L) }
}
