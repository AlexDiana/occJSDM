# Independent verification of the recoverability pilot. Recomputes every gate decision from
# species.csv with its own code, and checks seeds, inputs, status files and conditioning logs.
# Run after score:  run-r -f verify.R --args --study=DIR
args <- commandArgs(trailingOnly = TRUE)
study <- sub("^--study=", "", grep("^--study=", args, value = TRUE))
stopifnot(length(study) == 1L, dir.exists(study))
checks <- list()
add <- function(name, ok, detail = "") checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(ok), detail = as.character(detail), stringsAsFactors = FALSE)

cfg <- readRDS(file.path(study, "config.rds")); th <- cfg$thresholds
species <- read.csv(file.path(study, "scores", "species.csv"), stringsAsFactors = FALSE)
gate <- read.csv(file.path(study, "scores", "gate.csv"), stringsAsFactors = FALSE)
man <- read.csv(file.path(study, "input-manifest.csv"), stringsAsFactors = FALSE)

add("species rows equal communities x arrangements x 12", nrow(species) == nrow(man) * 2L * 12L, nrow(species))
add("input hashes unchanged since prepare", all(unname(tools::md5sum(file.path(study, "inputs", man$file))) == man$md5))
add("community seeds distinct", !anyDuplicated(c(man$seed_sites, man$seed_data)) > 0)
add("config predates every input", file.mtime(file.path(study, "config.rds")) < min(file.mtime(file.path(study, "inputs", man$file))))
add("thresholds in config are the spec's", identical(th, list(pattern_median = 0.6, pattern_floor = 0.4, pattern_count = 8L,
                                                               slope_low = 0.5, slope_high = 1.5, headroom_pp = 1.0, communities = 10L)))

# Independent gate recomputation, written without reference to score.R.
u <- species[species$arrangement == "uniform" & species$occupancy %in% c(0.25, 0.5, 0.75), ]
u$group <- ifelse(u$truth == "t2", u$class, "all")
per_level <- aggregate(u[c("contrast_correlation", "field_slope", "headroom_pp")], u[c("cell", "group", "replicate", "occupancy")], mean)
per_comm <- aggregate(per_level[c("contrast_correlation", "field_slope", "headroom_pp")], per_level[c("cell", "group", "replicate")], mean)
unres <- aggregate(list(unresolved = u$status == "unresolved"), u[c("cell", "group", "replicate")], any)
per_comm <- merge(per_comm, unres)
decide <- function(d) {
  if (any(d$unresolved)) return("unresolved")
  p <- d$contrast_correlation; s <- d$field_slope; h <- d$headroom_pp
  ci <- mean(h) + c(-1, 1) * qt(0.975, length(h) - 1) * sd(h) / sqrt(length(h))
  ok <- median(p) >= th$pattern_median && sum(p >= th$pattern_floor) >= th$pattern_count &&
    median(s) >= th$slope_low && median(s) <= th$slope_high && mean(h) >= th$headroom_pp && ci[1] > 0
  if (ok) "pass" else "fail"
}
recomputed <- do.call(rbind, lapply(split(per_comm, per_comm[c("cell", "group")], drop = TRUE),
                                    function(d) data.frame(cell = d$cell[1], group = d$group[1], decision = decide(d), communities = nrow(d))))
cmp <- merge(gate[c("cell", "group", "decision", "communities")], recomputed, by = c("cell", "group"), suffixes = c("_scored", "_recomputed"))
add("every gate row recomputed", nrow(cmp) == nrow(gate) && nrow(cmp) == nrow(recomputed), nrow(cmp))
add("gate decisions recomputed identically", all(cmp$decision_scored == cmp$decision_recomputed), paste(sum(cmp$decision_scored != cmp$decision_recomputed), "mismatches"))
add("ten communities behind every gate row", all(cmp$communities_scored == cfg$communities & cmp$communities_recomputed == cfg$communities))

# Map reading recomputed from the gate table.
map <- jsonlite::fromJSON(file.path(study, "scores", "map.json"))
t1n <- gate$n[gate$truth == "t1" & gate$range %in% 0.08 & gate$decision == "pass"]
t2n <- gate$n[gate$truth == "t2" & gate$class %in% "medium" & gate$decision == "pass"]
branch <- if (length(t1n) && length(t2n)) "both" else if (length(t1n)) "truth1_only" else if (length(t2n)) "truth2_only" else "stop"
add("map branch recomputed identically", identical(map$branch, branch), paste(map$branch, "vs", branch))

cond <- read.csv(file.path(study, "scores", "conditionals.csv"), stringsAsFactors = FALSE)
add("lattice conditional reconstruction within 1e-5", all(cond$error <= 1e-5), sprintf("max %.2e", max(cond$error)))
gen <- do.call(rbind, lapply(list.files(file.path(study, "logs"), pattern = "^generator-conditioning", full.names = TRUE), read.csv, stringsAsFactors = FALSE))
add("generator factor errors within 1e-5", all(gen$error <= 1e-5), sprintf("max %.2e", max(gen$error)))
add("generator conditioning logged for every community", length(unique(paste(gen$cell, gen$replicate))) == nrow(man))
fac <- read.csv(file.path(study, "scores", "factors.csv"), stringsAsFactors = FALSE)
add("prior factor errors within 1e-5", all(fac$error <= 1e-5), sprintf("max %.2e", max(fac$error)))

status_files <- grep("\\.resolution\\.json$", list.files(file.path(study, "status"), pattern = "\\.json$"), value = TRUE, invert = TRUE)
add("every job has a status file", length(status_files) == nrow(man) * 2L, length(status_files))
add("no job error files", length(list.files(file.path(study, "status"), pattern = "\\.error$")) == 0L)
jobs <- read.csv(file.path(study, "scores", "jobs.csv"), stringsAsFactors = FALSE)
add("every job scored", nrow(jobs) == nrow(man) * 2L, nrow(jobs))
add("unresolved status carried from jobs to species", all((jobs$status == "unresolved") == (tapply(species$status == "unresolved", species$job <- paste0(species$cell, "-rep", species$replicate, "-", species$arrangement), any)[jobs$job])))

# Prescribed resolution extension, if one was run (resolve.R).
if (file.exists(file.path(study, "resolution.json"))) {
  p <- jsonlite::fromJSON(file.path(study, "resolution.json"), simplifyVector = TRUE)
  prescribed <- sort(c(as.character(unlist(p$t1_jobs)), as.character(unlist(p$t2_jobs))))
  present <- sort(sub("\\.rds$", "", list.files(file.path(study, "jobs-resolution"), pattern = "\\.rds$")))
  add("resolution refits equal the prescription", identical(prescribed, present), sprintf("%d prescribed, %d present", length(prescribed), length(present)))
  add("resolution prescribed after first scoring and before refits", {
    created <- as.POSIXct(p$created); files <- file.path(study, "jobs-resolution", paste0(present, ".rds"))
    length(files) == 0L || all(file.mtime(files) > created) })
  add("no resolution error files", length(list.files(file.path(study, "status"), pattern = "\\.resolution\\.error$")) == 0L)
  carried_ok <- TRUE; refit_ok <- TRUE; multiplier_ok <- TRUE
  for (j in intersect(present, as.character(unlist(p$t2_jobs)))) {
    r <- readRDS(file.path(study, "jobs-resolution", paste0(j, ".rds"))); o <- readRDS(file.path(study, "jobs", paste0(j, ".rds")))
    refit <- r$resolution$refit_species; keep <- !(o$scores$species %in% refit)
    carried_ok <- carried_ok && isTRUE(all.equal(r$scores[keep, ], o$scores[keep, ], check.attributes = FALSE))
    refit_ok <- refit_ok && setequal(refit, as.character(unlist(p$t2_refits[[j]]))) && all(r$phases[o$scores$species %in% refit] == "resolution")
    multiplier_ok <- multiplier_ok && r$resolution$budget$niter == cfg$budget$ext_niter * p$multiplier && r$resolution$budget$nburn == cfg$budget$ext_nburn * p$multiplier
  }
  for (j in intersect(present, as.character(unlist(p$t1_jobs)))) {
    r <- readRDS(file.path(study, "jobs-resolution", paste0(j, ".rds")))
    refit_ok <- refit_ok && identical(r$phases, "resolution") && length(r$resolution$refit_species) == 12L
    multiplier_ok <- multiplier_ok && r$resolution$budget$niter == cfg$budget$ext_niter * p$multiplier
  }
  add("carried-over truth 2 species keep their original scores exactly", carried_ok)
  add("refit species match the prescription and are labelled resolution", refit_ok)
  add("resolution budget is the prescribed multiple of the extension", multiplier_ok)
  jobs_now <- read.csv(file.path(study, "scores", "jobs.csv"), stringsAsFactors = FALSE)
  add("jobs.csv marks exactly the prescribed jobs as resolved", setequal(jobs_now$job[jobs_now$resolved], prescribed))
}

out <- do.call(rbind, checks)
dir.create("dev/simstudy/spatial-recoverability-pilot/results", showWarnings = FALSE)
write.csv(out, "dev/simstudy/spatial-recoverability-pilot/results/verification.csv", row.names = FALSE)
print(out, right = FALSE)
if (!all(out$pass)) stop("Verification failed")
cat("Verification passed.\n")
