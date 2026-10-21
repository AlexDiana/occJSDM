# Stage 2 Spatial Representation Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fit five existing spatial representations on 40 of the recoverability pilot's communities and measure, against the pilot's oracle, which representation recovers the most held-out headroom under each truth.

**Architecture:** A standalone study under `dev/simstudy/spatial-stage2-comparison/` that reuses, verbatim where possible, the first comparison's adapter contract, basis, adapters, diagnostics and selection code (copied from branch `codex/spatial-comparison-plan` with recorded hashes), and the pilot's kernel and prediction helpers for the uninformed reference. New code is limited to the data bridge from pilot inputs, a config, a job driver, headroom-share scoring with paired intervals, verification and plotting. Fits run as forked single-threaded workers under the pinned R 4.5.0 runtime; job files go to the research volume; compact results are committed.

**Tech Stack:** R 4.5.0 pinned runtime, occJSDM (the pinned comparison library build), Hmsc and gllvm (provisioned in Task 1), posterior, testthat, jsonlite, `parallel::mclapply`.

**Spec:** `dev/superpowers/specs/2026-10-11-spatial-stage2-comparison-design.md`

## Global Constraints

- Markdown: one line per paragraph, no pipe tables, no em-dashes, escape `-\>` and `\~` in prose and never inside code.
- Interpreter: `dev/simstudy/spatial-stage2-comparison/run-r` only, never system `Rscript`. Single-threaded numerics through the launcher's environment.
- Nothing in `R/`, `src/`, priors, defaults or lessons changes.
- Inputs: the pilot's `inputs/<cell>-rep<k>.rds` for cells `t1-r0.08-n300`, `t1-r0.08-n900`, `t2-n300`, `t2-n900`, replicates 1 to 10, uniform arrangement; MD5s must match the pilot's `input-manifest.csv`. No new community, seed or selection.
- Oracle reference per dataset: pilot `jobs-resolution/<job>.rds` if present, else `jobs/<job>.rds`, with `<job> = <cell>-rep<k>-uniform`; file hash and screen status recorded.
- Arms: `occJSDM/env`, `occJSDM/gp`, `Hmsc/spatial2`, `Hmsc/spatial6`, `gllvm/spatial2`. Nothing tuned against the generating range or lattice truth.
- MCMC: four chains, 3,000 and 4,000; one extension to 8,000 and 16,000; screen Rhat at most 1.01, bulk and tail ESS at least 400 on the first comparison's truth-free training quantities; gllvm three starts plus three.
- Statuses `pass`, `unresolved`, `unstable`, `unsupported`, `no fit`; everything stays in every denominator.
- Attempt ceiling 400. Dispatch only if the all-extension projection at eight workers is under 48 wall-clock hours; otherwise stop and ask Doug. Support points at 900 sites, Hmsc spatial method and prediction thinning are frozen from the timing pilot before dispatch and recorded in the config.
- Meaningful difference 0.10 of oracle headroom; paired t intervals over communities; a comparison with fewer than six paired `pass` communities is inconclusive.
- Commit after each task on branch `claude/spatial-stage2-comparison`; commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.

## Review Focus

- An arm whose prediction is `unsupported` or `no fit` for one dataset: the headroom share must be `NA` for that dataset, the paired comparisons involving that arm must drop that community and report the count, and the cell must stay in the table. Test in Task 6.
- The env arm's RMSE equal to the oracle's RMSE for a species (zero headroom): the share must be `NA`, not infinite, and aggregation must skip it with a recorded count. Test in Task 6.
- A Bayesian attempt that returns but fails the screen twice: status `unresolved`, best prescribed output kept, excluded from the pre-specified interval but present in the descriptive table. Test in Task 7.
- A rerun of `fit` after interruption: completed attempts are reused from disk, never refitted, and selection reads only diagnostics and training predictions. Test in Task 7.
- The gp arm at 900 sites with fewer support points than sites: the recorded realised count must equal the frozen setting, and prediction must use the fitted support set, not the lattice. Test in Task 5.

---

### Task 1: Worktree, launcher, package provisioning and copied sources

**Files:**
- Create: `dev/simstudy/spatial-stage2-comparison/run-r`
- Create: `dev/simstudy/spatial-stage2-comparison/sources.csv`
- Create (copied): `dev/simstudy/spatial-stage2-comparison/upstream/{adapter-contract,basis,adapter-occjsdm,adapter-hmsc,adapter-gllvm,diagnostics,select}.R`
- Create (copied): `dev/simstudy/spatial-stage2-comparison/upstream/tests/{native-core-fixture,native-spatial-fixture,test-core-adapters,test-native-spatial}.R`
- Create: `dev/simstudy/spatial-stage2-comparison/test-runtime.R`
- Modify: `dev/superpowers/specs/2026-10-11-spatial-stage2-comparison-design.md` (one revision-record paragraph)

**Interfaces:**
- Produces: a launcher that loads occJSDM, Hmsc and gllvm; a `sources.csv` with upstream path, upstream commit and MD5 for every copied file.

- [ ] **Step 1: Create the worktree and branch from the pilot branch**

```bash
cd ~/src/occJSDM && git worktree add .worktrees/spatial-stage2-comparison -b claude/spatial-stage2-comparison claude/spatial-recoverability-pilot
cd .worktrees/spatial-stage2-comparison && git log --oneline -1
```
Expected: the pilot branch head (the stage 2 spec commit `591fe84` or later).

- [ ] **Step 2: Copy the upstream sources with recorded hashes**

```bash
mkdir -p dev/simstudy/spatial-stage2-comparison/upstream/tests
U=codex/spatial-comparison-plan; S=dev/simstudy/spatial-package-comparison; T=dev/simstudy/spatial-stage2-comparison/upstream
C=$(git rev-parse $U)
echo "file,upstream_path,upstream_commit,md5" > dev/simstudy/spatial-stage2-comparison/sources.csv
for f in adapter-contract basis adapter-occjsdm adapter-hmsc adapter-gllvm diagnostics select; do
  git show $U:$S/$f.R > $T/$f.R
  echo "upstream/$f.R,$S/$f.R,$C,$(md5 -q $T/$f.R)" >> dev/simstudy/spatial-stage2-comparison/sources.csv
done
for f in native-core-fixture native-spatial-fixture test-core-adapters test-native-spatial; do
  git show $U:$S/tests/$f.R > $T/tests/$f.R
  echo "upstream/tests/$f.R,$S/tests/$f.R,$C,$(md5 -q $T/tests/$f.R)" >> dev/simstudy/spatial-stage2-comparison/sources.csv
done
cat dev/simstudy/spatial-stage2-comparison/sources.csv
```

- [ ] **Step 3: Write the launcher**

The first comparison's launcher points `R_LIBS_USER` at a library under `~/Documents`, which this environment cannot read. The stage 2 launcher points at a study library that Task 1 Step 4 provisions.

```sh
#!/bin/sh
# Pinned R 4.5.0 for the stage 2 spatial representation comparison. Reuses the private
# R home and the occJSDM build from the first comparison; Hmsc and gllvm come from the
# study library STAGE2_LIBRARY (see test-runtime.R). Single-threaded numerics.
set -eu
RUN=/Users/douglasyu/src/occJSDM/dev/simstudy/results
export R_HOME="${STAGE2_R_HOME:-$RUN/spatial-package-comparison-20261008/runtime/R-home}"
export R_SHARE_DIR="$R_HOME/share" R_INCLUDE_DIR="$R_HOME/include" R_DOC_DIR="$R_HOME/doc"
export R_LIBS_USER="${STAGE2_LIBRARY:-$RUN/spatial-stage2-20261011/library}:$RUN/spatial-package-comparison-20261008/library:/Users/douglasyu/Library/R/arm64/4.5/library"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec "$R_HOME/bin/exec/R" --vanilla --slave "$@"
```

- [ ] **Step 4: Provision Hmsc and gllvm into the study library**

Doug chose Route B on 11 October 2026, and it has been done: Hmsc 3.3.7 and gllvm 2.0.15 were installed as CRAN binaries into `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-stage2-20261011/library` under the pinned R 4.5.0 (log in that run directory's `logs/install-packages.log`), with TMB 1.9.22 resolved from the user library. These are the same versions the first comparison used. The executor only verifies (Step 6) and records the route in `results/runtime.json` (Task 8). The two routes below are kept for the record.

Route A, the first comparison's library copied out by Doug. Ask Doug to run, in his own shell, `cp -R "/Users/douglasyu/Documents/Codex/2026-09-10/fam/work/lesson-n-environment-20260922/library" /Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-stage2-20261011/library` and to paste the output of `md5 -q` on `library/Hmsc/DESCRIPTION` and `library/gllvm/DESCRIPTION`. Then verify the copy loads (Step 6) and that the versions are Hmsc 3.3-7 and gllvm 2.0.15, as the first comparison recorded.

Route B, fresh installation under the pinned R if Route A is unavailable:
```bash
mkdir -p /Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-stage2-20261011/library
dev/simstudy/spatial-stage2-comparison/run-r -e 'lib <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-stage2-20261011/library"; install.packages(c("Hmsc","gllvm"), lib = lib, repos = "https://cloud.r-project.org", dependencies = TRUE, Ncpus = 4L); for (p in c("Hmsc","gllvm")) cat(p, as.character(packageVersion(p, lib.loc = lib)), "\n")'
```
Route B gives current CRAN versions, not the first comparison's; the config records the realised versions and the report says so.

- [ ] **Step 5: Write the runtime test**

```r
# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-runtime.R
library(testthat)
test_that("the pinned runtime loads every package the study needs", {
  expect_identical(R.version$major, "4"); expect_identical(R.version$minor, "5.0")
  for (p in c("occJSDM", "Hmsc", "gllvm", "posterior", "jsonlite", "testthat"))
    expect_true(requireNamespace(p, quietly = TRUE), info = p)
  expect_identical(as.character(packageVersion("occJSDM")), "0.1.0")
  expect_true(grepl("spatial-package-comparison-20261008/library", find.package("occJSDM")))
})
test_that("copied upstream sources match their recorded hashes", {
  s <- read.csv("dev/simstudy/spatial-stage2-comparison/sources.csv", stringsAsFactors = FALSE)
  expect_true(all(unname(tools::md5sum(file.path("dev/simstudy/spatial-stage2-comparison", s$file))) == s$md5))
  expect_equal(nrow(s), 11L)
})
cat("Runtime tests passed.\n")
```

- [ ] **Step 6: Run the runtime test**

Run: `chmod +x dev/simstudy/spatial-stage2-comparison/run-r && dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-runtime.R`
Expected: `Runtime tests passed.` If Hmsc or gllvm fails to load, stop and report to Doug; nothing else in this plan can proceed.

- [ ] **Step 7: Record the screen clarification in the spec**

Append one paragraph to the spec's end, before "Future work", under a heading `## Revision record`: "11 October 2026, at planning. The convergence screen for fitted arms uses the first comparison's truth-free training-side quantities as implemented in its `diagnostics.R`: per-species training probability mean, per-site training field, realised field SD, intercepts, coefficients, spatial prior SD and range (occJSDM), coefficients and the species spatial correlation at distance 0.5 (Hmsc), objective, gradient and native convergence (gllvm). The held-out probability average named in the design is not screened, because held-out prediction runs only after selection. No threshold changes." Keep it one line per paragraph.

- [ ] **Step 8: Commit**

```bash
git add dev/simstudy/spatial-stage2-comparison dev/superpowers/specs/2026-10-11-spatial-stage2-comparison-design.md
git commit -m "Scaffold the stage 2 comparison: launcher, provisioned packages and upstream sources"
```

---

### Task 2: Config and the data bridge from pilot inputs

**Files:**
- Create: `dev/simstudy/spatial-stage2-comparison/config.R`
- Create: `dev/simstudy/spatial-stage2-comparison/data.R`
- Create: `dev/simstudy/spatial-stage2-comparison/test-data.R`

**Interfaces:**
- Produces: `stage2_config()` list; `stage2_datasets(pilot_study)` data frame with `dataset_id, cell, truth, n, replicate, input_file, input_md5, oracle_file, oracle_md5`; `load_training(cm, dataset_id, input_md5)` returning a TrainingData list that passes `validate_training`; `load_prediction(cm)` returning the lattice PredictionData; `load_truth(cm)` returning `list(psi_lattice, field_lattice, psi_training, field_training, occupancy, species)`; `load_oracle_reference(path)` returning `list(probability_mean (1600 x 12), field_mean, screen, status, md5)`.

- [ ] **Step 1: Write the failing tests**

```r
# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-data.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("upstream/adapter-contract.R", "upstream/basis.R", "config.R", "data.R")) source(file.path(D, f))
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"

test_that("config fixes arms, cells, budgets and thresholds", {
  cfg <- stage2_config()
  expect_equal(cfg$arms, data.frame(package = c("occJSDM", "occJSDM", "Hmsc", "Hmsc", "gllvm"),
                                    arm = c("env", "gp", "spatial2", "spatial6", "spatial2"), stringsAsFactors = FALSE))
  expect_setequal(cfg$cells, c("t1-r0.08-n300", "t1-r0.08-n900", "t2-n300", "t2-n900"))
  expect_equal(cfg$mcmc$initial[["nburn"]], 3000L); expect_equal(cfg$mcmc$extended[["niter"]], 16000L)
  expect_equal(cfg$delta, 0.10); expect_equal(cfg$limits$attempts, 400L); expect_equal(cfg$wall_hours, 48)
  expect_true(is.na(cfg$frozen$support_points_900)); expect_true(is.na(cfg$frozen$hmsc_method))
})

test_that("the dataset table lists 40 datasets with matching input hashes and oracle files", {
  d <- stage2_datasets(PILOT)
  expect_equal(nrow(d), 40L); expect_equal(as.vector(table(d$cell)), rep(10L, 4))
  expect_true(all(file.exists(d$input_file))); expect_true(all(file.exists(d$oracle_file)))
  man <- read.csv(file.path(PILOT, "input-manifest.csv"), stringsAsFactors = FALSE)
  expect_true(all(d$input_md5 == man$md5[match(basename(d$input_file), man$file)]))
  expect_true(all(grepl("jobs-resolution", d$oracle_file[d$dataset_id %in% c("t1-r0.08-n900-rep2", "t2-n900-rep3")])))
  expect_false(any(grepl("jobs-resolution", d$oracle_file[d$cell == "t1-r0.08-n300" & d$replicate == 1L])))
})

test_that("training, prediction and truth are built from a pilot community and pass the contract", {
  d <- stage2_datasets(PILOT); row <- d[d$dataset_id == "t2-n300-rep1", ]
  cm <- readRDS(row$input_file)
  tr <- load_training(cm, row$dataset_id, row$input_md5)
  expect_true(validate_training(tr))
  expect_equal(dim(tr$y), c(300L, 12L)); expect_equal(tr$species_id, cm$species$species)
  expect_equal(tr$xy, unname(cm$points[cm$index$uniform, ]), ignore_attr = TRUE)
  pr <- load_prediction(cm)
  expect_true(validate_prediction_data(pr, tr)); expect_equal(nrow(pr$xy), 1600L)
  expect_false(any(c("y", "psi", "field") %in% names(pr)))
  truth <- load_truth(cm)
  expect_equal(dim(truth$psi_lattice), c(1600L, 12L)); expect_equal(truth$occupancy, cm$species$occupancy)
  expect_equal(colMeans(truth$psi_lattice), cm$species$occupancy, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("the oracle reference loads with its status and hash", {
  d <- stage2_datasets(PILOT); row <- d[d$dataset_id == "t1-r0.08-n900-rep2", ]
  o <- load_oracle_reference(row$oracle_file)
  expect_equal(dim(o$probability_mean), c(1600L, 12L)); expect_true(all(o$probability_mean > 0 & o$probability_mean < 1))
  expect_true(all(o$screen$status %in% c("pass", "unresolved"))); expect_equal(nchar(o$md5), 32L)
  expect_identical(o$md5, row$oracle_md5)
})
cat("Data tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-data.R`
Expected: error, `config.R` not found.

- [ ] **Step 3: Write config.R**

```r
# Frozen protocol for stage 2. Scientific changes need a dated spec amendment first.
stage2_config <- function() {
  list(schema_version = 1L,
       arms = data.frame(package = c("occJSDM", "occJSDM", "Hmsc", "Hmsc", "gllvm"),
                         arm = c("env", "gp", "spatial2", "spatial6", "spatial2"), stringsAsFactors = FALSE),
       cells = c("t1-r0.08-n300", "t1-r0.08-n900", "t2-n300", "t2-n900"), replicates = 1:10, arrangement = "uniform",
       gate_occupancy = c(0.25, 0.50, 0.75),
       mcmc = list(initial = c(nchain = 4L, nburn = 3000L, niter = 4000L, nthin = 1L),
                   extended = c(nchain = 4L, nburn = 8000L, niter = 16000L, nthin = 1L)),
       extension_factor = 4.43,
       diagnostics = list(rhat = 1.01, ess = 400, gradient = 0.01, objective_agreement = 0.1,
                          probability_agreement = 0.01, initial_starts = 3L, max_starts = 6L),
       gllvm = list(maxit = 12000L, method = "VA", family = "binomial", link = "logit"),
       support_points_300 = 300L,
       # Frozen from the timing pilot before dispatch (run.R --mode=freeze); NA until then.
       frozen = list(support_points_900 = NA_integer_, hmsc_method = NA_character_, hmsc_knots = NA_integer_, prediction_thin = NA_integer_),
       ladder = list(support_points_900 = c(900L, 600L, 450L, 300L), hmsc = list(c("Full", NA), c("GPP", 300L))),
       delta = 0.10, min_paired = 6L,
       comparisons = list(c("occJSDM/gp", "Hmsc/spatial2"), c("occJSDM/gp", "Hmsc/spatial6"), c("occJSDM/gp", "gllvm/spatial2"),
                          c("Hmsc/spatial6", "Hmsc/spatial2"), c("Hmsc/spatial2", "gllvm/spatial2")),
       limits = list(attempts = 400L, jobs = 200L), wall_hours = 48, workers = 8L, prediction_block = 25L)
}

# study_config() is the name the upstream diagnostics and selection code call.
study_config <- function() {
  cfg <- stage2_config()
  list(diagnostics = cfg$diagnostics, gllvm = cfg$gllvm, mcmc = cfg$mcmc, prediction_block = cfg$prediction_block)
}

arm_key <- function(package, arm) paste(package, arm, sep = "/")
```

- [ ] **Step 4: Write data.R**

```r
# Bridge from the pilot's community objects to the upstream TrainingData and PredictionData
# contracts, plus the oracle reference. Truth never enters TrainingData or PredictionData.
stage2_datasets <- function(pilot_study, config = stage2_config()) {
  man <- read.csv(file.path(pilot_study, "input-manifest.csv"), stringsAsFactors = FALSE)
  rows <- list()
  for (cell in config$cells) for (r in config$replicates) {
    file <- sprintf("%s-rep%d.rds", cell, r)
    m <- man[man$file == file, ]; stopifnot(nrow(m) == 1L)
    job <- sprintf("%s-rep%d-%s", cell, r, config$arrangement)
    resolved <- file.path(pilot_study, "jobs-resolution", paste0(job, ".rds"))
    oracle <- if (file.exists(resolved)) resolved else file.path(pilot_study, "jobs", paste0(job, ".rds"))
    stopifnot(file.exists(oracle))
    rows[[length(rows) + 1L]] <- data.frame(dataset_id = sprintf("%s-rep%d", cell, r), cell = cell, truth = substr(cell, 1, 2),
      n = m$n, replicate = r, input_file = file.path(pilot_study, "inputs", file), input_md5 = m$md5,
      oracle_file = oracle, oracle_md5 = unname(tools::md5sum(oracle)), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

site_ids <- function(n) sprintf("site%04d", seq_len(n))
lattice_ids <- function(n) sprintf("lattice%04d", seq_len(n))

load_training <- function(cm, dataset_id, input_md5) {
  idx <- cm$index$uniform
  xy <- unname(cm$points[idx, , drop = FALSE]); dimnames(xy) <- NULL
  y <- cm$z[idx, , drop = FALSE]; colnames(y) <- cm$species$species
  list(schema_version = 1L, dataset_id = dataset_id, site_id = site_ids(length(idx)), species_id = cm$species$species,
       y = y, environment = matrix(cm$environment[idx], ncol = 1L), xy = xy, input_hashes = c(input = input_md5))
}

load_prediction <- function(cm) {
  idx <- cm$index$lattice
  xy <- unname(cm$points[idx, , drop = FALSE]); dimnames(xy) <- NULL
  list(schema_version = 1L, site_id = lattice_ids(length(idx)), species_id = cm$species$species,
       environment = matrix(cm$environment[idx], ncol = 1L), xy = xy)
}

load_truth <- function(cm) {
  lat <- cm$index$lattice; idx <- cm$index$uniform
  list(psi_lattice = cm$psi[lat, , drop = FALSE], field_lattice = cm$spatial[lat, , drop = FALSE],
       psi_training = cm$psi[idx, , drop = FALSE], field_training = cm$spatial[idx, , drop = FALSE],
       eta_lattice = sweep(cm$environment[lat] %o% cm$B, 2, cm$B0, "+"),
       occupancy = cm$species$occupancy, species = cm$species$species, spatial_sd = cm$spatial_sd)
}

load_oracle_reference <- function(path) {
  r <- readRDS(path)
  list(probability_mean = r$prediction$probability_mean, field_mean = r$prediction$field_mean,
       screen = r$screen, status = if (all(r$screen$status == "pass")) "pass" else "unresolved",
       phases = r$phases, md5 = unname(tools::md5sum(path)), job = r$job)
}
```

- [ ] **Step 5: Run to verify it passes**

Run: same command. Expected: `Data tests passed.`

- [ ] **Step 6: Commit**

```bash
git add dev/simstudy/spatial-stage2-comparison/config.R dev/simstudy/spatial-stage2-comparison/data.R dev/simstudy/spatial-stage2-comparison/test-data.R
git commit -m "Add the stage 2 config and the data bridge from pilot communities"
```

---

### Task 3: Adapter adjustments for the new arms and frozen settings

**Files:**
- Modify: `dev/simstudy/spatial-stage2-comparison/upstream/basis.R` (arm list)
- Modify: `dev/simstudy/spatial-stage2-comparison/upstream/adapter-hmsc.R` (factor count and spatial method from the arm and settings)
- Modify: `dev/simstudy/spatial-stage2-comparison/upstream/adapter-occjsdm.R` (support points from settings)
- Modify: `dev/simstudy/spatial-stage2-comparison/upstream/diagnostics.R` (Hmsc arm test)
- Create: `dev/simstudy/spatial-stage2-comparison/test-adapters.R`
- Modify: `dev/simstudy/spatial-stage2-comparison/sources.csv` (mark modified files)

**Interfaces:**
- Consumes: `fit_*`/`predict_*` from upstream; `settings` now carries `support_points`, `hmsc_method`, `hmsc_knots`.
- Produces: `spatial_factors(arm)` returning 2L or 6L for `spatial2`/`spatial6`; `fit_basis` accepting `spatial6`.

- [ ] **Step 1: Make the four edits**

In `upstream/basis.R`, change the arm check and the two `%in% c("gp","spatial2")` tests:
```r
  if (!arm %in% c("env","trend","mem","dnn","gp","spatial2","spatial6")) stop("Unknown spatial arm")
  ...
  if (arm %in% c("gp","spatial2","spatial6")) b$training_matrix <- xy
  ...
  if (basis$arm %in% c("gp","spatial2","spatial6")) return(xy)
```

In `upstream/adapter-hmsc.R`, replace the level construction and the Lambda check:
```r
spatial_factors <- function(arm) { k <- suppressWarnings(as.integer(sub("^spatial", "", arm))); if (is.na(k)) stop("Not a spatial factor arm: ", arm); k }
...
  if(grepl("^spatial[0-9]+$", basis$arm)) {
    xy <- training$xy; rownames(xy) <- training$site_id; nf <- spatial_factors(basis$arm)
    method <- if (is.null(settings$hmsc_method)) "Full" else settings$hmsc_method
    # Hmsc::constructKnots takes a grid density, not a count; k-means with the attempt seed gives an exact, recorded knot set.
    level <- if (method == "Full") Hmsc::HmscRandomLevel(sData=xy,sMethod="Full") else {
      set.seed(settings$seed); knots <- stats::kmeans(xy,centers=settings$hmsc_knots,nstart=5L,iter.max=100L)$centers
      Hmsc::HmscRandomLevel(sData=xy,sMethod="GPP",sKnot=knots) }
    level <- Hmsc::setPriors(level,nfMin=nf,nfMax=nf)
    args$studyDesign <- data.frame(site=factor(training$site_id,levels=training$site_id))
    args$ranLevels <- list(site=level)
  }
  ...
  if(grepl("^spatial[0-9]+$", basis$arm)) stopifnot(all(vapply(native$postList,function(ch)
    all(vapply(ch,function(d) nrow(d$Lambda[[1]])==spatial_factors(basis$arm),logical(1))),logical(1))))
```
and in `hmsc_draw_extractor` replace both `basis$arm=="spatial2"` tests with `grepl("^spatial[0-9]+$", basis$arm)`. The GPP route's conditional prediction must use Hmsc's own `predictLatentFactor`, so when `settings$hmsc_method == "GPP"` the extractor calls `Hmsc::predictLatentFactor(unitsPred, units, postEta, postAlpha, rL)` per chain instead of `hmsc_conditional_draws`, and the conditioning string records "GPP".

In `upstream/adapter-gllvm.R`, replace `control.start=list(starting.val="res",n.init=1L,jitter.var=.1)` with `control.start=list(starting.val="random",n.init=1L,jitter.var=.1)`. The residual start fails with "Calculating starting values failed" on any dataset with a training site where no species is present, and 29 of the 80 pilot inputs at 300 and 900 sites have one; random starts seeded by the attempt seed succeed and make the three-start replication check meaningful. Record `starting_val="random"` in the fit's `extra` list.

In `upstream/adapter-occjsdm.R`, replace `n_supportpoints=100L` with `n_supportpoints=as.integer(settings$support_points)` and add after the fit `stopifnot(nrow(native$infos$list_Xs$X_tilde)==as.integer(settings$support_points))` for the gp arm only; record `support_points=settings$support_points` in the `extra` list.

In `upstream/diagnostics.R`, replace `if(f$arm=="spatial2")` with `if(grepl("^spatial[0-9]+$", f$arm))`.

Update `sources.csv`: add a column `modified` with `yes` for the four edited files and `no` otherwise, and recompute MD5s for all rows.

- [ ] **Step 2: Write the adapter tests**

```r
# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-adapters.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("upstream/adapter-contract.R", "upstream/basis.R", "upstream/adapter-occjsdm.R", "upstream/adapter-hmsc.R",
            "upstream/adapter-gllvm.R", "upstream/diagnostics.R", "upstream/select.R", "config.R", "data.R")) source(file.path(D, f))
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
d <- stage2_datasets(PILOT); row <- d[d$dataset_id == "t2-n300-rep1", ]
cm <- readRDS(row$input_file); tr <- load_training(cm, row$dataset_id, row$input_md5); pr <- load_prediction(cm)
# A 40-site, 4-species slice keeps every fit under a minute.
slice <- function(tr, sites = 1:40, species = c(4L, 6L, 9L, 12L)) {
  tr$site_id <- tr$site_id[sites]; tr$species_id <- tr$species_id[species]; tr$y <- tr$y[sites, species, drop = FALSE]
  tr$environment <- tr$environment[sites, , drop = FALSE]; tr$xy <- tr$xy[sites, , drop = FALSE]; tr }
small <- slice(tr); pr_small <- pr; pr_small$species_id <- small$species_id; pr_small <- subset_prediction(pr_small, 1:30)
tiny <- c(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L)
settings <- function(...) list(seed = 11L, job_id = 1L, attempt_id = 1L, source_hashes = c(x = "0"), mcmc = tiny,
                               gllvm = stage2_config()$gllvm, prediction_seed = 12L, ...)

test_that("spatial6 is a valid basis and Hmsc fits six factors", {
  b <- fit_basis(small, "spatial6"); expect_identical(b$arm, "spatial6"); expect_equal(dim(predict_basis(b, pr_small)), c(30L, 2L))
  f <- fit_hmsc(small, b, settings(hmsc_method = "Full"))
  expect_equal(nrow(f$native$postList[[1]][[1]]$Lambda[[1]]), 6L)
  p <- predict_hmsc(f, pr_small, b); expect_true(validate_prediction(p, pr_small))
})

test_that("Hmsc GPP route fits and predicts with recorded knots", {
  b <- fit_basis(small, "spatial2")
  f <- fit_hmsc(small, b, settings(hmsc_method = "GPP", hmsc_knots = 10L))
  expect_equal(nrow(f$native$ranLevels$site$sKnot), 10L)
  p <- predict_hmsc(f, pr_small, b); expect_true(validate_prediction(p, pr_small)); expect_true(grepl("GPP", p$conditioning))
})

test_that("occJSDM gp honours the support-point setting and records it", {
  b <- fit_basis(small, "gp")
  f <- fit_occjsdm(small, b, settings(support_points = 20L))
  expect_equal(nrow(f$native$infos$list_Xs$X_tilde), 20L); expect_equal(f$support_points, 20L)
  p <- predict_occjsdm(f, pr_small, b); expect_true(validate_prediction(p, pr_small))
  expect_error(fit_occjsdm(small, b, settings(support_points = 41L)))
})

test_that("gllvm fits a dataset with an all-zero training site", {
  zero <- small; zero$y[1, ] <- 0L   # a site where no species is present
  bg <- fit_basis(zero, "spatial2"); fz <- fit_gllvm(zero, bg, settings())
  expect_identical(fz$starting_val, "random"); expect_true(is.finite(fz$native$objective))
})

test_that("env and gllvm arms still fit, and diagnostics run for each package", {
  fe <- fit_occjsdm(small, fit_basis(small, "env"), settings(support_points = NA)); de <- diagnose_fit(fe, small, fit_basis(small, "env"))
  expect_true(all(c("failed", "reason") %in% names(de)))
  bg <- fit_basis(small, "spatial2"); fg <- fit_gllvm(small, bg, settings()); dg <- diagnose_fit(fg, small, bg)
  expect_equal(nrow(dg), 1L); pg <- predict_gllvm(fg, pr_small, bg); expect_true(validate_prediction(pg, pr_small))
})
cat("Adapter tests passed.\n")
```

- [ ] **Step 3: Run the adapter tests and the copied upstream certification tests**

Run: `dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-adapters.R`
Expected: `Adapter tests passed.` Then run the copied upstream tests, adjusting their `source()` paths to `upstream/` and their config to `study_config()` from `config.R`:
`dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/upstream/tests/test-native-spatial.R` and `test-core-adapters.R`. Expected: both pass. If a copied test depends on sjSDM or the Python runtime, delete that test block and record the deletion in `sources.csv`'s `modified` column.

- [ ] **Step 4: Commit**

```bash
git add dev/simstudy/spatial-stage2-comparison
git commit -m "Adapt the upstream adapters for six-factor Hmsc, GPP, and configurable support points"
```

---

### Task 4: Scoring: headroom share, pattern, slope, aggregation

**Files:**
- Create: `dev/simstudy/spatial-stage2-comparison/score.R`
- Create: `dev/simstudy/spatial-stage2-comparison/test-score.R`

**Interfaces:**
- Consumes: prediction records (`p_full`, `delta_p`, `field_mean`, `link`), `load_truth`, `load_oracle_reference`, the pilot's `kernel.R` (`gauss_hermite`) and `predict.R` (`uninformed_reference`), sourced from `dev/simstudy/spatial-recoverability-pilot/`.
- Produces: `score_arm(prediction, env_prediction, oracle, truth, dataset_id, package, arm, status)` returning one row per species with columns `dataset_id, package, arm, key, species, occupancy, status, rmse_env, rmse_oracle, oracle_status, rmse, mae, headroom_share, pattern, field_slope, field_correlation` (training-site scores are taken in Task 8 from the stored training predictions and reported, never used in decisions); `aggregate_arm(scores, gate_occupancy)` giving one row per dataset and arm; `paired_comparisons(agg, config)` and `development_reading(pairs, config)`.

- [ ] **Step 1: Write the failing tests**

```r
# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-score.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("config.R", "score.R")) source(file.path(D, f))
set.seed(1)
truth <- list(psi_lattice = matrix(plogis(rnorm(1600 * 3)), 1600, 3), field_lattice = matrix(rnorm(1600 * 3), 1600, 3),
              eta_lattice = matrix(0, 1600, 3), occupancy = c(0.25, 0.5, 0.75), species = paste0("s", 1:3), spatial_sd = rep(1, 3),
              psi_training = matrix(0.5, 10, 3), field_training = matrix(0, 10, 3))
pred <- function(p, field = NULL, link = "logit") list(p_full = p, delta_p = p - 0.5, field_mean = if (is.null(field)) p * 0 else field, link = link,
                                                        site_id = paste0("l", 1:1600), species_id = truth$species)
env <- pred(matrix(0.5, 1600, 3)); oracle <- list(probability_mean = truth$psi_lattice, status = "pass")

test_that("headroom share is 1 at the oracle, 0 at the env arm, negative when worse, NA at zero headroom", {
  s1 <- score_arm(pred(truth$psi_lattice), env, oracle, truth, "d", "occJSDM", "gp", "pass")
  expect_equal(s1$headroom_share, rep(1, 3), tolerance = 1e-12)
  s0 <- score_arm(env, env, oracle, truth, "d", "occJSDM", "env", "pass"); expect_equal(s0$headroom_share, rep(0, 3))
  worse <- pred(matrix(0.9, 1600, 3)); sw <- score_arm(worse, env, oracle, truth, "d", "Hmsc", "spatial2", "pass")
  expect_true(all(sw$headroom_share < 0))
  flat_oracle <- list(probability_mean = matrix(0.5, 1600, 3), status = "pass")
  sz <- score_arm(pred(truth$psi_lattice), env, flat_oracle, truth, "d", "occJSDM", "gp", "pass"); expect_true(all(is.na(sz$headroom_share)))
})

test_that("slope and correlation are computed only for logit arms and pattern for every arm", {
  s <- score_arm(pred(truth$psi_lattice, 0.5 * truth$field_lattice, "logit"), env, oracle, truth, "d", "occJSDM", "gp", "pass")
  expect_equal(s$field_slope, rep(0.5, 3), tolerance = 1e-10); expect_equal(s$field_correlation, rep(1, 3), tolerance = 1e-10)
  h <- score_arm(pred(truth$psi_lattice, 0.5 * truth$field_lattice, "probit"), env, oracle, truth, "d", "Hmsc", "spatial2", "pass")
  expect_true(all(is.na(h$field_slope))); expect_false(any(is.na(h$pattern)))
})

test_that("unsupported fits give NA scores with their status and stay in the table", {
  s <- score_arm(NULL, env, oracle, truth, "d", "gllvm", "spatial2", "unsupported")
  expect_equal(nrow(s), 3L); expect_true(all(is.na(s$headroom_share))); expect_true(all(s$status == "unsupported"))
})

test_that("aggregation pools gate occupancies and records how many species were skipped", {
  s <- rbind(score_arm(pred(truth$psi_lattice), env, oracle, truth, "d1", "occJSDM", "gp", "pass"),
             score_arm(pred(truth$psi_lattice), env, oracle, truth, "d2", "occJSDM", "gp", "pass"))
  s$headroom_share[1] <- NA
  a <- aggregate_arm(s, c(0.25, 0.5, 0.75))
  expect_equal(nrow(a), 2L); expect_equal(a$skipped_species, c(1L, 0L)); expect_equal(a$headroom_share, c(1, 1))
})

test_that("paired comparisons apply the delta rule and the minimum pairing", {
  cfg <- stage2_config()
  mk <- function(arm, h, status = "pass") data.frame(cell = "t1-r0.08-n900", dataset_id = sprintf("d%d", seq_along(h)), package = sub("/.*", "", arm),
                                                      arm = sub(".*/", "", arm), key = arm, headroom_share = h, status = status, stringsAsFactors = FALSE)
  agg <- rbind(mk("occJSDM/gp", rep(0.5, 10)), mk("Hmsc/spatial2", rep(0.5, 10) + 0.15 + rnorm(10, 0, 0.02)),
               mk("Hmsc/spatial6", rep(0.5, 10) + 0.03), mk("gllvm/spatial2", rep(0.5, 10) + 0.15, status = c(rep("pass", 5), rep("unstable", 5))))
  p <- paired_comparisons(agg, cfg)
  r <- p[p$a == "occJSDM/gp" & p$b == "Hmsc/spatial2", ]; expect_identical(r$verdict, "b clearly better"); expect_equal(r$paired, 10L)
  r6 <- p[p$a == "occJSDM/gp" & p$b == "Hmsc/spatial6", ]; expect_identical(r6$verdict, "no clear difference")
  rg <- p[p$a == "occJSDM/gp" & p$b == "gllvm/spatial2", ]; expect_identical(rg$verdict, "inconclusive"); expect_equal(rg$paired, 5L)
})

test_that("the development reading covers the four branches", {
  v <- function(t1_gp_vs_factor, t2_gp_vs_factor) data.frame(
    cell = c("t1-r0.08-n900", "t2-n900"), a = "occJSDM/gp", b = "Hmsc/spatial2", verdict = c(t1_gp_vs_factor, t2_gp_vs_factor), stringsAsFactors = FALSE)
  expect_identical(development_reading(v("b clearly better", "no clear difference"), 0.6)$reading, "shared factors")
  expect_identical(development_reading(v("no clear difference", "a clearly better"), 0.6)$reading, "per-species ranges")
  expect_identical(development_reading(v("b clearly better", "a clearly better"), 0.6)$reading, "trade-off")
  expect_identical(development_reading(v("no clear difference", "no clear difference"), 0.3)$reading, "representation not the bottleneck")
  expect_identical(development_reading(v("no clear difference", "no clear difference"), 0.7)$reading, "no clear difference")
})
cat("Score tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-score.R`
Expected: error, `score.R` not found.

- [ ] **Step 3: Write score.R**

```r
# Scores for stage 2: headroom share against the pilot oracle, pattern, slope, aggregation,
# paired comparisons and the development reading. Requires config.R.
PILOT_CODE <- "dev/simstudy/spatial-recoverability-pilot"
source(file.path(PILOT_CODE, "kernel.R")); source(file.path(PILOT_CODE, "predict.R"))

rmse <- function(a, b) sqrt(mean((a - b)^2))

score_arm <- function(prediction, env_prediction, oracle, truth, dataset_id, package, arm, status) {
  S <- length(truth$species)
  rows <- lapply(seq_len(S), function(s) {
    psi <- truth$psi_lattice[, s]; unif <- uninformed_reference(truth$eta_lattice[, s], truth$spatial_sd[s])
    r_env <- rmse(env_prediction$p_full[, s], psi); r_or <- rmse(oracle$probability_mean[, s], psi)
    base <- data.frame(dataset_id = dataset_id, package = package, arm = arm, key = arm_key(package, arm), species = truth$species[s],
                       occupancy = truth$occupancy[s], status = status, rmse_env = r_env, rmse_oracle = r_or, oracle_status = oracle$status,
                       stringsAsFactors = FALSE)
    if (is.null(prediction)) return(cbind(base, rmse = NA_real_, mae = NA_real_, headroom_share = NA_real_, pattern = NA_real_,
                                          field_slope = NA_real_, field_correlation = NA_real_))
    p <- prediction$p_full[, s]; r <- rmse(p, psi)
    share <- if (r_env - r_or > 1e-12) (r_env - r) / (r_env - r_or) else NA_real_
    true_contrast <- psi - unif; pattern <- if (sd(prediction$delta_p[, s]) > 0) cor(prediction$delta_p[, s], true_contrast) else NA_real_
    slope <- corr <- NA_real_
    if (identical(prediction$link, "logit")) {
      fc <- prediction$field_mean[, s] - mean(prediction$field_mean[, s]); tc <- truth$field_lattice[, s] - mean(truth$field_lattice[, s])
      slope <- sum(fc * tc) / sum(tc^2); corr <- if (sd(fc) > 0) cor(fc, tc) else NA_real_
    }
    cbind(base, rmse = r, mae = mean(abs(p - psi)), headroom_share = share, pattern = pattern, field_slope = slope, field_correlation = corr)
  })
  do.call(rbind, rows)
}

aggregate_arm <- function(scores, gate_occupancy) {
  s <- scores[scores$occupancy %in% gate_occupancy, ]
  key <- c("dataset_id", "package", "arm", "key")
  metrics <- c("headroom_share", "pattern", "field_slope", "field_correlation", "rmse", "rmse_env", "rmse_oracle")
  by_level <- aggregate(s[metrics], s[c(key, "occupancy")], function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_)
  out <- aggregate(by_level[metrics], by_level[key], function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_)
  skipped <- aggregate(list(skipped_species = !is.finite(s$headroom_share)), s[key], sum)
  status <- aggregate(list(status = s$status), s[key], function(x) x[1])
  out <- merge(merge(out, skipped, by = key), status, by = key)
  out$cell <- sub("-rep[0-9]+$", "", out$dataset_id); out
}

paired_comparisons <- function(agg, config) {
  rows <- list()
  for (cell in unique(agg$cell)) for (cmp in config$comparisons) {
    a <- agg[agg$cell == cell & agg$key == cmp[1], ]; b <- agg[agg$cell == cell & agg$key == cmp[2], ]
    m <- merge(a, b, by = "dataset_id", suffixes = c("_a", "_b"))
    ok <- m$status_a == "pass" & m$status_b == "pass" & is.finite(m$headroom_share_a) & is.finite(m$headroom_share_b)
    d <- m$headroom_share_b[ok] - m$headroom_share_a[ok]; k <- length(d)
    est <- if (k) mean(d) else NA_real_; se <- if (k > 1) sd(d) / sqrt(k) else NA_real_
    lo <- est - qt(0.975, max(k - 1, 1)) * se; hi <- est + qt(0.975, max(k - 1, 1)) * se
    verdict <- if (k < config$min_paired) "inconclusive" else if (est >= config$delta && lo > 0) "b clearly better" else
      if (est <= -config$delta && hi < 0) "a clearly better" else "no clear difference"
    rows[[length(rows) + 1L]] <- data.frame(cell = cell, a = cmp[1], b = cmp[2], paired = k, excluded = nrow(m) - k,
                                            estimate = est, lower = lo, upper = hi, verdict = verdict, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# `pairs` holds gp against each factor arm by cell; `max_share` is the largest headroom share any arm reached in any cell.
development_reading <- function(pairs, max_share) {
  gp <- pairs[pairs$a == "occJSDM/gp" & pairs$b %in% c("Hmsc/spatial2", "Hmsc/spatial6", "gllvm/spatial2"), ]
  t1 <- gp[grepl("^t1", gp$cell) & grepl("n900", gp$cell), ]; t2 <- gp[grepl("^t2", gp$cell) & grepl("n900", gp$cell), ]
  factors_better_t1 <- any(t1$verdict == "b clearly better"); gp_worse_t2 <- any(t2$verdict == "b clearly better")
  gp_better_t2 <- any(t2$verdict == "a clearly better"); gp_worse_t1 <- any(t1$verdict == "b clearly better")
  reading <- if (factors_better_t1 && !gp_better_t2 && !gp_worse_t2) "shared factors" else
    if (!gp_worse_t1 && gp_better_t2) "per-species ranges" else
    if (factors_better_t1 && gp_better_t2) "trade-off" else
    if (all(c(t1$verdict, t2$verdict) %in% c("no clear difference", "inconclusive")) && is.finite(max_share) && max_share < 0.5) "representation not the bottleneck" else
    "no clear difference"
  list(reading = reading, truth1_900 = t1, truth2_900 = t2, max_share = max_share)
}
```

- [ ] **Step 4: Run to verify it passes**

Run: same command. Expected: `Score tests passed.`

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/spatial-stage2-comparison/score.R dev/simstudy/spatial-stage2-comparison/test-score.R
git commit -m "Add stage 2 scoring: headroom share, paired comparisons and the development reading"
```

---

### Task 5: Job driver: attempts, selection, prediction, storage

**Files:**
- Create: `dev/simstudy/spatial-stage2-comparison/jobs.R`
- Create: `dev/simstudy/spatial-stage2-comparison/test-jobs.R`

**Interfaces:**
- Consumes: Tasks 2 to 4, upstream `fit_*`, `diagnose_fit`, `select_attempts`, `predict_*`.
- Produces: `job_table(study)` with one row per dataset and arm (200) and `job_id`; `attempt_settings(job, attempt_id, config, frozen)`; `run_attempt(study, job, attempt_id)` writing `attempts/<job_id>-<attempt_id>.rds` and returning a summary with `status`, `acceptable`, `objective`, `training_probability`; `select_job(study, job)` applying `select_attempts`; `predict_job(study, job)` writing `predictions/<job_id>.rds` with the upstream prediction record, the thinning used and timing; `run_job(study, job)` driving the attempt loop to a terminal status.

- [ ] **Step 1: Write the failing test**

```r
# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-jobs.R
library(testthat)
source("dev/simstudy/spatial-stage2-comparison/jobs.R")
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
study <- file.path(tempdir(), paste0("stage2-smoke-", as.integer(Sys.time())))
tiny <- list(initial = c(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L), extended = c(nchain = 2L, nburn = 40L, niter = 80L, nthin = 1L))

test_that("prepare writes config and the dataset table before any attempt", {
  stage2_prepare(study, PILOT, jobs_dir = file.path(study, "jobs-local"))
  expect_true(file.exists(file.path(study, "config.rds"))); d <- read.csv(file.path(study, "datasets.csv")); expect_equal(nrow(d), 40L)
  j <- job_table(study); expect_equal(nrow(j), 200L); expect_false(anyDuplicated(j$job_id) > 0)
  expect_error(stage2_prepare(study, PILOT, jobs_dir = file.path(study, "jobs-local")))
})

test_that("a Bayesian job runs initial, extends on screen failure, selects, predicts and stores", {
  j <- job_table(study); job <- j[j$dataset_id == "t2-n300-rep1" & j$key == "occJSDM/gp", ]
  frozen <- list(support_points_900 = 300L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L)
  r <- run_job(study, job, mcmc = tiny, frozen = frozen, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))
  expect_true(r$status %in% c("pass", "unresolved"))   # tiny budgets fail the screen; extension then follows
  expect_equal(length(list.files(file.path(study, "jobs-local", "attempts"), pattern = paste0("^", job$job_id, "-"))), 2L)
  p <- readRDS(file.path(study, "jobs-local", "predictions", paste0(job$job_id, ".rds")))
  expect_equal(dim(p$prediction$p_full), c(1600L, 4L)); expect_equal(p$prediction_thin, 10L)
  expect_identical(run_job(study, job, mcmc = tiny, frozen = frozen, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))$status, r$status)
  expect_equal(length(list.files(file.path(study, "jobs-local", "attempts"), pattern = paste0("^", job$job_id, "-"))), 2L)   # reused, not refitted
})

test_that("a gllvm job runs three starts and selects or marks unstable", {
  j <- job_table(study); job <- j[j$dataset_id == "t2-n300-rep1" & j$key == "gllvm/spatial2", ]
  frozen <- list(support_points_900 = 300L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 1L)
  r <- run_job(study, job, mcmc = tiny, frozen = frozen, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))
  expect_true(r$status %in% c("pass", "unstable", "no fit"))
  expect_true(length(list.files(file.path(study, "jobs-local", "attempts"), pattern = paste0("^", job$job_id, "-"))) >= 3L)
})
cat("Job tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-jobs.R`
Expected: error, `jobs.R` not found.

- [ ] **Step 3: Write jobs.R**

```r
# Job driver for stage 2: attempts, selection, prediction and storage. One job is one dataset and arm.
D2 <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("upstream/adapter-contract.R", "upstream/basis.R", "upstream/adapter-occjsdm.R", "upstream/adapter-hmsc.R",
            "upstream/adapter-gllvm.R", "upstream/diagnostics.R", "upstream/select.R", "config.R", "data.R", "score.R")) source(file.path(D2, f))

stage2_prepare <- function(study, pilot_study, jobs_dir, config = stage2_config()) {
  dir.create(study, showWarnings = FALSE, recursive = TRUE)
  stopifnot(!file.exists(file.path(study, "config.rds")))
  dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)
  for (d in c("attempts", "predictions")) dir.create(file.path(jobs_dir, d), showWarnings = FALSE)
  if (!file.exists(file.path(study, "jobs"))) file.symlink(normalizePath(jobs_dir), file.path(study, "jobs"))
  for (d in c("status", "scores", "logs")) dir.create(file.path(study, d), showWarnings = FALSE)
  datasets <- stage2_datasets(pilot_study, config)
  config$pilot_study <- pilot_study; config$created <- Sys.time(); config$git <- system("git rev-parse HEAD", intern = TRUE)
  config$r_version <- R.version.string
  config$packages <- sapply(c("occJSDM", "Hmsc", "gllvm", "posterior"), function(p) as.character(packageVersion(p)))
  config$sources <- tools::md5sum(list.files(D2, pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
  saveRDS(config, file.path(study, "config.rds"))
  write.csv(datasets, file.path(study, "datasets.csv"), row.names = FALSE)
  invisible(config)
}

job_table <- function(study) {
  cfg <- readRDS(file.path(study, "config.rds")); d <- read.csv(file.path(study, "datasets.csv"), stringsAsFactors = FALSE)
  jobs <- merge(d, cfg$arms, by = NULL)
  jobs$key <- arm_key(jobs$package, jobs$arm)
  jobs$job_id <- sprintf("%s--%s-%s", jobs$dataset_id, jobs$package, jobs$arm)
  jobs$kind <- ifelse(jobs$package == "gllvm", "point", "bayesian")
  jobs[order(jobs$job_id), ]
}

stage2_seed <- function(job_id, attempt_id, kind) {
  b <- as.numeric(utf8ToInt(paste0("stage2:", kind, ":", job_id)))
  as.integer((sum(b * seq_along(b)) * 7919 + sum(b) + 97 * attempt_id) %% 2147483646 + 1)
}

attempt_settings <- function(job, attempt_id, config, frozen, mcmc = config$mcmc) {
  sched <- if (job$kind == "bayesian" && attempt_id == 2L) mcmc$extended else mcmc$initial
  list(seed = stage2_seed(job$job_id, attempt_id, "fit"), prediction_seed = stage2_seed(job$job_id, 1L, "prediction"),
       job_id = job$job_id, attempt_id = attempt_id, source_hashes = config$sources, mcmc = sched, gllvm = config$gllvm,
       support_points = if (job$n == 900L) frozen$support_points_900 else config$support_points_300,
       hmsc_method = frozen$hmsc_method, hmsc_knots = frozen$hmsc_knots)
}

slice_training <- function(tr, sites, species) {
  if (is.null(sites)) return(tr)
  tr$site_id <- tr$site_id[sites]; tr$species_id <- tr$species_id[species]; tr$y <- tr$y[sites, species, drop = FALSE]
  tr$environment <- tr$environment[sites, , drop = FALSE]; tr$xy <- tr$xy[sites, , drop = FALSE]; tr
}

run_attempt <- function(study, job, attempt_id, config, frozen, mcmc, slice_sites = NULL, slice_species = NULL) {
  path <- file.path(study, "jobs", "attempts", sprintf("%s-%d.rds", job$job_id, attempt_id))
  if (file.exists(path)) return(readRDS(path)$summary)
  cm <- readRDS(job$input_file); tr <- slice_training(load_training(cm, job$dataset_id, job$input_md5), slice_sites, slice_species)
  basis <- fit_basis(tr, job$arm); settings <- attempt_settings(job, attempt_id, config, frozen, mcmc)
  t0 <- Sys.time()
  fit <- tryCatch(switch(job$package, occJSDM = fit_occjsdm(tr, basis, settings), Hmsc = fit_hmsc(tr, basis, settings), gllvm = fit_gllvm(tr, basis, settings)),
                  error = function(e) structure(list(message = conditionMessage(e)), class = "fit_error"))
  if (inherits(fit, "fit_error")) {
    summary <- list(attempt_id = attempt_id, status = "failed", acceptable = FALSE, error = fit$message, elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"))
    saveRDS(list(summary = summary, fit = NULL, diagnostics = NULL), path); return(summary)
  }
  diag <- diagnose_fit(fit, tr, basis, study_config())
  training_probability <- if (job$kind == "point") { p <- predict_gllvm(fit, c(tr[c("schema_version", "site_id", "species_id", "environment", "xy")]), basis); p$p_full } else NULL
  summary <- list(attempt_id = attempt_id, status = "complete", acceptable = !any(diag$failed), objective = if (job$kind == "point") fit$native$objective else NA_real_,
                  training_probability = training_probability, elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"), warnings = fit$warnings)
  saveRDS(list(summary = summary, fit = fit, diagnostics = diag, basis = basis, training = tr), path)
  summary
}

select_job <- function(study, job, config) {
  files <- list.files(file.path(study, "jobs", "attempts"), pattern = paste0("^", job$job_id, "-[0-9]+\\.rds$"), full.names = TRUE)
  summaries <- lapply(files, function(f) readRDS(f)$summary)
  select_attempts(summaries, study_config(), kind = job$kind)
}

predict_job <- function(study, job, selected_attempt, frozen) {
  out <- file.path(study, "jobs", "predictions", paste0(job$job_id, ".rds"))
  if (file.exists(out)) return(readRDS(out))
  a <- readRDS(file.path(study, "jobs", "attempts", sprintf("%s-%d.rds", job$job_id, selected_attempt)))
  cm <- readRDS(job$input_file); pr <- load_prediction(cm); pr$species_id <- a$training$species_id
  thin <- if (is.na(frozen$prediction_thin)) 1L else frozen$prediction_thin
  fit <- a$fit
  if (job$kind == "bayesian" && thin > 1L) fit <- thin_fit(fit, thin)
  t0 <- Sys.time()
  prediction <- switch(job$package, occJSDM = predict_occjsdm(fit, pr, a$basis), Hmsc = predict_hmsc(fit, pr, a$basis), gllvm = predict_gllvm(fit, pr, a$basis))
  record <- list(job_id = job$job_id, selected_attempt = selected_attempt, prediction = prediction, prediction_thin = thin,
                 elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"), finished = Sys.time())
  tmp <- paste0(out, ".tmp"); saveRDS(record, tmp); file.rename(tmp, out); record
}

# Keeps every thin-th retained draw of a Bayesian fit for prediction only; selection and diagnostics used all draws.
thin_fit <- function(fit, thin) {
  if (fit$package == "occJSDM") {
    js <- fit$native$results_output$jsdm_output; keep <- seq(1L, dim(js$B0_output)[2], by = thin)
    for (nm in c("B0_output", "B_output", "Bs_output", "sigmab_output", "sigmabs_output", "idx_ls_output")) if (!is.null(js[[nm]])) {
      d <- dim(js[[nm]]); it <- length(d) - 1L
      js[[nm]] <- if (it == 1L) js[[nm]][keep, , drop = FALSE] else if (it == 2L) js[[nm]][, keep, , drop = FALSE] else js[[nm]][, , keep, , drop = FALSE]
    }
    fit$native$results_output$jsdm_output <- js
  } else if (fit$package == "Hmsc") {
    fit$native$postList <- lapply(fit$native$postList, function(ch) ch[seq(1L, length(ch), by = thin)])
  }
  fit$prediction_thin <- thin; fit
}

run_job <- function(study, job, config = readRDS(file.path(study, "config.rds")), frozen = config$frozen, mcmc = config$mcmc,
                    slice_sites = NULL, slice_species = NULL) {
  repeat {
    sel <- select_job(study, job, config)
    if (sel$status == "pending") { for (id in sel$next_attempts) run_attempt(study, job, id, config, frozen, mcmc, slice_sites, slice_species); next }
    break
  }
  status <- switch(sel$status, selected = "pass", unresolved_diagnostics = "unresolved", unstable = "unstable", failed = "no fit", budget_amendment_required = "no fit")
  prediction <- if (!is.na(sel$selected_attempt)) tryCatch(predict_job(study, job, sel$selected_attempt, frozen),
                                                         error = function(e) { status <<- "unsupported"; list(error = conditionMessage(e)) }) else NULL
  result <- list(job_id = job$job_id, status = status, selection = sel, prediction_error = prediction$error, finished = Sys.time())
  writeLines(jsonlite::toJSON(result[c("job_id", "status")], auto_unbox = TRUE), file.path(study, "status", paste0(job$job_id, ".json")))
  result
}
```

- [ ] **Step 4: Run to verify it passes**

Run: same command. Expected: `Job tests passed.` The smoke fits use a 40-site, 4-species slice; `run_job` takes `slice_*` only for tests and the production driver never passes them.

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/spatial-stage2-comparison/jobs.R dev/simstudy/spatial-stage2-comparison/test-jobs.R
git commit -m "Add the stage 2 job driver with attempts, selection, prediction and storage"
```

---

### Task 6: Study driver: prepare, timing, freeze, fit, score

**Files:**
- Create: `dev/simstudy/spatial-stage2-comparison/run.R`
- Create: `dev/simstudy/spatial-stage2-comparison/test-run.R`

**Interfaces:**
- Produces: command-line modes `prepare`, `timing`, `freeze`, `fit`, `score`; `stage2_timing(study, ladder_step)`; `stage2_freeze(study, settings)`; `stage2_fit(study, workers)`; `stage2_score(study)` writing `scores/species.csv`, `scores/aggregate.csv`, `scores/pairs.csv`, `scores/reading.json`, `scores/jobs.csv`, `scores/attempts.csv`.

- [ ] **Step 1: Write the failing test**

```r
# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-run.R
library(testthat)
source("dev/simstudy/spatial-stage2-comparison/run.R")
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
study <- file.path(tempdir(), paste0("stage2-run-", as.integer(Sys.time())))
tiny <- list(initial = c(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L), extended = c(nchain = 2L, nburn = 40L, niter = 80L, nthin = 1L))

test_that("freeze records settings once and fit refuses to run before it", {
  stage2_prepare(study, PILOT, jobs_dir = file.path(study, "jobs-local"))
  expect_error(stage2_fit(study, workers = 1L))   # frozen settings are NA
  stage2_freeze(study, list(support_points_900 = 600L, hmsc_method = "GPP", hmsc_knots = 300L, prediction_thin = 10L), timing = list(note = "test"))
  cfg <- readRDS(file.path(study, "config.rds")); expect_equal(cfg$frozen$support_points_900, 600L)
  expect_error(stage2_freeze(study, list(support_points_900 = 900L, hmsc_method = "Full", hmsc_knots = NA, prediction_thin = 1L), timing = list()))
})

test_that("fit on a subset runs every arm, and score writes every table with statuses carried", {
  j <- job_table(study); subset <- j$job_id[j$dataset_id %in% c("t2-n300-rep1", "t2-n300-rep2")]
  r <- stage2_fit(study, workers = 2L, jobs = subset, mcmc = tiny, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))
  expect_equal(length(list.files(file.path(study, "status"))), 10L)
  stage2_score(study, species_subset = c(4L, 6L, 9L, 12L))
  sp <- read.csv(file.path(study, "scores", "species.csv")); expect_equal(nrow(sp), 10L * 4L)
  expect_true(all(sp$status %in% c("pass", "unresolved", "unstable", "unsupported", "no fit")))
  agg <- read.csv(file.path(study, "scores", "aggregate.csv")); expect_equal(nrow(agg), 10L)
  pairs <- read.csv(file.path(study, "scores", "pairs.csv")); expect_true(all(pairs$verdict == "inconclusive"))   # two communities only
  expect_true(file.exists(file.path(study, "scores", "reading.json")))
  jobs <- read.csv(file.path(study, "scores", "jobs.csv")); expect_true(all(c("status", "attempts", "selected_attempt", "fit_seconds", "prediction_seconds") %in% names(jobs)))
})
cat("Run tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-run.R`
Expected: error, `run.R` not found.

- [ ] **Step 3: Write run.R**

```r
# Study driver for stage 2. Always through run-r:
#   run-r -f run.R --args --mode=prepare --study=DIR --pilot=PILOT_DIR --jobs=VOLUME_DIR
#   run-r -f run.R --args --mode=timing  --study=DIR --step=1        (ladder step 1, 2, ...)
#   run-r -f run.R --args --mode=freeze  --study=DIR --step=1        (freeze the settings of a timing step that passed)
#   run-r -f run.R --args --mode=fit     --study=DIR --workers=8
#   run-r -f run.R --args --mode=score   --study=DIR
source("dev/simstudy/spatial-stage2-comparison/jobs.R")

ladder_settings <- function(config, step) {
  sp <- config$ladder$support_points_900; hm <- config$ladder$hmsc
  grid <- expand.grid(sp = seq_along(sp), hm = seq_along(hm))
  grid <- grid[order(grid$sp + grid$hm, grid$sp), ]
  stopifnot(step >= 1L, step <= nrow(grid))
  g <- grid[step, ]
  list(support_points_900 = sp[g$sp], hmsc_method = hm[[g$hm]][1], hmsc_knots = as.integer(hm[[g$hm]][2]), prediction_thin = 10L)
}

stage2_timing <- function(study, step = 1L, workers = 5L) {
  cfg <- readRDS(file.path(study, "config.rds")); frozen <- ladder_settings(cfg, step)
  jobs <- job_table(study); pick <- jobs[jobs$n == 900L & jobs$replicate == 1L, ]   # 2 truths x 5 arms = 10 jobs
  pick$job_id <- paste0(pick$job_id, "--timing", step)
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(pick)), function(i) {
    r <- run_job(study, pick[i, ], config = cfg, frozen = frozen)
    a <- list.files(file.path(study, "jobs", "attempts"), pattern = paste0("^", pick$job_id[i], "-"), full.names = TRUE)
    fit_secs <- sum(vapply(a, function(f) readRDS(f)$summary$elapsed_seconds, numeric(1)))
    p <- file.path(study, "jobs", "predictions", paste0(pick$job_id[i], ".rds"))
    data.frame(job_id = pick$job_id[i], key = pick$key[i], truth = pick$truth[i], status = r$status, attempts = length(a), fit_seconds = fit_secs,
               prediction_seconds = if (file.exists(p)) readRDS(p)$elapsed_seconds else NA_real_, stringsAsFactors = FALSE)
  }, mc.cores = workers, mc.preschedule = FALSE)
  timing <- do.call(rbind, res)
  # Upper bound: every job priced at its arm's 900-site initial cost; Bayesian arms times the extension factor; plus prediction.
  per_arm <- aggregate(cbind(fit_seconds, prediction_seconds) ~ key, timing, max)
  per_arm$bayesian <- !grepl("^gllvm", per_arm$key)
  per_arm$job_upper <- ifelse(per_arm$bayesian, per_arm$fit_seconds * cfg$extension_factor, per_arm$fit_seconds * 2) + per_arm$prediction_seconds
  n_jobs <- table(jobs$key)
  cpu_hours <- sum(per_arm$job_upper * as.numeric(n_jobs[per_arm$key])) / 3600
  projection <- list(step = step, settings = frozen, timing = timing, per_arm = per_arm, projected_cpu_hours_upper = cpu_hours,
                     projected_wall_hours_upper = cpu_hours / cfg$workers, limit_wall_hours = cfg$wall_hours,
                     passes = cpu_hours / cfg$workers < cfg$wall_hours, wall_seconds = as.numeric(Sys.time() - t0, units = "secs"))
  writeLines(jsonlite::toJSON(projection, auto_unbox = TRUE, pretty = TRUE, digits = 4, na = "null"), file.path(study, sprintf("timing-step%d.json", step)))
  projection
}

stage2_freeze <- function(study, settings, timing) {
  cfg <- readRDS(file.path(study, "config.rds"))
  stopifnot(all(is.na(unlist(cfg$frozen))), setequal(names(settings), names(cfg$frozen)))
  cfg$frozen <- settings; cfg$frozen_at <- Sys.time(); cfg$frozen_from <- timing
  saveRDS(cfg, file.path(study, "config.rds")); invisible(cfg)
}

stage2_fit <- function(study, workers = 8L, jobs = NULL, mcmc = NULL, slice_sites = NULL, slice_species = NULL) {
  cfg <- readRDS(file.path(study, "config.rds"))
  if (any(is.na(unlist(cfg$frozen)[c("support_points_900", "hmsc_method", "prediction_thin")]))) stop("Settings are not frozen; run timing and freeze first")
  all_jobs <- job_table(study); if (!is.null(jobs)) all_jobs <- all_jobs[all_jobs$job_id %in% jobs, ]
  if (is.null(mcmc)) mcmc <- cfg$mcmc
  res <- parallel::mclapply(seq_len(nrow(all_jobs)), function(i) tryCatch(run_job(study, all_jobs[i, ], cfg, cfg$frozen, mcmc, slice_sites, slice_species)$status,
    error = function(e) { writeLines(conditionMessage(e), file.path(study, "status", paste0(all_jobs$job_id[i], ".error"))); "error" }),
    mc.cores = workers, mc.preschedule = FALSE)
  table(unlist(res))
}

stage2_score <- function(study, species_subset = NULL) {
  cfg <- readRDS(file.path(study, "config.rds")); jobs <- job_table(study)
  done <- jobs[file.exists(file.path(study, "status", paste0(jobs$job_id, ".json"))), ]
  status_of <- function(id) jsonlite::fromJSON(file.path(study, "status", paste0(id, ".json")))$status
  rows <- list(); job_rows <- list()
  for (ds in unique(done$dataset_id)) {
    cm <- readRDS(done$input_file[done$dataset_id == ds][1]); truth <- load_truth(cm)
    if (!is.null(species_subset)) for (nm in c("psi_lattice", "field_lattice", "eta_lattice", "psi_training", "field_training")) truth[[nm]] <- truth[[nm]][, species_subset, drop = FALSE]
    if (!is.null(species_subset)) { truth$occupancy <- truth$occupancy[species_subset]; truth$species <- truth$species[species_subset]; truth$spatial_sd <- truth$spatial_sd[species_subset] }
    oracle <- load_oracle_reference(done$oracle_file[done$dataset_id == ds][1])
    if (!is.null(species_subset)) oracle$probability_mean <- oracle$probability_mean[, species_subset, drop = FALSE]
    load_pred <- function(id) { p <- file.path(study, "jobs", "predictions", paste0(id, ".rds")); if (file.exists(p)) readRDS(p) else NULL }
    env_id <- done$job_id[done$dataset_id == ds & done$key == "occJSDM/env"]; env <- load_pred(env_id)
    for (i in which(done$dataset_id == ds)) {
      j <- done[i, ]; st <- status_of(j$job_id); p <- load_pred(j$job_id)
      rows[[length(rows) + 1L]] <- score_arm(if (is.null(p) || is.null(env)) NULL else p$prediction, if (is.null(env)) NULL else env$prediction, oracle, truth, ds, j$package, j$arm,
                                             if (is.null(env)) "unsupported" else st)
      a <- list.files(file.path(study, "jobs", "attempts"), pattern = paste0("^", j$job_id, "-[0-9]+\\.rds$"), full.names = TRUE)
      sel <- select_job(study, j, cfg)
      job_rows[[length(job_rows) + 1L]] <- data.frame(job_id = j$job_id, dataset_id = ds, cell = j$cell, key = j$key, status = st, attempts = length(a),
        selected_attempt = sel$selected_attempt, fit_seconds = sum(vapply(a, function(f) readRDS(f)$summary$elapsed_seconds, numeric(1))),
        prediction_seconds = if (is.null(p)) NA_real_ else p$elapsed_seconds, oracle_status = oracle$status, stringsAsFactors = FALSE)
    }
  }
  species <- do.call(rbind, rows); write.csv(species, file.path(study, "scores", "species.csv"), row.names = FALSE)
  agg <- aggregate_arm(species, cfg$gate_occupancy); write.csv(agg, file.path(study, "scores", "aggregate.csv"), row.names = FALSE)
  pairs <- paired_comparisons(agg, cfg); write.csv(pairs, file.path(study, "scores", "pairs.csv"), row.names = FALSE)
  max_share <- suppressWarnings(max(aggregate(headroom_share ~ cell + key, agg[agg$status == "pass", ], median)$headroom_share, na.rm = TRUE))
  writeLines(jsonlite::toJSON(development_reading(pairs, max_share), auto_unbox = TRUE, pretty = TRUE, na = "null"), file.path(study, "scores", "reading.json"))
  write.csv(do.call(rbind, job_rows), file.path(study, "scores", "jobs.csv"), row.names = FALSE)
  invisible(pairs)
}

if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  get <- function(key, default = NULL) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), args, value = TRUE)); if (length(v)) v else default }
  mode <- get("mode"); study <- get("study")
  switch(mode,
    prepare = stage2_prepare(study, get("pilot"), jobs_dir = get("jobs")),
    timing = print(stage2_timing(study, step = as.integer(get("step", "1")))[c("settings", "projected_cpu_hours_upper", "projected_wall_hours_upper", "passes")]),
    freeze = { s <- as.integer(get("step")); t <- jsonlite::fromJSON(file.path(study, sprintf("timing-step%d.json", s))); stopifnot(isTRUE(t$passes))
               stage2_freeze(study, ladder_settings(readRDS(file.path(study, "config.rds")), s), timing = t[c("step", "projected_wall_hours_upper")]) },
    fit = print(stage2_fit(study, workers = as.integer(get("workers", "8")))),
    score = print(stage2_score(study)),
    stop("Unknown mode ", mode))
  cat("Completed mode", mode, "at", format(Sys.time()), "\n")
}
```

- [ ] **Step 4: Run to verify it passes**

Run: same command. Expected: `Run tests passed.`

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/spatial-stage2-comparison/run.R dev/simstudy/spatial-stage2-comparison/test-run.R
git commit -m "Add the stage 2 study driver with timing ladder, freeze, fit and score"
```

---

### Task 7: Prepare, timing ladder, freeze and dispatch

**Files:**
- Create (ignored): `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-stage2-20261011/`
- Create (volume): `/Volumes/occJSDM Research/occJSDM/spatial-stage2-20261011/jobs/`

- [ ] **Step 1: Run every test file**

```bash
P=dev/simstudy/spatial-stage2-comparison
for t in runtime data adapters score jobs run; do $P/run-r -f $P/test-$t.R || exit 1; done
```
Expected: six "passed" lines.

- [ ] **Step 2: Prepare**

```bash
P=dev/simstudy/spatial-stage2-comparison
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-stage2-20261011
PILOT=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010
JOBS="/Volumes/occJSDM Research/occJSDM/spatial-stage2-20261011/jobs"
$P/run-r -f $P/run.R --args --mode=prepare --study=$STUDY --pilot=$PILOT --jobs="$JOBS"
wc -l $STUDY/datasets.csv
```
Expected: 41 lines (header plus 40).

- [ ] **Step 2b: Micro-timing estimate (minutes, not hours)**

Before the ladder, give Doug a run-time estimate cheap enough to make before the laptop is handed over. Add `--mode=microtiming` to `run.R`: for each of the five arms on `t1-r0.08-n900-rep1` and `t2-n900-rep1`, fit with one chain, 50 warm-up and 150 retained iterations (gllvm: one start with `maxit = 200`) under ladder step 1 settings, time the fit and a 25-site prediction block with every retained draw, and extrapolate: per-iteration cost times the full initial schedule (4 chains, 7,000 iterations) times the extension factor 4.43 for Bayesian arms, plus per-block prediction cost times 64 blocks at the ladder's prediction thinning. Price 300-site jobs at the 900-site cost (upper bound). Write `microtiming.json` with per-arm projected hours, the total upper bound in worker-hours and in wall-hours at eight workers, and the same for ladder steps 2 to 4 by rescaling the gp cost with the support-point count squared and the Hmsc cost by the knot count over 900 (recorded as rescaling assumptions, not measurements). Report the table to Doug before anything else runs. The micro fits are discarded, not study attempts.

```bash
$P/run-r -f $P/run.R --args --mode=microtiming --study=$STUDY
cat $STUDY/microtiming.json
```
Expected: about 10 to 30 minutes of wall time, dominated by the Hmsc six-factor and gp micro fits.

- [ ] **Step 3: Timing ladder**

```bash
nohup $P/run-r -f $P/run.R --args --mode=timing --study=$STUDY --step=1 > $STUDY/logs/timing-step1.log 2>&1 &
```
Step 1 is 900 support points and Hmsc Full. When the log ends with `Completed mode timing`, read `timing-step1.json`. If `passes` is true, freeze step 1. If not, run step 2 (600 support points, Full), then step 3 (900, GPP 300 knots), and so on down the ladder, freezing the first step that passes. If no step passes, stop and report the projections to Doug; do not reduce cells, arms, budgets or diagnostics.

```bash
$P/run-r -f $P/run.R --args --mode=freeze --study=$STUDY --step=1
```

- [ ] **Step 4: Dispatch**

```bash
nohup caffeinate -i $P/run-r -f $P/run.R --args --mode=fit --study=$STUDY --workers=8 > $STUDY/logs/fit.log 2>&1 &
```
Monitor with `ls $STUDY/status/*.json | wc -l` (target 200) and `ls $STUDY/status/*.error`. Relaunch after any interruption; completed attempts are reused. Record launch time, PID and power state in `$STUDY/logs/fit-notes.txt`.

- [ ] **Step 5: Score**

```bash
$P/run-r -f $P/run.R --args --mode=score --study=$STUDY
cat $STUDY/scores/reading.json; cat $STUDY/scores/pairs.csv
```
Expected: `species.csv` with 200 × 12 = 2,400 rows, `aggregate.csv` with 200 rows, `pairs.csv` with 4 cells × 5 comparisons = 20 rows, `reading.json` naming one of the five readings.

- [ ] **Step 6: Record and commit**

Add one revision-record paragraph to the spec: frozen settings and the ladder step, timing numbers, launch and completion times, status counts by arm. Commit.

---

### Task 8: Verification, plots, report and handover

**Files:**
- Create: `dev/simstudy/spatial-stage2-comparison/verify.R`
- Create: `dev/simstudy/spatial-stage2-comparison/plot.R`
- Create: `dev/simstudy/spatial-stage2-comparison/results/` (compact CSVs, `headroom-share.png`, `verification.csv`, `runtime.json`)
- Create: `dev/simstudy/spatial-stage2-comparison/REPORT.md`, `README.md`
- Modify: `TODO.md` (one line)

- [ ] **Step 1: Write verify.R**

Independent recomputation, written without reference to `score.R`'s functions: for every cell and comparison, recompute the paired estimate, interval and verdict from `species.csv` (aggregate species within occupancy level, pool 25, 50 and 75, pair within dataset, require both `pass`, apply delta 0.10 and minimum pairing 6) and compare with `pairs.csv`; recompute the development reading from `pairs.csv`; check `datasets.csv` hashes against the pilot manifest and the oracle files; check that `config.rds`'s frozen settings were set after the timing step's JSON and before the earliest attempt file in `jobs/attempts`; check 200 status files after excluding the `--timing` jobs' files, no `.error` files, attempts at most 400 excluding timing attempts and per job at most 2 for Bayesian and 6 for point; check every prediction record passed `validate_prediction`; write `results/verification.csv` and exit non-zero on any failure.

- [ ] **Step 2: Write plot.R**

Copy `scores/*.csv`, `scores/reading.json`, `timing-step*.json` and a `runtime.json` (R version, package versions, frozen settings, provisioning route) into `results/`. Draw `results/headroom-share.png`: one panel per cell, arms on the x axis in roster order, one point per community (jittered), the median as a bar, the oracle at 1 and the env arm at 0 as reference lines, status shown by point shape (pass filled, others open). Draw `results/pairs.png`: the 20 paired estimates with intervals, the delta band at plus and minus 0.10 shaded.

- [ ] **Step 3: Run verify and plot, inspect the PNGs**

```bash
$P/run-r -f $P/verify.R --args --study=$STUDY
$P/run-r -f $P/plot.R --args --study=$STUDY
```
Expected: `Verification passed.`; open both PNGs with the Read tool and check labels and clipping.

- [ ] **Step 4: Write REPORT.md**

Pedagogical, as Doug asked for the pilot report, one line per paragraph, every number from `results/`. Sections: the short answer (the development reading and what it rests on); the idea (what an arm is, what the headroom share measures and why the oracle is the ceiling, why truth 1 pools); design as run (datasets, arms, frozen settings and the ladder step with its timing, budgets, runtime); results (the headroom-share figure; one bullet per cell and arm with median share, pattern, slope where defined, status counts; the 20 paired comparisons with verdicts; the cross-truth misspecification table as bullets); diagnostics (statuses by arm, extension counts, gllvm instability, how unresolved cells affected the intervals); the development reading and what stage 3 would target; the verbatim mechanism paragraph from the pilot spec; limits (oracle knows intercepts and slopes; probit arms unscored on slope; exponential against squared-exponential kernels; GPP if used; ten communities); decisions for Doug; files; glossary.

- [ ] **Step 5: README, TODO line, commit, review**

README: purpose, results pointer, the exact reproduction commands from Task 7, file list. TODO.md: one line under the spatial section with the reading and the report path. Commit, then request a whole-branch review against `main` with superpowers:requesting-code-review and fix findings. Hand Doug `REPORT.md` and the two PNGs; he reviews rendered output.

---

## Deviations expected during execution

- If Route B provisioning is used, Hmsc and gllvm versions differ from the first comparison's; the config and report record the realised versions.
- If the timing ladder freezes GPP for Hmsc, prediction uses Hmsc's own `predictLatentFactor` rather than the first comparison's certified Full conditional; the conditioning string records it and the report says the approximation differs.
- Prediction thinning (default 10 in the ladder) departs from the first comparison's unthinned prediction; selection and diagnostics still use every draw. Recorded in the config.
