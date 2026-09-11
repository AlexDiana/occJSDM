# Portable reproduction of run_binary.R (six repeats) and run_binary30.R
# (30 repeats). Use a PRECOMPILED package snapshot; this runner never builds
# or patches it. The source path selects baseline/intermediate/blocked code.
#
# Rscript dev/simstudy/validate_spatial_binary.R --source=/path/to/package \
#   --out=/path/to/new-output --label=blocked30-range6 --grid-index=6 \
#   --repeats=30 --burn=500 --iter=500 --chains=2
#
# Six repeats reproduce the original single RNG stream exactly. More than
# six reuse its coordinates and GP innovations, add environmental values
# using a separate stream, and use a separate outcome stream. Outcomes are
# not nested between the six- and 30-repeat designs. Both use 80 locations,
# eight species, unit field SD, fixed coefficients and default model priors.
# --generate-only=true saves data and provenance without fitting.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Required: --source=PATH --out=PATH --grid-index=1..10 --burn=N --iter=N\n",
      "Optional: --label=NAME --repeats=6 --chains=2 --knots=79 --generate-only=false\n",
      "The source must be precompiled. The output directory must be empty.\n", sep = "")
  quit(status = 0)
}
option <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", name)
  if (!length(hit)) {
    if (is.null(default)) stop("Missing --", name, "=...")
    return(default)
  }
  substring(hit, nchar(name) + 4L)
}
known <- c("source", "out", "label", "grid-index", "burn", "iter", "chains",
            "knots", "repeats", "generate-only")
if (any(!sub("^--([^=]+)=.*$", "\\1", args) %in% known)) stop("Unknown option; use --help")
integer_option <- function(name, default = NULL, minimum = 1L) {
  value <- suppressWarnings(as.numeric(option(name, default)))
  if (length(value) != 1L || !is.finite(value) || value != floor(value) || value < minimum)
    stop("Invalid integer option: ", name)
  as.integer(value)
}
source_dir <- normalizePath(option("source"), mustWork = TRUE)
out_dir <- option("out")
label <- option("label", "spatial-binary")
grid_index <- integer_option("grid-index")
nburn <- integer_option("burn", minimum = 0L)
niter <- integer_option("iter")
chains <- integer_option("chains", "2")
knots <- integer_option("knots", "79")
repeats <- integer_option("repeats", "6", minimum = 6L)
generate_only <- option("generate-only", "false")
stopifnot(grid_index <= 10L, knots < 80L, generate_only %in% c("true", "false"))
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
  stop("Output directory must be empty: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
pkgload::load_all(source_dir, compile = FALSE, quiet = TRUE, helpers = FALSE)
RcppParallel::setThreadOptions(numThreads = 1)

# Hash source and the actual loaded DLL. A snapshot without its own Git
# metadata must not inherit an unrelated enclosing repository's revision.
files <- c(file.path(source_dir, c("DESCRIPTION", "NAMESPACE")),
           list.files(file.path(source_dir, "R"), pattern = "\\.R$", full.names = TRUE),
           list.files(file.path(source_dir, "src"), pattern = "\\.(cpp|h|so|dll|dylib)$|^Makevars", full.names = TRUE))
dll <- getLoadedDLLs()[["occJSDM"]][["path"]]
files <- c(files, dll)
hashes <- tools::md5sum(files)
manifest <- data.frame(file = ifelse(startsWith(names(hashes), paste0(source_dir, "/")),
  substring(names(hashes), nchar(source_dir) + 2L), paste0("loaded-dll/", basename(names(hashes)))),
  md5 = unname(hashes), stringsAsFactors = FALSE)
git <- function(command) {
  if (!file.exists(file.path(source_dir, ".git"))) return(NULL)
  suppressWarnings(system2("git", c("-C", shQuote(source_dir), command), stdout = TRUE, stderr = TRUE))
}
snapshot_record <- file.path(source_dir, "source-provenance.rds")
runner <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")])
mcmc <- list(nchain = chains, nburn = nburn, niter = niter, nthin = 1L)
parameters <- list(n_factors = 0L, n_lattrait = 0L, n_supportpoints = knots)
provenance <- list(label = label, source = source_dir, loaded_dll = dll,
  head = git(c("rev-parse", "HEAD")), diff = git("diff"), source_hashes = manifest,
  snapshot_provenance = if (file.exists(snapshot_record)) readRDS(snapshot_record) else NULL,
  runner_hash = tools::md5sum(runner),
  command = commandArgs(), started = Sys.time(), session = sessionInfo(),
  settings = list(grid_index = grid_index, repeats = repeats, parameters = parameters,
                   mcmc = mcmc, priors = "Unmodified defaults of the selected source"))
write.csv(manifest, file.path(out_dir, "source-hashes.csv"), row.names = FALSE)
saveRDS(provenance, file.path(out_dir, "provenance.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))

data_seed <- 817000L + grid_index
fit_seed <- 918000L + grid_index
nlocations <- 80L; S <- 8L; n <- nlocations * repeats
location <- rep(seq_len(nlocations), each = repeats)
range_grid <- seq(.01, .3, length.out = 10)
true_range <- range_grid[grid_index]
set.seed(data_seed)
coordinates <- matrix(runif(nlocations * 2), nlocations, 2)
raw_coordinates <- coordinates[location, ]
standardized_coordinates <- scale(raw_coordinates)
unique_standardized <- standardized_coordinates[seq(1, n, by = repeats), ]
if (repeats == 6L) {
  X_raw <- rnorm(n)
  GP_innovations <- matrix(rnorm(nlocations * S), nlocations, S)
} else {
  X_first6 <- matrix(rnorm(nlocations * 6L), nlocations, 6L, byrow = TRUE)
  GP_innovations <- matrix(rnorm(nlocations * S), nlocations, S)
  set.seed(827000L + grid_index)
  X_extra <- matrix(rnorm(nlocations * (repeats - 6L)), nlocations, repeats - 6L)
  X_raw <- as.vector(t(cbind(X_first6, X_extra)))
}
X <- as.numeric(scale(X_raw))
d2 <- as.matrix(dist(unique_standardized))^2
K <- exp(-d2 / (2 * true_range^2)) + diag(1e-10, nlocations)
field_locations <- t(chol(K)) %*% GP_innovations
field <- field_locations[location, ]
B0 <- c(-.4, .3, -.2, .5, .2, -.3, .1, -.1)
B <- c(.3, -.3, .2, -.2, .4, -.4, .1, -.1)
eta <- matrix(B0, n, S, byrow = TRUE) + X %o% B + field
psi <- plogis(eta)
if (repeats > 6L) set.seed(837000L + grid_index)
y <- matrix(rbinom(n * S, 1, psi), n, S)
colnames(y) <- paste0("species", seq_len(S))
data <- list(info = data.frame(Site = seq_len(n), environment = X_raw,
  longitude = raw_coordinates[, 1], latitude = raw_coordinates[, 2]), OTU = y, traits = NULL)
truth <- list(range = true_range, grid_index = grid_index, field_sd = 1, field = field,
  field_locations = field_locations, eta = eta, psi = psi, B0 = B0, B = B,
  standardized_coordinates = standardized_coordinates, standardized_environment = X,
  location = location, data_seed = data_seed, fit_seed = fit_seed,
  environment_seed = if (repeats > 6L) 827000L + grid_index else NULL,
  outcome_seed = if (repeats > 6L) 837000L + grid_index else NULL,
  GP_innovations = GP_innovations, raw_coordinates = coordinates, repeats = repeats)
saveRDS(list(data = data, truth = truth), file.path(out_dir, "data-truth.rds"))
provenance$seeds <- truth[c("data_seed", "fit_seed", "environment_seed", "outcome_seed")]
saveRDS(provenance, file.path(out_dir, "provenance.rds"))
if (generate_only == "true") quit(status = 0)

set.seed(fit_seed)
fit_warnings <- character()
elapsed <- system.time(fit <- withCallingHandlers(suppressMessages(runOccJSDM(
  data, occCovariates = "environment", spatCovariates = c("longitude", "latitude"),
  listParams = parameters, MCMCparams = mcmc)), warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }))
saveRDS(fit, file.path(out_dir, "fit.rds"))
writeLines(fit_warnings, file.path(out_dir, "fit-warnings.txt"))
stopifnot(max(abs(unname(fit$Xs) - unname(standardized_coordinates))) < 1e-12,
          max(abs(as.vector(fit$X_psi) - X)) < 1e-12,
          fit$infos$ps == knots, fit$infos$n_factors == 0L)
if (!is.null(fit$infos$gt)) stopifnot(fit$infos$gt == 0L)
jsdm <- fit$results_output$jsdm_output
idx <- jsdm$idx_ls_output
bases <- suppressMessages(precomputeSORmatrices(fit$infos$l_s_grid, fit$infos$list_Xs))
field_chain <- psi_chain <- array(0, c(n, S, chains))
for (ch in seq_len(chains)) for (it in seq_len(niter)) {
  se <- KsBproduct(matrix(bases$Ks_all[, , idx[it, ch]], nrow = n),
                   matrix(jsdm$Bs_output[, , it, ch], nrow = knots), fit$infos$list_Xs$Xs_centers)
  eta_draw <- sweep(se + X %o% jsdm$B_output[1, , it, ch], 2, jsdm$B0_output[, it, ch], "+")
  field_chain[, , ch] <- field_chain[, , ch] + se / niter
  psi_chain[, , ch] <- psi_chain[, , ch] + plogis(eta_draw) / niter
}
field_mean <- apply(field_chain, c(1, 2), mean)
psi_mean <- apply(psi_chain, c(1, 2), mean)
if (!is.null(fit$results_output$psi_output))
  stopifnot(max(abs(psi_mean - fit$results_output$psi_output)) < 1e-12)
score <- function(estimate, truth) c(bias = mean(estimate - truth),
  mae = mean(abs(estimate - truth)), rmse = sqrt(mean((estimate - truth)^2)),
  correlation = cor(as.vector(estimate), as.vector(truth)),
  estimate_mean = mean(estimate), truth_mean = mean(truth))
metrics <- rbind(field = score(field_mean, field), probability = score(psi_mean, psi))
for (term in c("field", "probability")) for (ch in seq_len(chains)) {
  row <- if (term == "field") score(field_chain[, , ch], field) else score(psi_chain[, , ch], psi)
  metrics <- rbind(metrics, row)
  rownames(metrics)[nrow(metrics)] <- paste0(term, "_chain", ch)
}
stratum <- cut(psi, breaks = seq(0, 1, .2), include.lowest = TRUE)
strata <- do.call(rbind, lapply(levels(stratum), function(s) {
  i <- stratum == s
  data.frame(stratum = s, n = sum(i), t(score(psi_mean[i], psi[i])))
}))
range_distribution <- do.call(rbind, lapply(seq_len(chains), function(ch) {
  counts <- tabulate(idx[, ch], nbins = length(range_grid))
  data.frame(chain = ch, range = fit$infos$l_s_grid, count = counts, proportion = counts / niter)
}))
range_summary <- do.call(rbind, lapply(seq_len(chains), function(ch) {
  v <- fit$infos$l_s_grid[idx[, ch]]
  data.frame(chain = ch, truth = true_range, mean = mean(v), median = median(v),
    lower = unname(quantile(v, .025)), upper = unname(quantile(v, .975)),
    upper_boundary = mean(idx[, ch] == length(range_grid)), transitions = sum(diff(idx[, ch]) != 0))
}))
write.csv(metrics, file.path(out_dir, "metrics.csv"))
write.csv(strata, file.path(out_dir, "probability-strata.csv"), row.names = FALSE)
write.csv(range_distribution, file.path(out_dir, "range-distribution.csv"), row.names = FALSE)
write.csv(range_summary, file.path(out_dir, "range-summary.csv"), row.names = FALSE)
saveRDS(list(field_mean = field_mean, field_chain = field_chain, psi_mean = psi_mean,
  psi_chain = psi_chain, metrics = metrics, strata = strata, range_distribution = range_distribution,
  range_summary = range_summary, elapsed = elapsed, nburn = nburn, niter = niter),
  file.path(out_dir, "validation.rds"))
provenance$finished <- Sys.time()
provenance$elapsed <- elapsed
saveRDS(provenance, file.path(out_dir, "provenance.rds"))
print(range_summary); print(metrics); print(strata)
