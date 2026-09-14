# Portable reproduction of the continuous spatial range validation.
# Use a PRECOMPILED package snapshot; this runner never builds or patches it.
# The source path selects baseline/intermediate/blocked code. The label only
# identifies the run and does not select or change the sampler.
#
# Original 100-location design (cell 1/2/3 selects range grid index 4/6/8):
# Rscript dev/simstudy/validate_spatial_continuous.R --source=/path/to/package \
#   --out=/path/to/new-output --label=blocked --cell=2 --burn=400 --iter=400
# Informative repeated-coordinate design, sharing RNG streams across cells:
# Rscript dev/simstudy/validate_spatial_continuous.R --source=/path/to/package \
#   --out=/path/to/new-output --label=blocked-repeat10 --cell=2 --burn=500 \
#   --iter=500 --repeats=10 --tau=1 --shared-seed=true --family=repeat10
#
# Each repeated row has its own Site identity; coordinates, environmental
# value and GP field are shared within location, while measurement noise is
# independent. Coordinates are standardized over all observed rows before
# generating the GP. Default priors, unit field SD, 12 species and 99 knots
# match the original runner. --generate-only=true skips the fit.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Required: --source=PATH --out=PATH --cell=1..3 --burn=N --iter=N\n",
      "Optional: --label=NAME --family=NAME --repeats=1 --tau=0.1\n",
      "          --shared-seed=false --chains=2 --knots=99 --generate-only=false\n",
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
known <- c("source", "out", "label", "family", "cell", "burn", "iter", "chains",
           "knots", "repeats", "tau", "shared-seed", "generate-only")
if (any(!sub("^--([^=]+)=.*$", "\\1", args) %in% known)) stop("Unknown option; use --help")
integer_option <- function(name, default = NULL, minimum = 1L) {
  value <- suppressWarnings(as.numeric(option(name, default)))
  if (length(value) != 1L || !is.finite(value) || value != floor(value) || value < minimum)
    stop("Invalid integer option: ", name)
  as.integer(value)
}
source_dir <- normalizePath(option("source"), mustWork = TRUE)
out_dir <- option("out")
label <- option("label", "spatial-continuous")
family <- option("family", "")
cell <- integer_option("cell")
burn <- integer_option("burn", minimum = 0L)
retained <- integer_option("iter")
chains <- integer_option("chains", "2")
knots <- integer_option("knots", "99")
repeats <- integer_option("repeats", "1")
noise_sd <- suppressWarnings(as.numeric(option("tau", ".1")))
shared_seed <- option("shared-seed", "false")
generate_only <- option("generate-only", "false")
stopifnot(cell <= 3L, knots < 100L, length(noise_sd) == 1L,
          is.finite(noise_sd), noise_sd > 0,
          shared_seed %in% c("true", "false"), generate_only %in% c("true", "false"))
shared_seed <- shared_seed == "true"
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
  stop("Output directory must be empty: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
pkgload::load_all(source_dir, compile = FALSE, quiet = TRUE, helpers = FALSE)
RcppParallel::setThreadOptions(numThreads = 1)

# Hash source and the actual loaded DLL. Preserve a supplied snapshot's
# provenance rather than attributing it to an unrelated enclosing Git repo.
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
mcmc <- list(nchain = chains, nburn = burn, niter = retained, nthin = 1L)
parameters <- list(n_factors = 0L, n_lattrait = 0L, n_supportpoints = knots)
provenance <- list(label = label, source = source_dir, loaded_dll = dll,
  head = git(c("rev-parse", "HEAD")), diff = git("diff"), source_hashes = manifest,
  snapshot_provenance = if (file.exists(snapshot_record)) readRDS(snapshot_record) else NULL,
  runner_hash = tools::md5sum(runner),
  command = commandArgs(), started = Sys.time(), session = sessionInfo(),
  settings = list(cell = cell, family = family, repeats = repeats, tau = noise_sd,
    shared_seed = shared_seed, parameters = parameters, mcmc = mcmc,
    priors = "Unmodified defaults of the selected source"))
write.csv(manifest, file.path(out_dir, "source-hashes.csv"), row.names = FALSE)
saveRDS(provenance, file.path(out_dir, "provenance.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))

grid <- seq(.01, .3, length.out = 10)
true_range <- grid[c(4, 6, 8)[cell]]
seed_cell <- if (shared_seed) 1L else cell
dataset_seed <- 202609110L + seed_cell
fit_seed <- 202619110L + seed_cell
cell_tag <- if (nzchar(family)) paste(family, cell, sep = "-") else as.character(cell)
set.seed(dataset_seed)
n_locations <- 100L; S <- 12L
index <- rep(seq_len(n_locations), each = repeats)
n <- length(index)
first <- seq(1, n, by = repeats)
raw_coords <- matrix(runif(n_locations * 2), n_locations, 2)
coords <- scale(raw_coords[index, , drop = FALSE])
env <- as.numeric(scale(rnorm(n_locations)[index]))
intercept <- seq(-.8, .8, length.out = S)
beta <- rep(c(-.5, 0, .5), length.out = S)
K <- K2(coords[first, , drop = FALSE], coords[first, , drop = FALSE], 1, true_range) +
  diag(1e-5, n_locations)
spatial_locations <- t(chol(K)) %*% matrix(rnorm(n_locations * S), n_locations, S)
spatial <- spatial_locations[index, , drop = FALSE]
eta <- matrix(intercept, n, S, byrow = TRUE) + outer(env, beta) + spatial
y <- eta + matrix(rnorm(n * S, sd = noise_sd), n, S)
colnames(y) <- paste0("sp", seq_len(S))
data <- list(info = data.frame(Site = seq_len(n), env = env, x = coords[, 1], y = coords[, 2]), OTU = y)
truth <- list(spatial = spatial, eta = eta, intercept = intercept, beta = beta, tau = noise_sd,
  coordinates = coords, range = true_range, n = n, S = S,
  location_index = index, n_locations = n_locations, repeats = repeats)
input <- list(data = data, truth = truth, dataset_seed = dataset_seed, fit_seed = fit_seed)
saveRDS(input, file.path(out_dir, "data-truth.rds"))
provenance$seeds <- list(dataset_seed = dataset_seed, fit_seed = fit_seed)
saveRDS(provenance, file.path(out_dir, "provenance.rds"))
if (generate_only == "true") quit(status = 0)

set.seed(fit_seed)
fit_warnings <- character()
elapsed <- system.time(fit <- withCallingHandlers(suppressMessages(runOccJSDM(
  input$data, listParams = parameters, occCovariates = "env", spatCovariates = c("x", "y"),
  MCMCparams = mcmc)), warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }))
saveRDS(fit, file.path(out_dir, "fit.rds"))
writeLines(fit_warnings, file.path(out_dir, "fit-warnings.txt"))
stopifnot(max(abs(fit$Xs - truth$coordinates)) < 1e-12,
          fit$infos$ps == knots, fit$infos$n_factors == 0L)
ro <- fit$results_output$jsdm_output
summaries <- suppressMessages(precomputeSORmatrices(fit$infos$l_s_grid, fit$infos$list_Xs))
field <- matrix(0, truth$n, truth$S)
eta <- field
chain_fields <- array(0, c(dim(field), chains))
for (ch in seq_len(chains)) for (it in seq_len(retained)) {
  index <- ro$idx_ls_output[it, ch]
  Ks <- matrix(summaries$Ks_all[, , index], nrow = nrow(field))
  field_draw <- KsBproduct(Ks, matrix(ro$Bs_output[, , it, ch], nrow = knots),
                          fit$infos$list_Xs$Xs_centers)
  field <- field + field_draw / (retained * chains)
  chain_fields[, , ch] <- chain_fields[, , ch] + field_draw / retained
  eta <- eta + (matrix(ro$B0_output[, it, ch], nrow(field), ncol(field), byrow = TRUE) +
    fit$X_psi %*% matrix(ro$B_output[, , it, ch], nrow = 1) + field_draw) / (retained * chains)
}
metric <- function(estimated, truth) c(bias = mean(estimated - truth),
  rmse = sqrt(mean((estimated - truth)^2)), cor = cor(c(estimated), c(truth)),
  centered_rmse = sqrt(mean((scale(estimated, scale = FALSE) - scale(truth, scale = FALSE))^2)))
range_draws <- matrix(grid[ro$idx_ls_output], retained, chains)
range_table <- table(factor(ro$idx_ls_output, levels = seq_along(grid)),
                     rep(seq_len(chains), each = retained))
metrics <- list(version = label, cell = cell, family = family, cell_tag = cell_tag,
  true_range = true_range, elapsed_seconds = unname(elapsed[["elapsed"]]),
  spatial = metric(field, truth$spatial), eta = metric(eta, truth$eta),
  range_mean = colMeans(range_draws), range_median = apply(range_draws, 2, median),
  range_table = range_table, range_quantiles = apply(range_draws, 2, quantile, c(.025, .5, .975)),
  range_transitions = apply(range_draws, 2, function(x) sum(diff(x) != 0)),
  sigmabs = mean(ro$sigmabs_output), tau = mean(ro$tau_output),
  intercept_bias = mean(apply(ro$B0_output, 1, mean) - truth$intercept),
  chain_field_rmse = vapply(seq_len(chains), function(ch)
    metric(chain_fields[, , ch], truth$spatial)["rmse"], numeric(1)))
saveRDS(list(metrics = metrics, field = field, eta = eta, range_draws = range_draws,
             chain_fields = chain_fields), file.path(out_dir, "validation.rds"))
scores <- rbind(spatial = metrics$spatial, eta = metrics$eta)
for (ch in seq_len(chains)) {
  scores <- rbind(scores, metric(chain_fields[, , ch], truth$spatial))
  rownames(scores)[nrow(scores)] <- paste0("spatial_chain", ch)
}
write.csv(scores, file.path(out_dir, "metrics.csv"))
range_summary <- data.frame(chain = seq_len(chains), truth = true_range,
  mean = metrics$range_mean, median = metrics$range_median,
  lower = metrics$range_quantiles[1, ], upper = metrics$range_quantiles[3, ],
  upper_boundary = colMeans(ro$idx_ls_output == length(grid)), transitions = metrics$range_transitions)
write.csv(range_summary, file.path(out_dir, "range-summary.csv"), row.names = FALSE)
range_distribution <- do.call(rbind, lapply(seq_len(chains), function(ch)
  data.frame(chain = ch, range = grid, count = as.vector(range_table[, ch]),
             proportion = as.vector(range_table[, ch]) / retained)))
write.csv(range_distribution, file.path(out_dir, "range-distribution.csv"), row.names = FALSE)
write.csv(data.frame(sigmabs_mean = metrics$sigmabs, tau_mean = metrics$tau,
                    intercept_bias = metrics$intercept_bias),
          file.path(out_dir, "parameter-summary.csv"), row.names = FALSE)
provenance$finished <- Sys.time()
provenance$elapsed <- elapsed
saveRDS(provenance, file.path(out_dir, "provenance.rds"))
print(metrics)
