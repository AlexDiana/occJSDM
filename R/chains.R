# Parallel MCMC chains: core budget, per-chain seeds, the single-chain
# function, the PSOCK driver and the combination helpers.
# Design: dev/superpowers/specs/2026-10-13-parallel-chains-design.md.

# Holds the prepared setup list while chains run, in the parent and in every
# worker. Package level, so that a worker can reach it by name after the task
# function arrives as a reference to the namespace.
chain_env <- new.env(parent = emptyenv())

# `available` is an argument so that the warning can be tested without
# depending on the operating system: parallel::detectCores() is documented to
# return NA when the count is unknown, and then no warning is possible.
resolve_ncores <- function(MCMCparams, available = parallel::detectCores()) {
  ncores <- get_param(MCMCparams, "ncores", 1)
  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores < 1 || ncores != floor(ncores)) {
    stop("MCMCparams$ncores must be a single whole number of at least 1", call. = FALSE)
  }
  ncores <- as.integer(ncores)
  if (!is.na(available) && ncores > available) {
    warning(sprintf("MCMCparams$ncores = %d exceeds the %d cores detected on this machine",
                    ncores, available), call. = FALSE)
  }
  ncores
}

# ncores is the total core budget: workers = min(nchain, ncores), each with
# floor(ncores / workers) RcppParallel threads. ncores = 1 is today's fit and
# leaves the session's thread setting alone ("ambient").
plan_chain_budget <- function(nchain, ncores) {
  workers <- as.integer(min(nchain, ncores))
  list(ncores = as.integer(ncores),
       workers = workers,
       threads_per_worker = max(1L, as.integer(floor(ncores / workers))),
       thread_setting = if (ncores == 1L) "ambient" else "budget")
}

ambient_rcpp_threads <- function() {
  env <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "")
  if (nzchar(env)) as.integer(env) else as.integer(RcppParallel::defaultNumThreads())
}

# Merge two Welford accumulators (running mean, M2 and mean likelihood per
# observation, each over n draws) exactly, with the standard pairwise formula.
# Chains are merged in chain order by merge_waic_parts(), so the result cannot
# depend on which worker produced which chain.
merge_waic_pair <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  n <- a$n + b$n
  delta <- b$mean_log - a$mean_log
  list(M2 = a$M2 + b$M2 + delta^2 * a$n * b$n / n,
       mean_log = a$mean_log + delta * b$n / n,
       mean_lik = (a$n * a$mean_lik + b$n * b$mean_lik) / n,
       n = n)
}

merge_waic_parts <- function(parts) {
  Reduce(merge_waic_pair, parts, accumulate = FALSE)
}

# Per-chain sums of the retained latent states, added in chain order. A name
# that is NULL or absent on either side comes out NULL, so a shape mismatch
# between chains cannot produce a silently partial sum. NULL entries are kept
# in the result (`[<-` with list(NULL), not `[[<-` with NULL), so the names
# stay the same however many parts are folded in.
add_latent_sums <- function(total, part) {
  if (is.null(part)) return(total)
  if (is.null(total)) return(part)
  nms <- union(names(total), names(part))
  out <- stats::setNames(vector("list", length(nms)), nms)
  for (nm in nms) {
    if (!is.null(total[[nm]]) && !is.null(part[[nm]])) out[[nm]] <- total[[nm]] + part[[nm]]
  }
  out
}
